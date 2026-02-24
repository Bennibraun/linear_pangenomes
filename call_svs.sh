#!/bin/bash



################################################################################
# SV Calling Pipeline for Apis mellifera ONT Data
# 
# This pipeline performs:
# 1. Read-based SV calling (Sniffles2, cuteSV)
# 2. Assembly-based SV calling (minimap2/paftools, Assemblytics, SyRI)
# 3. Per-sample merging of callsets
# 4. Pan-sample catalog generation
#
# Requirements:
# - minimap2
# - samtools
# - sniffles (v2.x)
# - cuteSV
# - k8 (javascript shell for paftools.js)
# - paftools.js (from minimap2/misc/)
# - MUMmer (nucmer, delta-filter, show-coords)
# - Assemblytics
# - SyRI
# - SURVIVOR
# - Jasmine
# - bcftools
################################################################################

# Usage: ./2.call_SVs_and_filter.sh

if [[ -z "${SNAKEMAKE_CONDA_PREFIX:-}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate apis_ont
fi

# On Fiji, bcftools and samtools don't have the necessary libraries linked.
# Instead load the available modules for these.
if command -v module >/dev/null 2>&1; then
    module load samtools bcftools
fi


# Configuration
THREADS=${THREADS:-8}
REFERENCE=${REFERENCE:-"/Users/bebr1814/scratch/chuong_data/reference/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna"}
SURVIVOR_EXEC=${SURVIVOR_EXEC:-"SURVIVOR"}
MIN_SV_SIZE=${MIN_SV_SIZE:-50}
MIN_READ_SUPPORT=${MIN_READ_SUPPORT:-3}
BREAKPOINT_SLOP=${BREAKPOINT_SLOP:-1000}  # For SURVIVOR merge
JASMINE_SLOP=${JASMINE_SLOP:-500}      # For Jasmine merge

# Directories
READS_DIR=${READS_DIR:-"/Users/bebr1814/scratch/chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/fastq"}        # Directory containing ONT fastq files
ASSEMBLY_DIR=${ASSEMBLY_DIR:-"/Users/bebr1814/scratch/chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/bee_paper/flye_assembly"}    # Directory containing Flye assembly fastas
OUTPUT_DIR=${OUTPUT_DIR:-"/Users/bebr1814/scratch/chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/bee_paper/sv_calls"}
READ_ALIGN_DIR="${OUTPUT_DIR}/read_alignments"
ASM_ALIGN_DIR="${OUTPUT_DIR}/assembly_alignments"
READ_CALLS_DIR="${OUTPUT_DIR}/read_sv_calls"
ASM_CALLS_DIR="${OUTPUT_DIR}/assembly_sv_calls"
MERGED_CALLS_DIR="${OUTPUT_DIR}/merged_per_sample"
CATALOG_DIR="${OUTPUT_DIR}/pan_sample_catalog"

cd ${OUTPUT_DIR}

# Create output directories
mkdir -p ${READ_ALIGN_DIR} ${ASM_ALIGN_DIR} ${READ_CALLS_DIR} ${ASM_CALLS_DIR}
mkdir -p ${MERGED_CALLS_DIR} ${CATALOG_DIR}

if [[ -n "${SAMPLES_FILE:-}" && -f "${SAMPLES_FILE}" ]]; then
    mapfile -t SAMPLES < "${SAMPLES_FILE}"
else
    SAMPLES=(
        florida1.barcode10
        florida2.barcode1
        florida3.barcode2
        florida4.barcode3
        thailand1.barcode12
        thailand2.barcode4
        thailand4.barcode5
        tokyo1.barcode11
        tokyo2.barcode6
        tokyo3.barcode7
    )
fi

# SAMPLES=(
#     tokyo3.barcode7
# )

################################################################################
# STEP 1: READ-BASED SV CALLING
################################################################################
echo "========================================"
echo "STEP 1: Read-based SV calling"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"
    
    # Define input files
    READS="${READS_DIR}/${SAMPLE}.fastq"  # Adjust extension as needed
    
    # Output files
    SAM="${READ_ALIGN_DIR}/${SAMPLE}.sam"
    BAM="${READ_ALIGN_DIR}/${SAMPLE}.bam"
    SORTED_BAM="${READ_ALIGN_DIR}/${SAMPLE}.sorted.bam"
    
    # Step 1a: Align ONT reads to reference with minimap2
    echo "  Aligning reads to reference with minimap2..."
    minimap2 -ax map-ont \
        -t ${THREADS} \
        --MD \
        ${REFERENCE} \
        ${READS} > ${SAM}
    
    # Convert to BAM, sort, and index
    echo "  Converting to BAM and sorting..."
    samtools view -@ ${THREADS} -bS ${SAM} > ${BAM}
    samtools sort -@ ${THREADS} -o ${SORTED_BAM} ${BAM}
    samtools index -@ ${THREADS} ${SORTED_BAM}
    
    # Clean up intermediate files
    rm ${SAM} ${BAM}
    
    # Step 1b: Call SVs with Sniffles2
    echo "  Calling SVs with Sniffles2..."
    sniffles \
        --input ${SORTED_BAM} \
        --vcf ${READ_CALLS_DIR}/${SAMPLE}.sniffles2.vcf \
        --reference ${REFERENCE} \
        --threads ${THREADS} \
        --minsvlen ${MIN_SV_SIZE} \
        --minsupport ${MIN_READ_SUPPORT}
    
    # Step 1c: Call SVs with cuteSV
    echo "  Calling SVs with cuteSV..."
    # Create temporary directory for cuteSV
    mkdir -p ${READ_CALLS_DIR}/${SAMPLE}_cutesv_temp
    
    cuteSV \
        ${SORTED_BAM} \
        ${REFERENCE} \
        ${READ_CALLS_DIR}/${SAMPLE}.cutesv.vcf \
        ${READ_CALLS_DIR}/${SAMPLE}_cutesv_temp \
        --threads ${THREADS} \
        --min_size ${MIN_SV_SIZE} \
        --min_support ${MIN_READ_SUPPORT} \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --genotype
    
    # Clean up cuteSV temp directory
    rm -rf ${READ_CALLS_DIR}/${SAMPLE}_cutesv_temp
    
    echo "  Completed read-based SV calling for ${SAMPLE}"
done

################################################################################
# STEP 2: ASSEMBLY-BASED SV CALLING
################################################################################
echo "========================================"
echo "STEP 2: Assembly-based SV calling"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"
    
    # Define input files
    ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.assembly.fasta"  # Adjust extension as needed
    
    # Output files
    ASM_PAF="${ASM_ALIGN_DIR}/${SAMPLE}.paf"
    ASM_SAM="${ASM_ALIGN_DIR}/${SAMPLE}.sam"
    ASM_BAM="${ASM_ALIGN_DIR}/${SAMPLE}.bam"
    ASM_SORTED_BAM="${ASM_ALIGN_DIR}/${SAMPLE}.sorted.bam"
    
    # Step 2a: Align assembly to reference with minimap2
    echo "  Aligning assembly to reference..."
    # Generate both PAF (for paftools) and SAM (for SyRI)
    minimap2 -cx asm5 \
        -t ${THREADS} \
        --cs \
        ${REFERENCE} \
        ${ASSEMBLY} > ${ASM_PAF}
    
    # minimap2 -ax asm5 \
    #     -t ${THREADS} \
    #     --eqx \
    #     ${REFERENCE} \
    #     ${ASSEMBLY} > ${ASM_SAM}
    
    # Convert SAM to sorted BAM for SyRI
    # samtools view -@ ${THREADS} -bS ${ASM_SAM} > ${ASM_BAM}
    # samtools sort -@ ${THREADS} -o ${ASM_SORTED_BAM} ${ASM_BAM}
    # samtools index -@ ${THREADS} ${ASM_SORTED_BAM}
    # rm ${ASM_BAM}
    
    # Step 2b: Call SVs with paftools
    echo "  Calling SVs with paftools..."
    sort -k6,6 -k8,8n ${ASM_PAF} > ${ASM_ALIGN_DIR}/${SAMPLE}.sorted.paf
    paftools.js call \
        -f ${REFERENCE} \
        -L ${MIN_SV_SIZE} \
        ${ASM_ALIGN_DIR}/${SAMPLE}.sorted.paf > ${ASM_CALLS_DIR}/${SAMPLE}.paftools.vcf

    # # Step 2c: Call SVs with SyRI
    # echo "  Calling SVs with SyRI..."
    # syri \
    #     -c ${ASM_SAM} \
    #     -r ${REFERENCE} \
    #     -q ${ASSEMBLY} \
    #     -F S \
    #     --prefix ${SAMPLE}_syri \
    #     --dir ${ASM_CALLS_DIR}

    
    # # SyRI outputs multiple files; the VCF is syri.vcf
    # mv ${ASM_CALLS_DIR}/syri.vcf ${ASM_CALLS_DIR}/${SAMPLE}.syri.vcf || true
    
    # Clean up intermediate files
    # rm -f ${ASM_SAM}
    
    echo "  Completed assembly-based SV calling for ${SAMPLE}"
done

################################################################################
# STEP 3: MERGE CALLSETS PER SAMPLE
################################################################################
echo "========================================"
echo "STEP 3: Merging callsets per sample"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Merging callsets for sample: ${SAMPLE}"
    
    # Step 3a: Merge read-based calls (Sniffles2 + cuteSV)
    echo "  Merging read-based calls..."
    
    # Create file list for SURVIVOR
    READ_VCF_LIST="${MERGED_CALLS_DIR}/${SAMPLE}_read_vcfs.txt"
    echo "${READ_CALLS_DIR}/${SAMPLE}.sniffles2.vcf" > ${READ_VCF_LIST}
    echo "${READ_CALLS_DIR}/${SAMPLE}.cutesv.vcf" >> ${READ_VCF_LIST}
    
    # Merge with SURVIVOR
    # Parameters: max_distance min_support type_aware strand_aware size_estimate min_size output
    $SURVIVOR_EXEC merge \
        ${READ_VCF_LIST} \
        ${BREAKPOINT_SLOP} \
        1 \
        1 \
        1 \
        0 \
        ${MIN_SV_SIZE} \
        ${MERGED_CALLS_DIR}/${SAMPLE}.read_based.merged.vcf
    
    # Step 3b: Merge assembly-based calls
    echo "  Merging assembly-based calls..."
    
    # Note: Assemblytics outputs BED format, need to convert or handle separately
    # For this pipeline, we'll merge only VCF-format calls from paftools and SyRI
    ASM_VCF_LIST="${MERGED_CALLS_DIR}/${SAMPLE}_asm_vcfs.txt"
    echo "${ASM_CALLS_DIR}/${SAMPLE}.paftools.vcf" > ${ASM_VCF_LIST}
    
    $SURVIVOR_EXEC merge \
        ${ASM_VCF_LIST} \
        ${BREAKPOINT_SLOP} \
        1 \
        1 \
        1 \
        0 \
        ${MIN_SV_SIZE} \
        ${MERGED_CALLS_DIR}/${SAMPLE}.assembly_based.merged.vcf
    
    # Step 3c: Merge read-based and assembly-based callsets
    echo "  Merging read and assembly callsets..."
    
    FINAL_VCF_LIST="${MERGED_CALLS_DIR}/${SAMPLE}_all_vcfs.txt"
    echo "${MERGED_CALLS_DIR}/${SAMPLE}.read_based.merged.vcf" > ${FINAL_VCF_LIST}
    echo "${MERGED_CALLS_DIR}/${SAMPLE}.assembly_based.merged.vcf" >> ${FINAL_VCF_LIST}
    
    # Require support from both read-based OR assembly-based (min_support=1)
    # OR require support from both (min_support=2) for higher confidence
    # Using min_support=1 for comprehensive catalog
    $SURVIVOR_EXEC merge \
        ${FINAL_VCF_LIST} \
        ${BREAKPOINT_SLOP} \
        1 \
        1 \
        1 \
        0 \
        ${MIN_SV_SIZE} \
        ${MERGED_CALLS_DIR}/${SAMPLE}.high_confidence.vcf
    
    echo "  Completed merging for ${SAMPLE}"
done

################################################################################
# STEP 4: GENERATE PAN-SAMPLE SV CATALOG
################################################################################
echo "========================================"
echo "STEP 4: Generating pan-sample SV catalog"
echo "========================================"

# Step 4a: Create list of all sample VCFs
SAMPLE_VCF_LIST="${CATALOG_DIR}/all_samples.txt"
> ${SAMPLE_VCF_LIST}  # Clear file
for SAMPLE in "${SAMPLES[@]}"; do
    echo "${MERGED_CALLS_DIR}/${SAMPLE}.high_confidence.vcf" >> ${SAMPLE_VCF_LIST}
done

# Step 4b: Merge with SURVIVOR
echo "  Merging all samples with SURVIVOR..."
$SURVIVOR_EXEC merge \
    ${SAMPLE_VCF_LIST} \
    ${BREAKPOINT_SLOP} \
    1 \
    1 \
    1 \
    1 \
    ${MIN_SV_SIZE} \
    ${CATALOG_DIR}/pan_sample_catalog.survivor.vcf

# Step 4c: Merge with Jasmine (alternative/complementary approach)
echo "  Merging all samples with Jasmine..."
jasmine \
    file_list=${SAMPLE_VCF_LIST} \
    out_file=${CATALOG_DIR}/pan_sample_catalog.jasmine.vcf \
    max_dist=${JASMINE_SLOP} \
    --output_genotypes \
    --normalize_type

################################################################################
# STEP 5: GENERATE SUMMARY STATISTICS
################################################################################
echo "========================================"
echo "STEP 5: Generating summary statistics"
echo "========================================"

# Generate stats for the final catalogs
echo "SURVIVOR catalog statistics:" > ${CATALOG_DIR}/catalog_stats.txt
$SURVIVOR_EXEC stats ${CATALOG_DIR}/pan_sample_catalog.survivor.vcf >> ${CATALOG_DIR}/catalog_stats.txt

echo "" >> ${CATALOG_DIR}/catalog_stats.txt
echo "Jasmine catalog statistics:" >> ${CATALOG_DIR}/catalog_stats.txt
bcftools stats ${CATALOG_DIR}/pan_sample_catalog.jasmine.vcf >> ${CATALOG_DIR}/catalog_stats.txt

# Count SVs per sample
echo "" >> ${CATALOG_DIR}/catalog_stats.txt
echo "SV counts per sample (high-confidence callset):" >> ${CATALOG_DIR}/catalog_stats.txt
for SAMPLE in "${SAMPLES[@]}"; do
    COUNT=$(grep -v "^#" ${MERGED_CALLS_DIR}/${SAMPLE}.high_confidence.vcf | wc -l)
    echo "${SAMPLE}: ${COUNT}" >> ${CATALOG_DIR}/catalog_stats.txt
done

# Extract support vector from SURVIVOR output to analyze SV sharing
echo "" >> ${CATALOG_DIR}/catalog_stats.txt
echo "Extracting SV support patterns..." >> ${CATALOG_DIR}/catalog_stats.txt
grep -oP 'SUPP_VEC=\K[^,;]+' ${CATALOG_DIR}/pan_sample_catalog.survivor.vcf | \
    sed 's/./& /g' > ${CATALOG_DIR}/sv_support_matrix.txt


################################################################################
# PIPELINE COMPLETE
################################################################################
echo "========================================"
echo "Pipeline complete!"
echo "========================================"
echo "Output directories:"
echo "  Read alignments: ${READ_ALIGN_DIR}"
echo "  Assembly alignments: ${ASM_ALIGN_DIR}"
echo "  Read-based SV calls: ${READ_CALLS_DIR}"
echo "  Assembly-based SV calls: ${ASM_CALLS_DIR}"
echo "  Per-sample merged calls: ${MERGED_CALLS_DIR}"
echo "  Pan-sample catalogs: ${CATALOG_DIR}"
echo ""
echo "Key output files:"
echo "  SURVIVOR catalog: ${CATALOG_DIR}/pan_sample_catalog.survivor.vcf"
echo "  Jasmine catalog: ${CATALOG_DIR}/pan_sample_catalog.jasmine.vcf"
echo "  Summary statistics: ${CATALOG_DIR}/catalog_stats.txt"
echo "  Support matrix: ${CATALOG_DIR}/sv_support_matrix.txt"


## Post-processing

# Extract the sequences of SV insertions - including a 200bp flank on either end
mkdir -p $OUTPUT_DIR/augref
python extract_sv_sequences.py -v ${CATALOG_DIR}/pan_sample_catalog.survivor.vcf -r $REFERENCE --flank ${FLANK} -o $OUTPUT_DIR/augref/extracted_flanked_sv_seqs.fasta

# Deduplicate sequences to avoid redundancy in the augmented reference
cd-hit-est -i $OUTPUT_DIR/augref/extracted_flanked_sv_seqs.fasta -o $OUTPUT_DIR/augref/extracted_flanked_sv_seqs.dedup.fasta -c 0.95 -n 10

# Generate augmented reference
cat $REFERENCE $OUTPUT_DIR/augref/extracted_flanked_sv_seqs.dedup.fasta > $OUTPUT_DIR/augref/augmented_reference.fasta

