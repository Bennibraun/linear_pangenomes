#!/bin/bash


# Usage: ./4.align_wgs.sh

eval "$(conda shell.bash hook)"
conda activate apis_ont


# This script will align the WGS samples listed in /Users/bebr1814/projects/bee_paper/data/public_wgs/samples.txt to these:
# - The linear pangenome from the A. mellifera data
augref=/Users/bebr1814/projects/bee_paper/data/sv_calls/augref/augmented_reference.fasta
# - The A. mellifera reference genome (Amel_HAv3.1) i.e. conspecific
conspec=/Users/bebr1814/projects/bee_paper/reference/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna
# - The A. cerana reference genome (Acera_v1.0) i.e. heterospecific
hetspec=/Users/bebr1814/projects/bee_paper/reference/AcerK_1.0/GCF_029169275.1_AcerK_1.0_genomic.fna
# - The cactus graph built from the A. mellifera data
cactus=/Users/bebr1814/projects/chuong_rotation/data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/cactus/output/apis_mellifera_10_assemblies.gbz

# Samples to align
sample_list=/Users/bebr1814/projects/bee_paper/data/public_wgs/samples.txt

outdir=/Users/bebr1814/projects/bee_paper/data/public_wgs/wgs_alignments

mkdir -p $outdir

mkdir -p logs

# index the references if not already indexed
if [ ! -f "${augref}.bwt" ]; then
    bwa index $augref
fi
if [ ! -f "${conspec}.bwt" ]; then
    bwa index $conspec
fi
if [ ! -f "${hetspec}.bwt" ]; then
    bwa index $hetspec
fi


# BWA Mem for fastas and VG Giraffe for the graph. Giraffe will use default parameters, while BWA Mem will allow more multimapping than usual due to the nature of the augmented reference.

# Loop through samples, submit one sbatch per sample that runs all 4 alignments sequentially
for fq1 in $(cat $sample_list); do
    fq2="${fq1/_1/_2}"
    sample=$(basename "$fq1" _1.fastq.gz)
    
    # Check if both R1 and R2 exist
    if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
        echo "Warning: Skipping $sample (R1 or R2 not found)"
        continue
    fi
    
    # Submit sbatch for this sampleg
    sbatch --job-name="${sample}_align" \
           --output="logs/${sample}_align.%j.out" \
           --error="logs/${sample}_align.%j.err" \
           --cpus-per-task=4 \
           --mem=20G \
           --time=11:00:00 \
           --partition=short \
           --wrap="eval \"\$(conda shell.bash hook)\"; conda activate apis_ont; module load samtools; \
                   mkdir -p $outdir/${sample}; \
                   # Augmented reference alignment
                   bwa mem -k 17 -O 5,5 -E 2,2 -Y -t 4 $augref $fq1 $fq2 | samtools sort -@ 4 -o $outdir/${sample}/${sample}.augref.bam; samtools index $outdir/${sample}/${sample}.augref.bam; \
                   # Conspecific reference alignment
                   bwa mem -k 17 -O 5,5 -E 2,2 -Y -t 4 $conspec $fq1 $fq2 | samtools sort -@ 4 -o $outdir/${sample}/${sample}.conspec.bam; samtools index $outdir/${sample}/${sample}.conspec.bam; \
                   # Heterospecific reference alignment
                   bwa mem -k 17 -O 5,5 -E 2,2 -Y -t 4 $hetspec $fq1 $fq2 | samtools sort -@ 4 -o $outdir/${sample}/${sample}.hetspec.bam; samtools index $outdir/${sample}/${sample}.hetspec.bam; \
                   # Graph alignment
                   vg giraffe -Z $cactus -t 4 -f $fq1 -f $fq2 -p --rescue-attempts 0 > $outdir/${sample}/${sample}.cactus.gam;"
done

