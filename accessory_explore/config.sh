# Paths for the accessory-explore analysis. Edit these, then run the numbered
# scripts. This analysis is SEPARATE from the Snakemake pipeline: it only READS
# the pipeline's results (via PIPELINE_RESULTS) and writes everything into its
# own OUTDIR. Nothing is written into the pipeline tree.

# Where the finished Snakemake pipeline wrote its results (the dir containing
# sv_calls/, variants/, etc.). Absolute path -- this analysis lives elsewhere.
PIPELINE_RESULTS="/Users/bebr1814/scratch/cavefish/snakemake/results"

# Where THIS analysis writes its outputs. Local to the analysis dir.
OUTDIR="output"

# Protein DB built by 00_get_proteome.sh (in this analysis dir).
DIAMOND_DB="danio_rerio.dmnd"

# Sample -> population manifest (same one the pipeline used).
MANIFEST="/Users/bebr1814/scratch/cavefish/snakemake/reads_manifest.tsv"

# --- derived input paths (usually no need to edit) ---
ACCESSORY_FASTA="$PIPELINE_RESULTS/sv_calls/augref/extracted_flanked_sv_seqs.dedup.fasta"
MERGED_VCF="$PIPELINE_RESULTS/variants/bcftools/augref/combined/merged.vcf.gz"
