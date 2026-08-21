#!/bin/bash
#SBATCH --job-name=acc_overlay
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=acc_overlay_%j.log
# Step 2: overlay population presence/absence onto candidate contigs + hits.
# Run after 01_mask_and_search.slurm.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate accessory-explore

OUTDIR="results/accessory_explore"
python 02_population_overlay.py \
  results/variants/bcftools/augref/combined/merged.vcf.gz \
  "$OUTDIR/candidate_diamond.tsv" \
  "$OUTDIR/repeat_fraction.tsv" \
  reads_manifest.tsv \
  "$OUTDIR"
