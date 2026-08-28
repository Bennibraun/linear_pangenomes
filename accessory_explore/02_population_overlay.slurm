#!/bin/bash
#SBATCH --job-name=acc_overlay
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=acc_overlay_%j.log
# Step 2: overlay population presence/absence onto candidate contigs + hits.
# Paths from config.sh. Run after 01_mask_and_search.slurm.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate accessory-explore
source config.sh

python 02_population_overlay.py \
  "$MERGED_VCF" \
  "$OUTDIR/candidate_diamond.tsv" \
  "$OUTDIR/repeat_fraction.tsv" \
  "$MANIFEST" \
  "$OUTDIR"
