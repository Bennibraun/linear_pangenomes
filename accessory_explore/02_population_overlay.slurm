#!/bin/bash
#SBATCH --job-name=acc_overlay
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=acc_overlay_%j.log
# ---------------------------------------------------------------------------
# Accessory exploration, step 2: overlay population presence/absence onto the
# candidate contigs + protein hits from step 1. Run AFTER 01_mask_and_search.slurm.
#   sbatch 02_population_overlay.slurm
# ---------------------------------------------------------------------------
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate accessory-explore

# ======================= CONFIG (edit these) ===============================
MERGED_VCF="results/variants/bcftools/augref/combined/merged.vcf.gz"
MANIFEST="reads_manifest.tsv"
OUTDIR="results/accessory_explore"
PRESENT_FRAC=0.5
# ===========================================================================

python 02_population_overlay.py \
  --vcf "$MERGED_VCF" \
  --diamond "$OUTDIR/candidate_diamond.tsv" \
  --repeat-frac "$OUTDIR/repeat_fraction.tsv" \
  --manifest "$MANIFEST" \
  --outdir "$OUTDIR" \
  --present-frac "$PRESENT_FRAC"

echo "[$(date)] overlay done. See:"
echo "  $OUTDIR/accessory_population_overlay.tsv     (all contigs, ranked)"
echo "  $OUTDIR/accessory_candidates_shortlist.tsv   (genic + structured)"
