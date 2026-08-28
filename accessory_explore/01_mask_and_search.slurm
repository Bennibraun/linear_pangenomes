#!/bin/bash
#SBATCH --job-name=acc_explore
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=acc_explore_%j.log
# Step 1: repeat-mask (dustmasker) the accessory contigs, then DIAMOND blastx
# the non-repetitive ones. Paths come from config.sh. sbatch this.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate accessory-explore
source config.sh

FASTA="$ACCESSORY_FASTA"
REPEAT_FRAC=0.50                     # contig is "repetitive" if >= this fraction masked
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUTDIR"
seqkit stats -a "$FASTA" | tee "$OUTDIR/fasta_stats.txt"

# Per-contig masked fraction. dustmasker acclist gives ">contig start end";
# sum masked bp per contig, join to contig lengths.
dustmasker -in "$FASTA" -outfmt acclist \
  | sed 's/^>//' | awk '{m[$1]+=$3-$2} END{for(c in m)print c"\t"m[c]}' \
  | sort > "$OUTDIR/masked_bp.tsv"
seqkit fx2tab -nl "$FASTA" | sort > "$OUTDIR/lengths.tsv"
join -a1 -e0 -o '1.1 1.2 2.2' "$OUTDIR/lengths.tsv" "$OUTDIR/masked_bp.tsv" \
  | awk 'BEGIN{OFS="\t";print "contig","length","masked_bp","masked_frac"}
         {print $1,$2,$3,($2?$3/$2:0)}' > "$OUTDIR/repeat_fraction.tsv"

# Split candidate (searchable) vs repetitive.
awk -F'\t' -v t="$REPEAT_FRAC" 'NR>1 && $4<t {print $1}' \
  "$OUTDIR/repeat_fraction.tsv" > "$OUTDIR/candidate_contigs.txt"
echo "candidates (< $REPEAT_FRAC masked): $(wc -l < "$OUTDIR/candidate_contigs.txt")"
seqkit grep -f "$OUTDIR/candidate_contigs.txt" "$FASTA" > "$OUTDIR/candidate.fasta"

# Protein search.
diamond blastx -q "$OUTDIR/candidate.fasta" -d "$DIAMOND_DB" \
  -o "$OUTDIR/candidate_diamond.tsv" -p "$THREADS" --very-sensitive -e 1e-5 -k 5 \
  -f 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore

echo "done -> $OUTDIR/{repeat_fraction,candidate_diamond}.tsv"
