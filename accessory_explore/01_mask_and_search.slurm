#!/bin/bash
#SBATCH --job-name=acc_explore
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=acc_explore_%j.log
# ---------------------------------------------------------------------------
# Accessory-sequence exploration, step 1: repeat-mask + protein search.
#
# What it does, in order:
#   1. dustmasker on the accessory FASTA -> per-contig low-complexity/repeat
#      fraction (fast first-pass repeat filter; no repeat MODELLING).
#   2. Split contigs into "repetitive" (>= REPEAT_FRAC masked) and "candidate"
#      (< REPEAT_FRAC masked, i.e. mostly-unique sequence worth searching).
#   3. DIAMOND blastx the CANDIDATE contigs against a protein DB to find coding
#      potential / best hits (TEs will hit transposase/gag-pol and get flagged
#      in the analysis step).
#
# NOT a Snakemake rule -- run by hand for the paper. Edit the CONFIG block, then:
#   sbatch 01_mask_and_search.slurm
# ---------------------------------------------------------------------------
set -euo pipefail

# ======================= CONFIG (edit these) ===============================
# Activate the standalone env (adjust to your conda setup).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate accessory-explore

# Accessory-only FASTA: the deduped flanked SV_* insertion sequences (NOT the
# full augmented_reference, which also contains conspec).
ACCESSORY_FASTA="results/sv_calls/augref/extracted_flanked_sv_seqs.dedup.fasta"

# Output directory for all exploration products.
OUTDIR="results/accessory_explore"

# Protein DB for DIAMOND. Provide ONE of:
#   (a) a prebuilt DIAMOND db (.dmnd)  -> set DIAMOND_DB
#   (b) a protein FASTA (.faa/.fasta)  -> set PROTEIN_FASTA, script builds the db
# Recommended: a teleost/vertebrate proteome, e.g. zebrafish (Danio rerio)
# UniProt reference proteome, or a broad set like UniRef90. Bigger DB = slower.
DIAMOND_DB=""                      # e.g. /path/to/danio_rerio.dmnd
PROTEIN_FASTA=""                   # e.g. /path/to/danio_rerio.faa  (built if DIAMOND_DB empty)

# Contig is "repetitive" if this fraction or more of its bases are dustmasked.
REPEAT_FRAC=0.50

# DIAMOND sensitivity / thresholds.
THREADS=${SLURM_CPUS_PER_TASK:-16}
DIAMOND_EVALUE=1e-5
DIAMOND_SENS="--very-sensitive"    # --sensitive is faster; --very-sensitive catches diverged hits
# ===========================================================================

mkdir -p "$OUTDIR"

echo "[$(date)] accessory FASTA: $ACCESSORY_FASTA"
seqkit stats -a "$ACCESSORY_FASTA" | tee "$OUTDIR/accessory_fasta_stats.txt"

# Per-contig length table (contig, length). Used everywhere downstream.
seqkit fx2tab -nl "$ACCESSORY_FASTA" \
  | awk 'BEGIN{OFS="\t"; print "contig","length"} {print $1,$2}' \
  > "$OUTDIR/contig_lengths.tsv"

# --- 1. dustmasker -> masked intervals (acclist format: >contig start end) --
echo "[$(date)] dustmasker ..."
dustmasker -in "$ACCESSORY_FASTA" -outfmt acclist \
  > "$OUTDIR/dust.acclist" 2> "$OUTDIR/dustmasker.log" || {
    echo "dustmasker failed; see $OUTDIR/dustmasker.log" >&2; exit 1; }

# acclist lines: ">contigname\tstart\tend" (0-based, half-open). Sum masked bp
# per contig, then divide by contig length to get masked fraction.
python3 - "$OUTDIR/dust.acclist" "$OUTDIR/contig_lengths.tsv" "$OUTDIR/repeat_fraction.tsv" <<'PY'
import sys
acclist, lentsv, out = sys.argv[1:4]
masked = {}
with open(acclist) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        name = parts[0].lstrip(">")
        try:
            s, e = int(parts[1]), int(parts[2])
        except ValueError:
            continue
        masked[name] = masked.get(name, 0) + max(0, e - s)
lengths = {}
with open(lentsv) as fh:
    next(fh, None)
    for line in fh:
        c, l = line.rstrip("\n").split("\t")
        lengths[c] = int(l)
with open(out, "w") as fh:
    fh.write("contig\tlength\tmasked_bp\tmasked_frac\n")
    for c, L in lengths.items():
        m = masked.get(c, 0)
        fh.write(f"{c}\t{L}\t{m}\t{m/L if L else 0:.4f}\n")
print(f"wrote {out} for {len(lengths)} contigs")
PY

# --- 2. split candidate vs repetitive by masked fraction --------------------
echo "[$(date)] splitting at masked_frac >= $REPEAT_FRAC ..."
awk -F'\t' -v thr="$REPEAT_FRAC" 'NR>1 && $4>=thr {print $1}' \
  "$OUTDIR/repeat_fraction.tsv" > "$OUTDIR/repetitive_contigs.txt"
awk -F'\t' -v thr="$REPEAT_FRAC" 'NR>1 && $4<thr  {print $1}' \
  "$OUTDIR/repeat_fraction.tsv" > "$OUTDIR/candidate_contigs.txt"
n_rep=$(wc -l < "$OUTDIR/repetitive_contigs.txt")
n_cand=$(wc -l < "$OUTDIR/candidate_contigs.txt")
echo "  repetitive (>= $REPEAT_FRAC masked): $n_rep"
echo "  candidate  (<  $REPEAT_FRAC masked): $n_cand"

seqkit grep -f "$OUTDIR/candidate_contigs.txt" "$ACCESSORY_FASTA" \
  > "$OUTDIR/candidate_contigs.fasta"

# --- 3. DIAMOND blastx of candidate contigs ---------------------------------
if [ -z "$DIAMOND_DB" ]; then
    if [ -z "$PROTEIN_FASTA" ]; then
        echo "ERROR: set DIAMOND_DB or PROTEIN_FASTA in the CONFIG block." >&2
        echo "Masking is done ($OUTDIR/repeat_fraction.tsv); skipping DIAMOND." >&2
        exit 2
    fi
    echo "[$(date)] building DIAMOND db from $PROTEIN_FASTA ..."
    DIAMOND_DB="$OUTDIR/proteins.dmnd"
    diamond makedb --in "$PROTEIN_FASTA" -d "$OUTDIR/proteins" -p "$THREADS"
fi

echo "[$(date)] DIAMOND blastx ($DIAMOND_SENS, e<=$DIAMOND_EVALUE) ..."
# outfmt 6 with the fields the analysis script expects.
diamond blastx \
  -q "$OUTDIR/candidate_contigs.fasta" \
  -d "$DIAMOND_DB" \
  -o "$OUTDIR/candidate_diamond.tsv" \
  -p "$THREADS" \
  $DIAMOND_SENS \
  -e "$DIAMOND_EVALUE" \
  -k 5 \
  -f 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore

echo "[$(date)] done. Key outputs in $OUTDIR:"
echo "  repeat_fraction.tsv       per-contig masked fraction"
echo "  candidate_contigs.fasta   non-repetitive contigs searched"
echo "  candidate_diamond.tsv     protein hits for candidate contigs"
echo "Next: run 02_population_overlay.py to cross hits with population presence/absence."
