#!/bin/bash
# ---------------------------------------------------------------------------
# Download a protein DB for the DIAMOND search (step 1) and build the .dmnd.
#
# Default: zebrafish (Danio rerio) UniProt reference proteome UP000000437.
# Run on a LOGIN/TRANSFER node (needs internet; compute nodes are often
# air-gapped). Then set DIAMOND_DB="<here>/danio_rerio.dmnd" in
# 01_mask_and_search.slurm.
#
#   bash 00_get_proteome.sh                 # zebrafish, full proteome
#   bash 00_get_proteome.sh UP000000437 danio_rerio   # explicit
#
# Why paginated: UniProt's `stream` endpoint drops large downloads mid-transfer
# (HTTP/2 stream errors). The `search` endpoint paginates via a Link: rel="next"
# header and is reliable. This script follows that cursor, appending pages, and
# is resumable-ish (delete the .faa to restart cleanly).
# ---------------------------------------------------------------------------
set -euo pipefail

PROTEOME="${1:-UP000000437}"     # zebrafish reference proteome
NAME="${2:-danio_rerio}"
PAGE_SIZE=500
OUT="${NAME}.faa"

# reviewed-only? set REVIEWED=true to get just curated Swiss-Prot (smaller,
# faster, fewer genes). Default false = full proteome (recommended for us).
REVIEWED="${REVIEWED:-false}"

if [ "$REVIEWED" = "true" ]; then
    query="%28proteome%3A${PROTEOME}%29%20AND%20%28reviewed%3Atrue%29"
else
    query="%28proteome%3A${PROTEOME}%29"
fi

url="https://rest.uniprot.org/uniprotkb/search?query=${query}&format=fasta&size=${PAGE_SIZE}"

echo "[$(date)] downloading proteome $PROTEOME -> $OUT (reviewed=$REVIEWED)"
: > "$OUT"                       # truncate/create
page=0
while [ -n "$url" ]; do
    page=$((page + 1))
    hdrs=$(mktemp)
    # --retry handles transient drops; -f fails on HTTP errors so we notice.
    curl -sf --retry 5 --retry-delay 3 --max-time 300 \
         -D "$hdrs" "$url" >> "$OUT"
    n=$(grep -c '^>' "$OUT" || true)
    echo "  page $page  (cumulative sequences: $n)"
    # follow Link: <...>; rel="next"
    url=$(grep -i '^link:' "$hdrs" | sed -nE 's/.*<([^>]+)>; *rel="next".*/\1/p' || true)
    rm -f "$hdrs"
done

nseq=$(grep -c '^>' "$OUT")
echo "[$(date)] done: $nseq sequences in $OUT"

# Build the DIAMOND db if diamond is available (e.g. after `conda activate
# accessory-explore`). Otherwise just leave the FASTA and point PROTEIN_FASTA
# at it in step 1.
if command -v diamond >/dev/null 2>&1; then
    echo "[$(date)] building DIAMOND db -> ${NAME}.dmnd"
    diamond makedb --in "$OUT" -d "$NAME" -p "${SLURM_CPUS_PER_TASK:-8}"
    echo "Set in 01_mask_and_search.slurm:  DIAMOND_DB=\"$(pwd)/${NAME}.dmnd\""
else
    echo "diamond not on PATH; skipped makedb."
    echo "Either 'conda activate accessory-explore' and rerun, or set"
    echo "  PROTEIN_FASTA=\"$(pwd)/${OUT}\"  in 01_mask_and_search.slurm."
fi
