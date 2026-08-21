#!/bin/bash
# ---------------------------------------------------------------------------
# Download a protein DB for the DIAMOND search (step 1) and build the .dmnd.
#
# Default: zebrafish (Danio rerio) UniProt reference proteome UP000000437
# (~68,855 sequences as of 2026_02). Run on a LOGIN/TRANSFER node (needs
# internet; compute nodes are often air-gapped). Then set
# DIAMOND_DB="<here>/danio_rerio.dmnd" in 01_mask_and_search.slurm.
#
#   bash 00_get_proteome.sh                 # zebrafish, full proteome
#   REVIEWED=true bash 00_get_proteome.sh   # curated Swiss-Prot only (smaller)
#
# Notes:
#  - Uses the paginated `search` endpoint (follows Link: rel="next"); the
#    `stream` endpoint drops large downloads mid-transfer -- don't use it.
#  - Deliberately NOT `set -e` and NOT `curl -s`: earlier this failed silently.
#    Now every page prints its HTTP status, and a failed page is retried and
#    then reported loudly instead of killing the script without a message.
#  - Resumable: re-running continues from where it left off (keeps a cursor
#    file). Delete ${NAME}.faa and ${NAME}.cursor to start fresh.
# ---------------------------------------------------------------------------

PROTEOME="${1:-UP000000437}"
NAME="${2:-danio_rerio}"
PAGE_SIZE=500
OUT="${NAME}.faa"
CURSOR_FILE="${NAME}.cursor"
REVIEWED="${REVIEWED:-false}"

if [ "$REVIEWED" = "true" ]; then
    query="%28proteome%3A${PROTEOME}%29%20AND%20%28reviewed%3Atrue%29"
else
    query="%28proteome%3A${PROTEOME}%29"
fi
base="https://rest.uniprot.org/uniprotkb/search?query=${query}&format=fasta&size=${PAGE_SIZE}"

# Resume support: if a cursor was saved, start from it and append; else fresh.
if [ -f "$CURSOR_FILE" ] && [ -s "$OUT" ]; then
    url=$(cat "$CURSOR_FILE")
    echo "[$(date)] resuming from saved cursor ($(grep -c '^>' "$OUT") seqs so far)"
else
    url="$base"
    : > "$OUT"
    rm -f "$CURSOR_FILE"
    echo "[$(date)] downloading proteome $PROTEOME -> $OUT (reviewed=$REVIEWED)"
fi

page=0
while [ -n "$url" ]; do
    page=$((page + 1))
    hdrs=$(mktemp)
    body=$(mktemp)

    # -f: HTTP errors -> nonzero exit. NO -s, so curl prints network errors.
    # --retry covers transient drops; long timeouts for slow login-node links.
    http=$(curl -f -w '%{http_code}' \
                --retry 6 --retry-delay 5 --retry-all-errors \
                --connect-timeout 30 --max-time 600 \
                -D "$hdrs" -o "$body" "$url")
    rc=$?

    if [ $rc -ne 0 ]; then
        echo ""
        echo "ERROR: curl failed on page $page (exit $rc, http=$http)." >&2
        echo "  URL: $url" >&2
        case $rc in
          6)  echo "  -> DNS failure. On an HPC node you may need a proxy:" >&2
              echo "       export https_proxy=http://<your-cluster-proxy>:<port>" >&2 ;;
          7)  echo "  -> connection refused/blocked. Login node may lack outbound" >&2
              echo "     internet; try a dedicated transfer/data node, or set https_proxy." >&2 ;;
          28) echo "  -> timed out. Network is very slow; re-run to resume from cursor." >&2 ;;
          *)  echo "  -> see 'man curl' exit code $rc." >&2 ;;
        esac
        echo "Progress saved: $(grep -c '^>' "$OUT" 2>/dev/null || echo 0) seqs in $OUT." >&2
        echo "Re-run the script to resume." >&2
        rm -f "$hdrs" "$body"
        exit $rc
    fi

    # Guard: HTTP 200 but empty body (the silent-nothing case you hit).
    nbytes=$(wc -c < "$body")
    if [ "$nbytes" -eq 0 ]; then
        echo "WARNING: page $page returned 0 bytes (http=$http). Stopping." >&2
        rm -f "$hdrs" "$body"
        break
    fi

    cat "$body" >> "$OUT"
    n=$(grep -c '^>' "$OUT")
    total=$(grep -i '^x-total-results:' "$hdrs" | tr -d '\r' | awk '{print $2}')
    echo "  page $page  http=$http  cumulative=$n${total:+ / $total}"

    # follow Link: <...>; rel="next"; save it so we can resume.
    url=$(grep -i '^link:' "$hdrs" | sed -nE 's/.*<([^>]+)>; *rel="next".*/\1/p' | tr -d '\r')
    if [ -n "$url" ]; then echo "$url" > "$CURSOR_FILE"; else rm -f "$CURSOR_FILE"; fi
    rm -f "$hdrs" "$body"
done

nseq=$(grep -c '^>' "$OUT")
echo "[$(date)] done: $nseq sequences in $OUT"
if [ "$nseq" -eq 0 ]; then
    echo "ERROR: 0 sequences downloaded. See errors above." >&2
    exit 1
fi

if command -v diamond >/dev/null 2>&1; then
    echo "[$(date)] building DIAMOND db -> ${NAME}.dmnd"
    diamond makedb --in "$OUT" -d "$NAME" -p "${SLURM_CPUS_PER_TASK:-8}"
    echo "Set in 01_mask_and_search.slurm:  DIAMOND_DB=\"$(pwd)/${NAME}.dmnd\""
else
    echo "diamond not on PATH; skipped makedb."
    echo "Either 'conda activate accessory-explore' and rerun, or set"
    echo "  PROTEIN_FASTA=\"$(pwd)/${OUT}\"  in 01_mask_and_search.slurm."
fi
