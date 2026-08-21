#!/bin/bash
# Get a protein DB for the DIAMOND search (zebrafish proteome from Ensembl).
# Run on a node with internet, then set DIAMOND_DB in 01_mask_and_search.slurm.
set -euo pipefail

curl -O https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
gunzip -f Danio_rerio.GRCz11.pep.all.fa.gz
diamond makedb --in Danio_rerio.GRCz11.pep.all.fa -d danio_rerio -p "${SLURM_CPUS_PER_TASK:-8}"

echo "Done. Set in 01_mask_and_search.slurm:  DIAMOND_DB=\"$(pwd)/danio_rerio.dmnd\""
