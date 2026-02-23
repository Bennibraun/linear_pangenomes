#!/bin/bash

# This script runs simple QC on input nanopore data and creates a flye assembly.
# It then performs several QC steps on the resulting assembly.

# Usage: ./1.assemblies.sh <input_fastq>
# e.g. ./1.assemblies.sh /Users/bebr1814/scratch/chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/fastq/tokyo3.barcode7.fastq


eval "$(conda shell.bash hook)"
conda activate nanopore_asm

cd /Users/bebr1814/projects/bee_paper/data

# ---------------------
# Config
# ---------------------
GENOME_SIZE="225m"
THREADS=1
LINEAGE="hymenoptera_odb10"


# ---------------------
# Input handling
# ---------------------
FASTQ="$1"
BASENAME=$(basename "$FASTQ")
PREFIX="${BASENAME%.fastq}"
PREFIX="${PREFIX%.fq}"

WORKDIR=$(pwd)

echo "Input FASTQ: $FASTQ"
echo "Prefix: $PREFIX"
echo "Working directory: $WORKDIR"
echo "Threads: $THREADS"
date
hostname

########################
# NanoStat
########################
echo "[$(date +"%b %d %H:%M:%S")] Running NanoStat..."

NANOSTAT_OUT="nanostat"
mkdir -p "$NANOSTAT_OUT"

NanoStat \
  --fastq "$FASTQ" \
  --name "${PREFIX}_nanostat.txt" \
  --outdir "$NANOSTAT_OUT" \
  --threads "$THREADS"

echo "[$(date +"%b %d %H:%M:%S")] NanoStat done."

########################
# NanoPlot
########################
echo "[$(date +"%b %d %H:%M:%S")] Running NanoPlot..."

NANOPLOT_OUT="nanoplot/${PREFIX}"
mkdir -p "$NANOPLOT_OUT"

NanoPlot \
  --fastq "$FASTQ" \
  --prefix "$PREFIX" \
  --outdir "$NANOPLOT_OUT" \
  --threads "$THREADS"

echo "[$(date +"%b %d %H:%M:%S")] NanoPlot done."

########################
# Flye assembly
########################
echo "[$(date +"%b %d %H:%M:%S")] Running Flye..."

FLYE_OUT="flye_assembly/${PREFIX}"
mkdir -p "$FLYE_OUT"

# flye \
#   --nano-hq "$FASTQ" \
#   --out-dir "$FLYE_OUT" \
#   --genome-size "$GENOME_SIZE" \
#   --threads "$THREADS"

ASSEMBLY="${FLYE_OUT}/${PREFIX}.assembly.fasta"
ln -s "$(pwd)/${FLYE_OUT}/assembly.fasta" "$ASSEMBLY"

echo "[$(date +"%b %d %H:%M:%S")] Flye done."
echo "Assembly: $ASSEMBLY"

########################
# QUAST
########################
echo "[$(date +"%b %d %H:%M:%S")] Running QUAST..."

QUAST_OUT="quast/${PREFIX}"
mkdir -p "$QUAST_OUT"

quast \
  "$ASSEMBLY" \
  -t "$THREADS" \
  -o "$QUAST_OUT"

echo "[$(date +"%b %d %H:%M:%S")] QUAST done."

########################
# BUSCO
########################
echo "[$(date +"%b %d %H:%M:%S")] Running BUSCO..."

BUSCO_OUT="busco/${PREFIX}"
mkdir -p busco

busco \
  -i "$ASSEMBLY" \
  -m genome \
  --lineage_dataset "$LINEAGE" \
  -c "$THREADS" \
  -o $BUSCO_OUT \
  -f

echo "[$(date +"%b %d %H:%M:%S")] BUSCO done."

########################
# Done
########################
echo "[$(date +"%b %d %H:%M:%S")] Pipeline complete."
