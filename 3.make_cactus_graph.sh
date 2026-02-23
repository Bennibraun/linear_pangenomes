#!/bin/bash
#SBATCH --job-name=cactus_apis
#SBATCH --output=cactus_apis.%j.out
#SBATCH --error=cactus_apis.%j.err
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --partition=highmem
#SBATCH -N 1

# Example script for running cactus - currently uses the A. mellifera data stored on the Chuong lab share

module load singularity

cd /Users/bebr1814/projects/chuong_rotation/data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/cactus

# singularity pull docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.9

# https://github.com/ComparativeGenomicsToolkit/cactus/blob/v2.9.9/doc/pangenome.md
# cactus-pangenome <jobStorePath> <seqFile> --outDir <output directory> --outName <output file prefix> --reference <reference sample name>
# consider specifying the chromosomes with --refContigs
singularity exec --bind /Users/bebr1814/scratch/chuong_data:/chuong_data /Users/bebr1814/software/cactus/cactus_v2.9.9.sif cactus-pangenome /chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/cactus/work /chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/cactus/seqfile.txt --reference "amelhav31" --collapse --outDir /chuong_data/20250211_1135_P2S-00613-A_PBA08559_66fdcccd/uncorrected_fastq/cactus/output --outName "cactus_test_1" --giraffe --vcf --gfa --gbz --maxCores 8


