# Linear Pangenomes

A Snakemake pipeline for building pangenome graphs from long-read assemblies and analyzing population genetics metrics across multiple references.

## Overview

This pipeline produces a pangenome graph from multiple de novo ONT assemblies, aligns short-read sequences to multiple reference genomes (linear and graph-based), and computes population genetics statistics including FST, nucleotide diversity, allelic balance, and runs of homozygosity.

## Inputs

- **Long reads (ONT)**: FASTQ files for de novo assembly
- **Short reads (Illumina)**: Paired-end FASTQ files for alignment
- **Reference genomes**: FASTA files for augmented reference, conspecific, and heterospecific comparisons
- **Cactus seqfile**: Text file defining samples for pangenome graph construction
- **Population definitions**: Sample lists for FST/AFS population comparisons

## Outputs

### Assembly Stage
- `results/assemblies/flye/{sample}/assembly.fasta` - de novo assemblies
- `results/assemblies/nanostat/` - read statistics
- `results/assemblies/quast/` - assembly quality metrics
- `results/assemblies/busco/` - completeness assessment

### SV Calling
- `results/sv_calls/pan_sample_catalog/` - merged SV catalogs (SURVIVOR, Jasmine)
- `results/sv_calls/augref/augmented_reference.fasta` - reference with SV insertions

### Pangenome Graph
- `results/cactus/{name}.gbz` - graph bundle (giraffe index)
- `results/cactus/{name}.gfa.gz` - graph in GFA format
- `results/cactus/{name}.vcf.gz` - structural variants in VCF

### Alignments
- `results/wgs_alignments/{sample}/{sample}.{ref}.bam` - BWA alignments (augref, conspec, hetspec)
- `results/wgs_alignments/{sample}/{sample}.cactus.gam` - VG Giraffe alignments to graph

### Variants
- `results/variants/gatk/{ref}/{sample}.vcf.gz` - SNVs per reference
- `results/variants/gatk/{ref}/merged.vcf.gz` - merged VCFs per reference
- `results/variants/vg/{sample}.vcf.gz` - graph-native variants

### Population Genetics
- `results/align_metrics/alignment_metrics.tsv` - mapping statistics across all samples/references
- `results/fst_afs/{ref}/*.weir.fst` - FST between population pairs
- `results/fst_afs/afs/{ref}.afs.tsv` - allele frequency spectrum
- `results/pi/{ref}.windowed.pi` - nucleotide diversity (windowed)
- `results/allelic_balance/{ref}/{sample}.allelic_balance.tsv` - ref/alt ratios at heterozygous sites
- `results/roh/f_roh_summary.tsv` - inbreeding coefficient per sample

## Dependencies

All tools are installed via conda environments defined in `envs/`. Main dependencies:

- **Assembly**: Flye, QUAST, BUSCO, NanoStat, NanoPlot
- **SV calling**: minimap2, Sniffles2, cuteSV, SURVIVOR, Jasmine, cd-hit
- **Graph construction**: Cactus (via Singularity container)
- **Alignment**: BWA, VG (Giraffe)
- **Variant calling**: GATK4, VG
- **Population genetics**: vcftools, bcftools, pysam

## Usage

### Setup

1. Clone and configure:
```bash
git clone <repo>
cd linear_pangenomes
```

2. Edit `config/config.yaml`:
   - Set `assembly.samples` and `assembly.reads_dir`
   - Set `sv_calling.reference` and `align_wgs` reference paths
   - Set `cactus.seqfile` and `cactus.jobstore` paths
   - Define populations in `fst_afs.populations`

### Run

```bash
# Dry run to check graph
snakemake -n

# Run entire pipeline
snakemake -j 8 --use-conda

# Run specific stage
snakemake results/assemblies/flye/ -j 8 --use-conda
snakemake results/sv_calls/ -j 8 --use-conda
snakemake results/variants/ -j 8 --use-conda
```

### Configuration Notes

- **Genome lengths** for ROH calculation are auto-generated from reference `.fai` files
- **Reference indices** (BWA, samtools) are automatically created
- **Default settings** are tuned for honeybee-sized genomes (~225 Mb)
- **FLANK parameter** (default 200 bp) controls flanking sequence around SVs in augmented reference

## Pipeline Stages

1. **Assembly + QC**: Flye assembly with NanoStat, NanoPlot, QUAST, BUSCO
2. **SV Calling**: Multi-caller (Sniffles2, cuteSV) merged with SURVIVOR/Jasmine
3. **Graph Construction**: Cactus pangenome with augmented reference
4. **Alignment**: 3 linear references + graph-based alignment
5. **Metrics**: Alignment statistics, variant calling
6. **Population Genetics**: FST, π, allelic balance, F_ROH

## Citation

Built with Snakemake and Cactus. See individual tool papers for details.
