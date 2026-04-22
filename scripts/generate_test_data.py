#!/usr/bin/env python3
"""
Generate a tiny dummy test dataset for the linear-pangenomes pipeline.

Creates:
  test_data/
    references/
      conspec.fasta + .fai   (5 chromosomes, varying lengths)
      hetspec.fasta + .fai
    reads/
      <sample>_R1/R2.fastq.gz  (15 short-read samples, 5 per group × 3 groups)
      <sample>.fastq.gz         (6 long-read samples, 2 per group × 3 groups)
    reads_manifest.tsv
    config_test.yaml

Usage:
    python scripts/generate_test_data.py [--outdir test_data]
"""

import argparse
import gzip
import random
import subprocess
import sys
from pathlib import Path

# Reference: 5 chromosomes with realistic relative size differences
CHROMS = [
    ("chr1", 80_000),
    ("chr2", 60_000),
    ("chr3", 45_000),
    ("chr4", 30_000),
    ("chr5", 20_000),
]

READ_LEN    = 150    # short-read length (bp)
LONG_RL     = 5_000  # long-read length (bp)
N_SHORT     = 400    # read pairs per short sample
N_LONG      = 80     # reads per long sample
RANDOM_SEED = 42

# 3 groups × 5 short samples + 3 groups × 2 long samples
GROUPS = ["popA", "popB", "popC"]
N_SHORT_PER_GROUP = 5
N_LONG_PER_GROUP  = 2

# Per-group SNP divergence relative to conspec (popA is the reference-like group)
GROUP_DIVERGENCE = {"popA": 0.002, "popB": 0.008, "popC": 0.015}
# Within-sample additional noise on top of group background
SAMPLE_NOISE = 0.001

random.seed(RANDOM_SEED)

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def rand_seq(n: int) -> str:
    return "".join(random.choices("ACGT", k=n))


def mutate(seq: str, rate: float) -> str:
    bases = list(seq)
    for i in range(len(bases)):
        if random.random() < rate:
            bases[i] = random.choice("ACGT")
    return "".join(bases)


def write_fasta(path: Path, chroms: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in chroms:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def write_fai(fasta: Path, chroms: list[tuple[str, str]]) -> None:
    fai = Path(str(fasta) + ".fai")
    with open(fai, "w") as fh:
        offset = 0
        line_bases = 60
        line_bytes  = 61  # 60 bases + newline
        for name, seq in chroms:
            header_len = len(name) + 2  # ">" + name + "\n"
            seq_len    = len(seq)
            n_full     = seq_len // line_bases
            remainder  = seq_len %  line_bases
            fh.write(f"{name}\t{seq_len}\t{offset + header_len}\t{line_bases}\t{line_bytes}\n")
            offset += header_len + n_full * line_bytes + (remainder + 1 if remainder else 0)


def write_fastq_gz(path: Path, reads: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for i, (seq, qual) in enumerate(reads):
            fh.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def fake_qual(n: int, min_q: int = 30, max_q: int = 40) -> str:
    return "".join(chr(random.randint(min_q, max_q) + 33) for _ in range(n))


def sample_short_reads(
    chrom_seqs: list[str], n: int, rl: int
) -> tuple[list, list]:
    """Draw n read pairs uniformly across all chromosomes."""
    total_len = sum(len(s) for s in chrom_seqs)
    weights   = [len(s) / total_len for s in chrom_seqs]
    r1_reads, r2_reads = [], []
    for _ in range(n):
        chrom_seq = random.choices(chrom_seqs, weights=weights, k=1)[0]
        if len(chrom_seq) <= rl:
            continue
        start = random.randint(0, len(chrom_seq) - rl - 1)
        fwd = mutate(chrom_seq[start:start+rl], SAMPLE_NOISE)
        rev = mutate(revcomp(chrom_seq[start:start+rl]), SAMPLE_NOISE)
        r1_reads.append((fwd, fake_qual(rl)))
        r2_reads.append((rev, fake_qual(rl)))
    return r1_reads, r2_reads


def sample_long_reads(chrom_seqs: list[str], n: int, rl: int) -> list:
    total_len = sum(len(s) for s in chrom_seqs)
    weights   = [len(s) / total_len for s in chrom_seqs]
    reads = []
    for _ in range(n):
        chrom_seq = random.choices(chrom_seqs, weights=weights, k=1)[0]
        if len(chrom_seq) <= rl:
            rl = len(chrom_seq) - 1
        start = random.randint(0, len(chrom_seq) - rl - 1)
        seq = mutate(chrom_seq[start:start+rl], rate=0.02)
        reads.append((seq, fake_qual(rl, min_q=5, max_q=20)))
    return reads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="test_data")
    args = parser.parse_args()

    outdir    = Path(args.outdir)
    ref_dir   = outdir / "references"
    reads_dir = outdir / "reads"

    print(f"Generating test dataset in: {outdir.resolve()}")

    # ------------------------------------------------------------------
    # 1. Reference genomes (5 chromosomes each)
    # ------------------------------------------------------------------
    conspec_chroms = [(name, rand_seq(length)) for name, length in CHROMS]
    hetspec_chroms = [(name, mutate(seq, rate=0.03)) for name, seq in conspec_chroms]

    conspec_fa = ref_dir / "conspec.fasta"
    hetspec_fa = ref_dir / "hetspec.fasta"

    write_fasta(conspec_fa, conspec_chroms)
    write_fai(conspec_fa, conspec_chroms)

    write_fasta(hetspec_fa, hetspec_chroms)
    write_fai(hetspec_fa, hetspec_chroms)

    print(f"  References written ({len(CHROMS)} chromosomes each).")

    # Build per-group background sequences (diverged from conspec)
    group_chrom_seqs = {
        grp: [mutate(seq, GROUP_DIVERGENCE[grp]) for _, seq in conspec_chroms]
        for grp in GROUPS
    }

    # ------------------------------------------------------------------
    # 2. Short-read samples — 5 per group × 3 groups = 15 total
    # ------------------------------------------------------------------
    manifest_rows = []

    for grp in GROUPS:
        bg_seqs = group_chrom_seqs[grp]
        for i in range(1, N_SHORT_PER_GROUP + 1):
            sample = f"{grp}_S{i:02d}"
            r1_reads, r2_reads = sample_short_reads(bg_seqs, N_SHORT, READ_LEN)
            r1_path = reads_dir / f"{sample}_R1.fastq.gz"
            r2_path = reads_dir / f"{sample}_R2.fastq.gz"
            write_fastq_gz(r1_path, r1_reads)
            write_fastq_gz(r2_path, r2_reads)
            manifest_rows.append(
                f"{sample}\tshort\tILLUMINA\t{grp}"
                f"\t{r1_path.resolve()}\t{r2_path.resolve()}"
            )
        print(f"  Short reads: {grp} ({N_SHORT_PER_GROUP} samples)")

    # ------------------------------------------------------------------
    # 3. Long-read samples — 2 per group × 3 groups = 6 total
    #    Alternate ONT / PACBIO within each group
    # ------------------------------------------------------------------
    platforms = ["ONT", "PACBIO"]
    for grp in GROUPS:
        bg_seqs = group_chrom_seqs[grp]
        for i in range(1, N_LONG_PER_GROUP + 1):
            sample   = f"{grp}_L{i:02d}"
            platform = platforms[(i - 1) % len(platforms)]
            reads    = sample_long_reads(bg_seqs, N_LONG, LONG_RL)
            fq_path  = reads_dir / f"{sample}.fastq.gz"
            write_fastq_gz(fq_path, reads)
            manifest_rows.append(
                f"{sample}\tlong\t{platform}\t{grp}"
                f"\t{fq_path.resolve()}\t"
            )
        print(f"  Long reads:  {grp} ({N_LONG_PER_GROUP} samples)")

    # ------------------------------------------------------------------
    # 4. Reads manifest
    # ------------------------------------------------------------------
    manifest_path = outdir / "reads_manifest.tsv"
    header = "sample_id\tseq_type\tplatform\tgrouping\tfastq_r1\tfastq_r2"
    manifest_path.write_text(header + "\n" + "\n".join(manifest_rows) + "\n")
    print(f"  Manifest written: {manifest_path}")

    # ------------------------------------------------------------------
    # 5. Test config — all pairs among 3 groups
    # ------------------------------------------------------------------
    pop_pairs_yaml = "\n".join(
        f'  - ["{a}", "{b}"]'
        for i, a in enumerate(GROUPS)
        for b in GROUPS[i+1:]
    )

    config_text = f"""\
# Test configuration — generated by generate_test_data.py
# 5 chromosomes, 15 short-read samples (5×3 groups), 6 long-read samples (2×3 groups)

inputs:
  reads_manifest: "{manifest_path.resolve()}"

references:
  augref:
    fasta: "results_test/sv_calls/augref/augmented_reference.fasta"
  conspec:
    fasta: "{conspec_fa.resolve()}"
  hetspec:
    fasta: "{hetspec_fa.resolve()}"

population_pairs:
{pop_pairs_yaml}

assembly:
  outdir: "results_test/assemblies"
  lineage: "hymenoptera_odb10"
  threads: 2

sv_calling:
  outdir: "results_test/sv_calls"
  assembly_dir: "results_test/assemblies/hifiasm"
  threads: 2
  min_sv_size: 50
  min_read_support: 1
  breakpoint_slop: 500
  jasmine_slop: 300
  flank: 100

cactus:
  image: "docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.9"
  outname: "cactus_graph"
  max_cores: 2
  ref_contigs: ""
  extra_args: ""

align_wgs:
  outdir: "results_test/wgs_alignments"
  threads: 2

align_metrics:
  outdir: "results_test/align_metrics"
  include_gam: true
  min_mapq: 0
  min_baseq: 0
  threads: 2

variant_calling:
  outdir: "results_test/variants"
  bcftools:
    threads: 2
  vg:
    threads: 2

fst_afs:
  outdir: "results_test/fst_afs"
  window_size: 0
  window_step: 0
  afs_bins: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

pi:
  outdir: "results_test/pi"
  window_size: 5000
  window_step: 2500

allelic_balance:
  outdir: "results_test/allelic_balance"
  bins: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

pcangsd:
  outdir: "results_test/pcangsd"
  ld_window: 10
  ld_step: 5
  ld_r2: 0.3
  maf: 0.05
  n_pcs: 3

plot:
  outdir: "results_test/plots"

roh:
  outdir: "results_test/roh"
  genome_lengths: {{}}
  min_roh_length: 10000
  bcftools_args: ""
"""

    config_path = outdir / "config_test.yaml"
    config_path.write_text(config_text)
    print(f"  Config written: {config_path}")

    # ------------------------------------------------------------------
    # 6. Summary
    # ------------------------------------------------------------------
    n_short = len(GROUPS) * N_SHORT_PER_GROUP
    n_long  = len(GROUPS) * N_LONG_PER_GROUP
    chrom_summary = ", ".join(f"{n} ({l:,} bp)" for n, l in CHROMS)
    print()
    print("Done. To run the pipeline with this test dataset:")
    print(f"  snakemake --configfile {config_path.resolve()} --cores 4 --use-conda")
    print()
    print(f"Groups:        {', '.join(GROUPS)}")
    print(f"Short samples: {n_short} ({N_SHORT_PER_GROUP} per group)")
    print(f"Long  samples: {n_long} ({N_LONG_PER_GROUP} per group, alternating ONT/PACBIO)")
    print(f"Chromosomes:   {chrom_summary}")


if __name__ == "__main__":
    main()
