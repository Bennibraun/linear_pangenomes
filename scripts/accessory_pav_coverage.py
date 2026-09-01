"""Depth-aware presence/absence (PAV) of each accessory contig, per sample.

The old PAV called a contig "present" if a sample had >=1 non-missing genotype
on it. That is depth-confounded: accessory depth is ~2x higher in cave than
surface samples (cave-flavoured content vs a surface reference), so surface
samples can carry a sequence yet fail to genotype it and be scored absent -- a
coverage artifact, not real PAV (this is exactly how the OCA2 false positive
arose).

Here presence is a NORMALISED coverage call:

    norm_cov(contig, sample) = mean_depth(contig) / genome_baseline(sample)

where genome_baseline is the sample's own mean depth on the core chromosomes.
Dividing by each sample's own baseline cancels the per-sample depth differences
(including the cave/surface gap), so the threshold compares like with like. A
contig is present when a fraction >= `min_breadth` of its length is covered AND
its normalised coverage >= `min_norm_cov`.

Output: one row per (sample, contig) with raw and normalised coverage, breadth,
and a boolean `present`, plus a wide presence matrix (contigs x samples) for
downstream PAV PCA / FST.
"""
import subprocess
from pathlib import Path

import pandas as pd

bam          = snakemake.input.bam
fai          = snakemake.input.fai
out_long     = Path(snakemake.output.long_tsv)
out_matrix   = Path(snakemake.output.matrix_tsv)
sample       = snakemake.params.sample
min_breadth  = float(snakemake.params.min_breadth)
min_norm_cov = float(snakemake.params.min_norm_cov)

out_long.parent.mkdir(parents=True, exist_ok=True)

fai_rows = [ln.split("\t") for ln in Path(fai).read_text().splitlines() if ln.strip()]
acc = {r[0]: int(r[1]) for r in fai_rows if r[0].startswith(("SV_", "UNMAP_"))}
core_names = [r[0] for r in fai_rows if not r[0].startswith(("SV_", "UNMAP_"))]


def write_bed(names_lengths):
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as f:
        for name, length in names_lengths:
            f.write(f"{name}\t0\t{length}\n")
        return f.name


def per_contig_depth(bed_file):
    """samtools depth -a over a BED -> {contig: (sum_depth, covered_bases)}."""
    agg = {}
    with subprocess.Popen(
        ["samtools", "depth", "-a", "-b", bed_file, bam],
        stdout=subprocess.PIPE, text=True,
    ) as proc:
        for line in proc.stdout:
            c, _pos, d = line.rstrip("\n").split("\t")
            d = int(d)
            s, cov = agg.get(c, (0, 0))
            agg[c] = (s + d, cov + (1 if d > 0 else 0))
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)
    return agg


# Sample baseline: mean depth over core chromosomes.
fai_len = {r[0]: int(r[1]) for r in fai_rows}
core_bed = write_bed([(n, fai_len[n]) for n in core_names])
try:
    core_depth = per_contig_depth(core_bed)
finally:
    Path(core_bed).unlink(missing_ok=True)
core_sum = sum(s for s, _ in core_depth.values())
core_bp = sum(fai_len[n] for n in core_names)
baseline = (core_sum / core_bp) if core_bp else 0.0

# Accessory per-contig coverage.
acc_bed = write_bed(list(acc.items()))
try:
    acc_depth = per_contig_depth(acc_bed)
finally:
    Path(acc_bed).unlink(missing_ok=True)

rows = []
for contig, length in acc.items():
    sum_d, covered = acc_depth.get(contig, (0, 0))
    mean_d = sum_d / length if length else 0.0
    breadth = covered / length if length else 0.0
    norm_cov = (mean_d / baseline) if baseline else 0.0
    present = (breadth >= min_breadth) and (norm_cov >= min_norm_cov)
    rows.append({
        "sample": sample,
        "contig": contig,
        "length": length,
        "mean_depth": round(mean_d, 4),
        "baseline_depth": round(baseline, 4),
        "norm_cov": round(norm_cov, 4),
        "breadth": round(breadth, 4),
        "present": int(present),
    })

df = pd.DataFrame(rows)
df.to_csv(out_long, sep="\t", index=False)
# Wide matrix: contig x {this sample} presence (merged across samples later).
df[["contig", "present"]].rename(columns={"present": sample}).to_csv(
    out_matrix, sep="\t", index=False
)
print(f"{sample}: {int(df['present'].sum())}/{len(df)} accessory contigs present "
      f"(baseline {baseline:.2f}x)")
