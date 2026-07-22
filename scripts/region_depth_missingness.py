"""Per-sample depth, missingness, and mapping quality split core vs accessory
(SV_* contigs) on augref. This is the diagnostic for whether the core-vs-
accessory Tajima's D / FST gap reflects real shared variation or a calling
artifact (low coverage / poor mappability in accessory regions undercalls
rare variants, which mechanically inflates Tajima's D and can suppress FST).

Depth and MAPQ come from `samtools depth`/`samtools view` on the augref BAM,
split by contig using the .fai. Missingness comes from `bcftools query` GT
fields on the per-sample VCF, same split. One row per (sample, region).

Depth/MAPQ use a BED file (samtools -b/-L accepts one); missingness uses a
comma-joined -t/--targets argument (bcftools -R/-T require CHROM+POS columns
or a BED and reject a bare one-column contig list). Either way, contigs are
batched per region rather than spawning one subprocess per contig — augref
can carry hundreds to thousands of SV_* contigs.
"""
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

bam        = snakemake.input.bam
vcf        = snakemake.input.vcf
fai        = snakemake.input.fai
out_path   = Path(snakemake.output.tsv)
sample     = snakemake.params.sample

out_path.parent.mkdir(parents=True, exist_ok=True)

fai_rows = [line.split("\t") for line in Path(fai).read_text().splitlines() if line.strip()]
core_rows = [r for r in fai_rows if not r[0].startswith("SV_")]
acc_rows  = [r for r in fai_rows if r[0].startswith("SV_")]


def region_depth_mapq(fai_row_list):
    if not fai_row_list:
        return 0.0, 0.0
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
        for row in fai_row_list:
            f.write(f"{row[0]}\t0\t{row[1]}\n")
        bed_file = f.name
    try:
        total_sum, total_n = 0.0, 0
        depth_out = subprocess.run(
            ["samtools", "depth", "-a", "-b", bed_file, bam],
            capture_output=True, text=True, check=True,
        ).stdout
        for line in depth_out.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            total_sum += float(parts[2])
            total_n += 1

        mapq_sum, mapq_n = 0.0, 0
        mapq_out = subprocess.run(
            ["samtools", "view", "-F", "0x904", "-L", bed_file, bam],
            capture_output=True, text=True, check=True,
        ).stdout
        for line in mapq_out.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            mapq_sum += float(fields[4])
            mapq_n += 1
    finally:
        Path(bed_file).unlink(missing_ok=True)

    mean_depth = total_sum / total_n if total_n else 0.0
    mean_mapq  = mapq_sum / mapq_n if mapq_n else 0.0
    return mean_depth, mean_mapq


def region_missingness(fai_row_list):
    if not fai_row_list:
        return 0.0, 0
    # -t/--targets with a comma-separated region string (not -R/-T with a
    # file): both -R and -T require CHROM+POS columns or a BED and reject a
    # bare one-column contig list ("Could not parse ... using the columns
    # 1,2[,-1]"). -t/-r are documented to accept plain region names directly
    # in the argument. Passed as a single argv element (not through a shell),
    # so OS ARG_MAX (~2MB on Linux) is the only limit — thousands of contig
    # names fit comfortably under that.
    keep_contigs = ",".join(row[0] for row in fai_row_list)
    out = subprocess.run(
        ["bcftools", "query", "-t", keep_contigs, "-f", "[%GT]\n", vcf],
        capture_output=True, text=True, check=True,
    ).stdout
    n_sites = 0
    n_missing = 0
    for line in out.splitlines():
        gt = line.strip()
        if not gt:
            continue
        n_sites += 1
        if "." in gt:
            n_missing += 1
    missing_rate = n_missing / n_sites if n_sites else 0.0
    return missing_rate, n_sites


rows = []
for region, row_list in (("core", core_rows), ("accessory", acc_rows)):
    mean_depth, mean_mapq = region_depth_mapq(row_list)
    missing_rate, n_sites = region_missingness(row_list)
    rows.append({
        "sample": sample,
        "region": region,
        "n_contigs": len(row_list),
        "mean_depth": mean_depth,
        "mean_mapq": mean_mapq,
        "n_sites": n_sites,
        "missing_rate": missing_rate,
    })

pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
