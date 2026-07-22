"""Per-sample depth, missingness, and mapping quality split core vs accessory
(SV_* contigs) on augref. This is the diagnostic for whether the core-vs-
accessory Tajima's D / FST gap reflects real shared variation or a calling
artifact (low coverage / poor mappability in accessory regions undercalls
rare variants, which mechanically inflates Tajima's D and can suppress FST).

Depth and MAPQ come from `samtools depth`/`samtools view` on the augref BAM,
split by contig using the .fai. Missingness comes from `bcftools query` GT
fields on the per-sample VCF, same split. One row per (sample, region).
"""
import subprocess
from pathlib import Path

import pandas as pd

bam        = snakemake.input.bam
vcf        = snakemake.input.vcf
fai        = snakemake.input.fai
out_path   = Path(snakemake.output.tsv)
sample     = snakemake.params.sample

out_path.parent.mkdir(parents=True, exist_ok=True)

contigs = [line.split("\t")[0] for line in Path(fai).read_text().splitlines() if line.strip()]
core_contigs = [c for c in contigs if not c.startswith("SV_")]
acc_contigs  = [c for c in contigs if c.startswith("SV_")]


def region_depth_mapq(contig_list):
    if not contig_list:
        return 0.0, 0.0
    total_sum, total_n = 0.0, 0
    mapq_sum, mapq_n = 0.0, 0
    # Batch contigs through one samtools view/depth call each to avoid
    # spawning one subprocess per contig on refs with many SV_* entries.
    for c in contig_list:
        depth_out = subprocess.run(
            ["samtools", "depth", "-a", "-r", c, bam],
            capture_output=True, text=True, check=True,
        ).stdout
        for line in depth_out.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            total_sum += float(parts[2])
            total_n += 1

        mapq_out = subprocess.run(
            ["samtools", "view", "-F", "0x904", bam, c],
            capture_output=True, text=True, check=True,
        ).stdout
        for line in mapq_out.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            mapq_sum += float(fields[4])
            mapq_n += 1

    mean_depth = total_sum / total_n if total_n else 0.0
    mean_mapq  = mapq_sum / mapq_n if mapq_n else 0.0
    return mean_depth, mean_mapq


def region_missingness(contig_list):
    if not contig_list:
        return 0.0, 0
    regions_arg = ",".join(contig_list)
    out = subprocess.run(
        ["bcftools", "query", "-r", regions_arg, "-f", "[%GT]\n", vcf],
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
for region, clist in (("core", core_contigs), ("accessory", acc_contigs)):
    mean_depth, mean_mapq = region_depth_mapq(clist)
    missing_rate, n_sites = region_missingness(clist)
    rows.append({
        "sample": sample,
        "region": region,
        "n_contigs": len(clist),
        "mean_depth": mean_depth,
        "mean_mapq": mean_mapq,
        "n_sites": n_sites,
        "missing_rate": missing_rate,
    })

pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
