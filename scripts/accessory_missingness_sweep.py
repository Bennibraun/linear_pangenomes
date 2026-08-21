"""E3: accessory folded-SFS missingness sweep.

Recompute the accessory site-frequency spectrum at a series of max per-site
missingness cutoffs. If the accessory SFS's intermediate-frequency skew (which
drives the positive accessory Tajima's D) is a missingness artifact, it should
relax toward the core spectrum as the cutoff tightens; if it persists at low
missingness, some of it is real shared variation.

Consumes the already-built accessory region VCF; writes a long-format TSV:
cutoff, bin_low, bin_high, count, prop, n_sites_kept.

bcftools-only (F_MISSING + folded MAF from AF), no plink -- so the augref
accessory contig count is a non-issue here.
"""

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

in_vcf = snakemake.input.vcf
cutoffs = list(snakemake.params.cutoffs)
bins = list(snakemake.params.bins)
out_path = Path(snakemake.output.sweep)
out_path.parent.mkdir(parents=True, exist_ok=True)

edges = np.asarray(bins, dtype=float)


def folded_maf_table(max_missing):
    """Return (maf_array) for sites passing F_MISSING <= max_missing.
    Uses bcftools to emit per-site AF and F_MISSING, computes folded MAF."""
    # +fill-tags recomputes AF/F_MISSING from genotypes so we don't depend on
    # them being present/current in the INFO field.
    cmd = (
        f"bcftools +fill-tags {in_vcf} -Ou -- -t AF,F_MISSING "
        f"| bcftools query -f '%INFO/AF\\t%INFO/F_MISSING\\n'"
    )
    proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    afs = []
    for line in proc.stdout.splitlines():
        af_s, miss_s = line.split("\t")
        # multiallelic AF can be comma-list; take the first ALT allele.
        af_s = af_s.split(",")[0]
        try:
            af = float(af_s)
            miss = float(miss_s)
        except ValueError:
            continue
        if miss <= max_missing:
            afs.append(min(af, 1.0 - af))  # fold
    return np.asarray(afs, dtype=float)


rows = []
for cutoff in cutoffs:
    maf = folded_maf_table(cutoff)
    n_kept = len(maf)
    counts, _ = np.histogram(maf, bins=edges)
    total = counts.sum()
    for i in range(len(counts)):
        rows.append({
            "cutoff": cutoff,
            "bin_low": edges[i],
            "bin_high": edges[i + 1],
            "count": int(counts[i]),
            "prop": float(counts[i] / total) if total else np.nan,
            "n_sites_kept": n_kept,
        })

pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False, float_format="%.6f")
print(f"wrote {out_path} for cutoffs {cutoffs}")
