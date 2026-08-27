"""E3: accessory folded-SFS missingness sweep.

Recompute the accessory folded SFS at a series of max per-site missingness
cutoffs. If the accessory intermediate-frequency skew (which drives the positive
accessory Tajima's D) is a missingness artifact, it should relax toward the core
spectrum as the cutoff tightens; if it persists, some of it is real.

Reads the RAW merged VCF (all sites, pre-pruning/pre-MAF-filter) and subsets to
SV_* accessory contigs itself. Also emits the core (non-SV_*) spectrum once, as
the reference the accessory should relax toward. bcftools-only; parses the VCF
once (not once per cutoff).

Output columns: region, cutoff, bin_low, bin_high, count, prop, n_sites_kept.
(core has cutoff='none' -- reported unfiltered as the comparison baseline.)
"""

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

in_vcf = snakemake.input.vcf
cutoffs = list(snakemake.params.cutoffs)
edges = np.asarray(list(snakemake.params.bins), dtype=float)
out_path = Path(snakemake.output.sweep)
out_path.parent.mkdir(parents=True, exist_ok=True)


def af_missing():
    """One pass over the VCF: per-site (is_accessory, folded MAF, F_MISSING).
    +fill-tags recomputes AF/F_MISSING from genotypes so we don't rely on INFO."""
    cmd = (
        f"bcftools +fill-tags {in_vcf} -Ou -- -t AF,F_MISSING "
        f"| bcftools query -f '%CHROM\\t%INFO/AF\\t%INFO/F_MISSING\\n'"
    )
    proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    acc_maf, acc_miss, core_maf = [], [], []
    for line in proc.stdout.splitlines():
        chrom, af_s, miss_s = line.split("\t")
        try:
            af = float(af_s.split(",")[0])   # first ALT for multiallelics
            miss = float(miss_s)
        except ValueError:
            continue
        maf = min(af, 1.0 - af)
        if chrom.startswith("SV_"):
            acc_maf.append(maf)
            acc_miss.append(miss)
        else:
            core_maf.append(maf)
    return (np.array(acc_maf), np.array(acc_miss), np.array(core_maf))


def spectrum(maf):
    counts, _ = np.histogram(maf, bins=edges)
    total = counts.sum()
    return counts, total


acc_maf, acc_miss, core_maf = af_missing()

rows = []
# core baseline (unfiltered)
counts, total = spectrum(core_maf)
for i in range(len(counts)):
    rows.append({"region": "core", "cutoff": "none",
                 "bin_low": edges[i], "bin_high": edges[i + 1],
                 "count": int(counts[i]),
                 "prop": counts[i] / total if total else np.nan,
                 "n_sites_kept": len(core_maf)})
# accessory at each missingness cutoff
for cutoff in cutoffs:
    keep = acc_maf[acc_miss <= cutoff]
    counts, total = spectrum(keep)
    for i in range(len(counts)):
        rows.append({"region": "accessory", "cutoff": cutoff,
                     "bin_low": edges[i], "bin_high": edges[i + 1],
                     "count": int(counts[i]),
                     "prop": counts[i] / total if total else np.nan,
                     "n_sites_kept": len(keep)})

pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False, float_format="%.6f")
print(f"wrote {out_path}: core n={len(core_maf)}, accessory n={len(acc_maf)}")
