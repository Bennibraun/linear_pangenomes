import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

selection_file  = snakemake.input.selection
snp_coords_file = snakemake.input.snp_coords
out_outliers    = snakemake.output.outliers
fdr             = snakemake.params.fdr

# pcangsd --selection writes a 2D array of shape (n_sites, n_pcs) as plain
# text (since pcangsd v1.x; earlier versions wrote .npy). Each column is an
# independent chi-squared statistic with df=1.
scores = np.loadtxt(selection_file, ndmin=2)
coords = pd.read_csv(snp_coords_file, sep="\t", names=["chrom", "pos", "id"])

if scores.ndim == 1:
    scores = scores[:, None]

if scores.shape[0] != len(coords):
    raise ValueError(
        f"Selection scores ({scores.shape}) and snp_coords ({len(coords)}) row counts disagree"
    )

# Per-PC p-values, then aggregate. We take the minimum p across PCs for the
# "any-PC" outlier flag and Bonferroni-correct across PCs to keep the
# combined null calibrated. Per-PC chi2/p columns are also retained so users
# can inspect which axis drives a hit.
n_sites, n_pcs = scores.shape
pvals_per_pc   = stats.chi2.sf(scores, df=1)
chi2_max       = scores.max(axis=1)
pc_argmax      = scores.argmax(axis=1) + 1  # 1-indexed for readability
pmin           = pvals_per_pc.min(axis=1)
pcomb          = np.minimum(1.0, pmin * n_pcs)  # Bonferroni across PCs

# BH FDR on the across-PC combined p-values.
ranks    = np.argsort(pcomb)
sorted_p = pcomb[ranks]
n        = len(pcomb)
below    = sorted_p <= (np.arange(1, n + 1) / n) * fdr
cutoff   = sorted_p[below].max() if below.any() else 0.0

qvals_sorted = np.minimum.accumulate((sorted_p * n / np.arange(1, n + 1))[::-1])[::-1]
qvals_sorted = np.minimum(1.0, qvals_sorted)
qvals        = np.empty(n)
qvals[ranks] = qvals_sorted

coords["chi2"]      = chi2_max
coords["best_pc"]   = pc_argmax
coords["pval"]      = pcomb
coords["qval"]      = qvals
coords["outlier"]   = pcomb <= cutoff

for j in range(n_pcs):
    coords[f"chi2_pc{j+1}"] = scores[:, j]
    coords[f"pval_pc{j+1}"] = pvals_per_pc[:, j]

Path(out_outliers).parent.mkdir(parents=True, exist_ok=True)
coords.sort_values(["chrom", "pos"]).to_csv(out_outliers, sep="\t", index=False)
