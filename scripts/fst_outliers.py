import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

selection_file = snakemake.input.selection
snp_coords_file = snakemake.input.snp_coords
out_outliers   = snakemake.output.outliers
fdr            = snakemake.params.fdr

scores = np.load(selection_file)
coords = pd.read_csv(snp_coords_file, sep="\t", names=["chrom", "pos", "id"])

pvals = stats.chi2.sf(scores, df=1)

n        = len(pvals)
ranks    = np.argsort(pvals)
sorted_p = pvals[ranks]
below    = sorted_p <= (np.arange(1, n + 1) / n) * fdr
cutoff   = sorted_p[below].max() if below.any() else 0.0

qvals_sorted = np.minimum(1.0, sorted_p * n / np.arange(1, n + 1))
qvals        = np.empty(n)
qvals[ranks] = qvals_sorted

coords["chi2"]    = scores
coords["pval"]    = pvals
coords["qval"]    = qvals
coords["outlier"] = pvals <= cutoff

Path(out_outliers).parent.mkdir(parents=True, exist_ok=True)
coords.sort_values(["chrom", "pos"]).to_csv(out_outliers, sep="\t", index=False)
