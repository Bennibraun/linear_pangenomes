"""Pairwise mean-FST heatmap across all population pairs for one reference.

Each per-pair vcftools FST file is reduced to a single genome-wide mean FST
(core genome, SV_* excluded), then arranged into a symmetric population ×
population matrix. This is the population-structure summary: one glance shows
which populations are most differentiated, replacing N separate per-pair
manhattan plots. A tidy pairwise table is written alongside.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

fst_files = list(snakemake.input.fsts)
pairs     = [tuple(p) for p in snakemake.params.pairs]
out_pdf   = snakemake.output.pdf
out_png   = snakemake.output.png
out_table = snakemake.output.table
ref_name  = snakemake.params.ref

NA = ["-nan", "nan", "-NaN", "NaN", "-inf", "inf"]


def mean_fst(path):
    df = pd.read_csv(path, sep="\t", na_values=NA)
    fst_col = "WEIGHTED_FST" if "WEIGHTED_FST" in df.columns \
        else "WEIR_AND_COCKERHAM_FST"
    if fst_col not in df.columns:
        return np.nan
    df[fst_col] = pd.to_numeric(df[fst_col], errors="coerce")
    df = df[~df["CHROM"].astype(str).str.startswith("SV_")]
    vals = df[fst_col].dropna()
    # FST is bounded at 0; vcftools can emit small negatives from estimator
    # noise. Clip to 0 for a cleaner mean/heatmap.
    vals = vals.clip(lower=0)
    return float(vals.mean()) if len(vals) else np.nan


rows = []
for (p1, p2), path in zip(pairs, fst_files):
    rows.append({"pop1": p1, "pop2": p2, "mean_fst": mean_fst(path)})
tidy = pd.DataFrame(rows)
tidy.to_csv(out_table, sep="\t", index=False)

pops = sorted(set([p for pair in pairs for p in pair]))
mat = pd.DataFrame(np.nan, index=pops, columns=pops, dtype=float)
for _, r in tidy.iterrows():
    mat.loc[r["pop1"], r["pop2"]] = r["mean_fst"]
    mat.loc[r["pop2"], r["pop1"]] = r["mean_fst"]
# Set the diagonal via pandas — mat.values can be a read-only view, so
# np.fill_diagonal on it raises "underlying array is read-only".
for p in pops:
    mat.loc[p, p] = 0.0

fig, ax = plt.subplots(figsize=(1.1 * len(pops) + 3, 1.1 * len(pops) + 2))
sns.heatmap(mat, annot=True, fmt=".3f", cmap="viridis", square=True,
            linewidths=0.5, linecolor="white",
            cbar_kws={"label": "mean FST"}, ax=ax,
            annot_kws={"fontsize": 7})
ax.set_title(f"Pairwise mean FST — {ref_name}", fontweight="bold")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
