import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.outliers, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref
fdr      = snakemake.params.fdr

chroms        = sorted(df["chrom"].unique())
chrom_palette = sns.color_palette("tab20", len(chroms))
chrom_color   = {c: chrom_palette[i % len(chrom_palette)] for i, c in enumerate(chroms)}

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    sub = df[df["chrom"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub["pos"].max() - sub["pos"].min()) / 2
    offset += sub["pos"].max() - sub["pos"].min() + 1_000_000

df["x"] = df.apply(lambda r: offsets[r["chrom"]] + r["pos"], axis=1)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["chrom"] == chrom]
    ax.scatter(sub["x"], sub["chi2"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.6, linewidths=0)

outliers = df[df["outlier"]]
if len(outliers):
    ax.scatter(outliers["x"], outliers["chi2"],
               c="crimson", s=12, zorder=5, label="FDR < " + str(fdr))
    ax.axhline(outliers["chi2"].min(), color="crimson", lw=0.8, linestyle="--", alpha=0.7)
    ax.legend(frameon=False, fontsize=8)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Chi-squared score")
ax.set_title("Selection scan \u2014 " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
