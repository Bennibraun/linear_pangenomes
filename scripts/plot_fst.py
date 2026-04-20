import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.fst, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref
pop1     = snakemake.params.pop1
pop2     = snakemake.params.pop2

is_windowed = "BIN_START" in df.columns
pos_col     = "BIN_START" if is_windowed else "POS"
fst_col     = "WEIGHTED_FST" if is_windowed else "WEIR_AND_COCKERHAM_FST"

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {c: palette[i % len(palette)] for i, c in enumerate(chroms)}

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub[pos_col].max() - sub[pos_col].min()) / 2
    offset += sub[pos_col].max() - sub[pos_col].min() + 1_000_000

df["x"] = df.apply(lambda r: offsets[r["CHROM"]] + r[pos_col], axis=1)
df_plot = df[df[fst_col] >= 0]

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df_plot[df_plot["CHROM"] == chrom]
    ax.scatter(sub["x"], sub[fst_col],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("FST")
ax.set_ylim(bottom=0)
ax.set_title("FST \u2014 " + ref_name + ": " + pop1 + " vs " + pop2, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
