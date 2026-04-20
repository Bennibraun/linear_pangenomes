import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.pi, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {c: palette[i % len(palette)] for i, c in enumerate(chroms)}

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub["BIN_START"].max() - sub["BIN_START"].min()) / 2
    offset += sub["BIN_START"].max() - sub["BIN_START"].min() + 1_000_000

df["x"] = df.apply(lambda r: offsets[r["CHROM"]] + r["BIN_START"], axis=1)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    ax.scatter(sub["x"], sub["PI"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("\u03c0")
ax.set_title("Nucleotide diversity \u2014 " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
