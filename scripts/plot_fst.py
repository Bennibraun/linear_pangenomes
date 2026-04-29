import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref
pop1     = snakemake.params.pop1
pop2     = snakemake.params.pop2

# vcftools writes "-nan" (and sometimes "nan") for undefined FST values; these
# are not in pandas' default NA list, so declare them explicitly. Without this,
# the FST column becomes object-dtype and downstream comparisons/plotting
# either error or crawl.
df = pd.read_csv(
    snakemake.input.fst,
    sep="\t",
    na_values=["-nan", "nan", "-NaN", "NaN", "-inf", "inf"],
)

is_windowed = "BIN_START" in df.columns
pos_col     = "BIN_START" if is_windowed else "POS"
fst_col     = "WEIGHTED_FST" if is_windowed else "WEIR_AND_COCKERHAM_FST"

# Force numeric (in case any unrecognized sentinel slipped through) and drop
# rows where FST is NaN — these are uninformative sites with no genotypes
# called in one or both populations.
df[fst_col] = pd.to_numeric(df[fst_col], errors="coerce")
df = df.dropna(subset=[fst_col, pos_col])

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {c: palette[i % len(palette)] for i, c in enumerate(chroms)}

# Compute per-chrom offsets via groupby (vectorized, single pass).
chrom_min = df.groupby("CHROM")[pos_col].min()
chrom_max = df.groupby("CHROM")[pos_col].max()
chrom_span = (chrom_max - chrom_min).reindex(chroms)

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    span = chrom_span[chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + span / 2
    offset += span + 1_000_000

# Vectorized x-coordinate computation — replaces a row-wise apply that was
# the dominant cost on genome-wide windowed FST tables.
df["x"] = df["CHROM"].map(offsets).astype(float) + df[pos_col].astype(float)
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
