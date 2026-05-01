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

# Drop augref SV alt contigs and tiny scaffolds (see plot_pi.py / plot_fst.py).
df = df[~df["chrom"].astype(str).str.startswith("SV_")]
chrom_span_all = df.groupby("chrom")["pos"].agg(lambda s: s.max() - s.min())
keep_chroms    = chrom_span_all[chrom_span_all >= 1_000_000].index.tolist()
df             = df[df["chrom"].isin(keep_chroms)]

if df.empty:
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.text(0.5, 0.5, f"No qualifying chromosomes for selection scan — {ref_name}",
            ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)

chroms      = sorted(df["chrom"].unique())
palette_arr = np.array(sns.color_palette("tab20", len(chroms)))

chrom_grp  = df.groupby("chrom")["pos"]
chrom_min  = chrom_grp.min()
chrom_max  = chrom_grp.max()
chrom_span = (chrom_max - chrom_min).reindex(chroms)

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    span = chrom_span[chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + span / 2
    offset += span + 1_000_000

df["x"] = df["chrom"].map(offsets).astype(float) + df["pos"].astype(float)

# Thin to plotting resolution keeping max chi2 per bin so outlier peaks survive.
x_range  = df["x"].max() - df["x"].min()
bin_size = max(1.0, x_range / 4800.0)
df["_bin"] = (df["x"] / bin_size).astype(np.int64)
df = df.loc[df.groupby(["chrom", "_bin"])["chi2"].idxmax()].drop(columns="_bin")

chrom_to_idx = {c: i for i, c in enumerate(chroms)}
colors = palette_arr[df["chrom"].map(chrom_to_idx).values % len(palette_arr)]

fig, ax = plt.subplots(figsize=(16, 4))
ax.scatter(df["x"].values, df["chi2"].values,
           c=colors, s=4, alpha=0.6, linewidths=0, rasterized=True)

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
ax.set_title("Selection scan — " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
