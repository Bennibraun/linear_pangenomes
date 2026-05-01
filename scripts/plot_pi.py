import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.pi, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref

# Drop augref-style SV alt contigs and any short scaffolds, otherwise the
# augref plot is dominated by thousands of ~1 kb contigs each forced to
# occupy a 1 Mb gap on the x-axis.
df["PI"] = pd.to_numeric(df["PI"], errors="coerce")
df = df.dropna(subset=["PI"])

df = df[~df["CHROM"].astype(str).str.startswith("SV_")]
chrom_span_all = df.groupby("CHROM")["BIN_START"].agg(lambda s: s.max() - s.min())
keep_chroms    = chrom_span_all[chrom_span_all >= 1_000_000].index.tolist()
df             = df[df["CHROM"].isin(keep_chroms)]

if df.empty:
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.text(0.5, 0.5, f"No qualifying chromosomes for π — {ref_name}",
            ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)

chroms      = sorted(df["CHROM"].unique())
palette_arr = np.array(sns.color_palette("tab20", len(chroms)))

chrom_grp  = df.groupby("CHROM")["BIN_START"]
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

df["x"] = df["CHROM"].map(offsets).astype(float) + df["BIN_START"].astype(float)

x_range  = df["x"].max() - df["x"].min()
bin_size = max(1.0, x_range / 4800.0)
df["_bin"] = (df["x"] / bin_size).astype(np.int64)
df = df.loc[df.groupby(["CHROM", "_bin"])["PI"].idxmax()].drop(columns="_bin")

chrom_to_idx = {c: i for i, c in enumerate(chroms)}
colors = palette_arr[df["CHROM"].map(chrom_to_idx).values % len(palette_arr)]

fig, ax = plt.subplots(figsize=(16, 4))
ax.scatter(df["x"].values, df["PI"].values,
           c=colors, s=4, alpha=0.5, linewidths=0, rasterized=True)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("π")
ax.set_title("Nucleotide diversity — " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
