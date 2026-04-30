import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.tajima, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref

# vcftools writes "nan" for windows with too few SNPs; treat those as missing.
df["TajimaD"] = pd.to_numeric(df["TajimaD"], errors="coerce")
df = df.dropna(subset=["TajimaD"])

# Drop augref SV alt contigs and tiny scaffolds (same approach as plot_pi.py).
df = df[~df["CHROM"].astype(str).str.startswith("SV_")]
chrom_span_all = df.groupby("CHROM")["BIN_START"].agg(lambda s: s.max() - s.min())
keep_chroms    = chrom_span_all[chrom_span_all >= 1_000_000].index.tolist()
df             = df[df["CHROM"].isin(keep_chroms)]

if df.empty:
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.text(0.5, 0.5, f"No qualifying chromosomes for Tajima's D — {ref_name}",
            ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {c: palette[i % len(palette)] for i, c in enumerate(chroms)}

chrom_min  = df.groupby("CHROM")["BIN_START"].min()
chrom_max  = df.groupby("CHROM")["BIN_START"].max()
chrom_span = (chrom_max - chrom_min).reindex(chroms)

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    span = chrom_span[chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + span / 2
    offset += span + 1_000_000

df["x"] = df["CHROM"].map(offsets).astype(float) + df["BIN_START"].astype(float)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    ax.scatter(sub["x"], sub["TajimaD"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.axhline(0, color="grey", lw=0.5, linestyle="--")
ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Tajima's D")
ax.set_title("Tajima's D — " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
