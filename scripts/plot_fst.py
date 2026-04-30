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

df[fst_col] = pd.to_numeric(df[fst_col], errors="coerce")
df = df.dropna(subset=[fst_col, pos_col])

# Drop augref-style alternate contigs (SV_chrom_pos_sample_insN) so they don't
# stretch the genome-wide x-axis. These add ~1 Mb of empty gap each and there
# can be thousands of them. Conventional unplaced scaffolds (which on most
# refs are also short) are dropped by the size filter below.
df = df[~df["CHROM"].astype(str).str.startswith("SV_")]

# Drop tiny contigs (< 1 Mb of span) — keeps only true chromosomes and
# substantial scaffolds. Span is the spread of variant positions on the
# contig in this VCF, which is a reasonable proxy for "chromosome-scale".
chrom_span_all = df.groupby("CHROM")[pos_col].agg(lambda s: s.max() - s.min())
keep_chroms    = chrom_span_all[chrom_span_all >= 1_000_000].index.tolist()
df             = df[df["CHROM"].isin(keep_chroms)]

if df.empty:
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.text(0.5, 0.5, f"No qualifying chromosomes for {ref_name}: {pop1} vs {pop2}",
            ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {c: palette[i % len(palette)] for i, c in enumerate(chroms)}

# Compute per-chrom offsets via groupby (vectorized, single pass).
chrom_min  = df.groupby("CHROM")[pos_col].min()
chrom_max  = df.groupby("CHROM")[pos_col].max()
chrom_span = (chrom_max - chrom_min).reindex(chroms)

offsets, centers = {}, {}
offset = 0
for chrom in chroms:
    span = chrom_span[chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + span / 2
    offset += span + 1_000_000

df["x"] = df["CHROM"].map(offsets).astype(float) + df[pos_col].astype(float)
# Clip negative FST to 0 for visualization (negatives indicate within-pop
# diversity exceeds between-pop diversity due to sampling noise; they were
# previously dropped, which hid all low-differentiation sites and made the
# plot look uniformly differentiated).
df["fst_plot"] = df[fst_col].clip(lower=0)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    ax.scatter(sub["x"], sub["fst_plot"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("FST")
ax.set_ylim(bottom=0)
ax.set_title("FST — " + ref_name + ": " + pop1 + " vs " + pop2, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
