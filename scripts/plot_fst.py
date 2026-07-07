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

# Accessory (SV_*) contigs have no genome coordinate, so they can't be placed on
# the positional manhattan x-axis. Instead of discarding them, split them out and
# show their FST in a side panel — the accessory sequence is the reason augref
# exists and its differentiation should be visible.
is_acc = df["CHROM"].astype(str).str.startswith("SV_")
acc_fst = pd.to_numeric(df.loc[is_acc, fst_col], errors="coerce").clip(lower=0).dropna()
df = df[~is_acc]

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
palette_arr = np.array(sns.color_palette("tab20", len(chroms)))

# Compute per-chrom offsets (vectorized).
chrom_grp  = df.groupby("CHROM")[pos_col]
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

df["x"] = df["CHROM"].map(offsets).astype(float) + df[pos_col].astype(float)
df["fst_plot"] = df[fst_col].clip(lower=0)

# Thin to plotting resolution: keep the highest-FST point per pixel column so
# peaks are preserved exactly. At 16 in × 150 dpi ≈ 2400 horizontal pixels;
# using 2× oversampling (4800 bins) guarantees nothing is missed.
x_range  = df["x"].max() - df["x"].min()
bin_size = max(1.0, x_range / 4800.0)
df["_bin"] = (df["x"] / bin_size).astype(np.int64)
df = df.loc[df.groupby(["CHROM", "_bin"])["fst_plot"].idxmax()].drop(columns="_bin")

# Assign per-row colors from a pre-built array — single scatter call avoids
# the PathCollection-per-chromosome overhead that caused the 30+ min runtime.
chrom_to_idx = {c: i for i, c in enumerate(chroms)}
colors = palette_arr[df["CHROM"].map(chrom_to_idx).values % len(palette_arr)]

# Add a narrow accessory-FST side panel when the ref has accessory sequence.
if len(acc_fst):
    fig, (ax, ax_acc) = plt.subplots(
        1, 2, figsize=(17, 4), gridspec_kw={"width_ratios": [16, 1.2]},
        sharey=True,
    )
else:
    fig, ax = plt.subplots(figsize=(16, 4))

ax.scatter(df["x"].values, df["fst_plot"].values,
           c=colors, s=4, alpha=0.5, linewidths=0, rasterized=True)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("FST")
ax.set_ylim(bottom=0)
ax.set_title("FST — " + ref_name + ": " + pop1 + " vs " + pop2, fontweight="bold")
sns.despine(ax=ax)

if len(acc_fst):
    ax_acc.boxplot([acc_fst.values], widths=0.6, showfliers=False)
    ax_acc.scatter(np.random.normal(1, 0.06, len(acc_fst)), acc_fst.values,
                   s=4, alpha=0.4, color="#C44E52", linewidths=0)
    ax_acc.set_xticks([1])
    ax_acc.set_xticklabels(["accessory\n(SV_*)"], fontsize=7)
    ax_acc.set_title(f"n={len(acc_fst)}", fontsize=8)
    sns.despine(ax=ax_acc)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
