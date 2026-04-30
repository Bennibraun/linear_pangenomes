import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

sns.set_theme(style="ticks", font_scale=0.9)

df      = pd.read_csv(snakemake.input.metrics, sep="\t")
out_pdf = snakemake.output.pdf
out_png = snakemake.output.png

# ── numeric coercion (cactus emits "NA" for depth / mapq) ─────────────────────
df["mapping_rate_pct"] = pd.to_numeric(df["mapping_rate"], errors="coerce") * 100
df["mean_mapq"]        = pd.to_numeric(df["mean_mapq"],    errors="coerce")
df["mean_depth"]       = pd.to_numeric(df["mean_depth"],   errors="coerce")

# ── reference-type order: alpha, cactus last ───────────────────────────────────
all_refs    = df["alignment_type"].unique().tolist()
cactus_refs = [r for r in all_refs if "cactus" in r.lower()]
other_refs  = sorted(r for r in all_refs if "cactus" not in r.lower())
ref_order   = other_refs + cactus_refs

samples = sorted(df["sample"].unique())
palette = dict(zip(samples, sns.color_palette("tab10", len(samples))))
rng     = np.random.default_rng(0)


def _strip_panel(ax, col, ylabel, ylabel_extra=""):
    """One dot per (sample, ref), horizontal mean line, light jitter."""
    for xi, ref in enumerate(ref_order):
        sub = df.loc[df["alignment_type"] == ref, ["sample", col]].dropna(subset=[col])
        if sub.empty:
            continue
        jitter = rng.uniform(-0.18, 0.18, len(sub))
        for j, (_, row) in zip(jitter, sub.iterrows()):
            ax.scatter(xi + j, row[col],
                       color=palette[row["sample"]], s=28,
                       alpha=0.85, linewidths=0, zorder=3, clip_on=False)
        ax.hlines(sub[col].mean(), xi - 0.32, xi + 0.32,
                  colors="#222222", lw=1.8, zorder=4)

    ax.set_xticks(range(len(ref_order)))
    ax.set_xticklabels(ref_order, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel(ylabel + ylabel_extra, fontsize=9)
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    ax.tick_params(axis="both", which="major", labelsize=8)
    sns.despine(ax=ax, trim=True)


# ── figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(11, 4), constrained_layout=True)

_strip_panel(axes[0], "mapping_rate_pct", "Mapping rate", " (%)")
lo = max(0.0, df["mapping_rate_pct"].dropna().min() - 3.0)
axes[0].set_ylim(bottom=lo, top=min(100, df["mapping_rate_pct"].dropna().max() + 1.5))
axes[0].set_title("Mapping rate", fontsize=10, fontweight="bold", pad=6)

_strip_panel(axes[1], "mean_mapq", "Mean MAPQ")
axes[1].set_title("Mean mapping quality", fontsize=10, fontweight="bold", pad=6)

_strip_panel(axes[2], "mean_depth", "Mean depth", " (×)")
axes[2].set_title("Mean coverage depth", fontsize=10, fontweight="bold", pad=6)

# ── shared sample legend below the panels ─────────────────────────────────────
handles = [
    plt.Line2D([0], [0], marker="o", color="w",
               markerfacecolor=c, markersize=6, label=s)
    for s, c in palette.items()
]
fig.legend(
    handles=handles,
    title="Sample", title_fontsize=8,
    fontsize=7,
    loc="lower center",
    ncol=min(len(samples), 8),
    bbox_to_anchor=(0.5, -0.10),
    frameon=False,
    handletextpad=0.3,
    columnspacing=0.8,
)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close(fig)
