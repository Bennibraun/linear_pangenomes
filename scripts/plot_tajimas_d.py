"""Per-population Tajima's D summary.

Input is the long-format per-population Tajima's D table (CHROM, BIN_START,
N_SNPS, TajimaD, POP). Tajima's D is a demographic/selection summary read *per
population*, so we compare distributions across populations rather than plotting
position. Two panels:

  1. Violin + box of per-window Tajima's D for each population, with the D = 0
     neutral line. D < 0 → excess rare alleles (expansion / purifying or positive
     selection); D > 0 → excess intermediate-frequency alleles (contraction /
     balancing selection / structure).
  2. Genome-wide mean Tajima's D per population (bar), the headline number.

SV_* accessory contigs are excluded so the summary reflects the core genome.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.tajima, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref


def _empty(msg):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)


if df.empty or "POP" not in df.columns:
    _empty(f"No per-population Tajima's D data — {ref_name}")

df["TajimaD"] = pd.to_numeric(df["TajimaD"], errors="coerce")
df = df.dropna(subset=["TajimaD"])
df = df[~df["CHROM"].astype(str).str.startswith("SV_")]
if df.empty:
    _empty(f"No core-genome Tajima's D windows — {ref_name}")

# Order populations by median Tajima's D (most positive first).
med = df.groupby("POP")["TajimaD"].median().sort_values(ascending=False)
order = med.index.tolist()
colors = sns.color_palette("husl", len(order))
palette = dict(zip(order, colors))  # for the bar panel

fig, (ax_v, ax_b) = plt.subplots(
    1, 2, figsize=(max(10, 1.2 * len(order) + 6), 5),
    gridspec_kw={"width_ratios": [2, 1]},
)

# Panel 1: per-window Tajima's D distribution per population.
# palette as an ordered list (not hue=) so this works across seaborn versions.
sns.violinplot(data=df, x="POP", y="TajimaD", order=order,
               palette=colors, cut=0, inner="box",
               linewidth=0.8, ax=ax_v)
ax_v.axhline(0, color="black", lw=0.9, linestyle="--", label="D = 0 (neutral)")
ax_v.set_xlabel("")
ax_v.set_ylabel("Tajima's D (per window)")
ax_v.set_xticklabels(ax_v.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax_v.legend(frameon=False, fontsize=8)

# Panel 2: genome-wide mean Tajima's D per population.
mean_d = df.groupby("POP")["TajimaD"].mean().reindex(order)
ax_b.barh(range(len(order)), mean_d.values,
          color=[palette[p] for p in order], edgecolor="white")
ax_b.axvline(0, color="black", lw=0.8)
ax_b.set_yticks(range(len(order)))
ax_b.set_yticklabels(order, fontsize=8)
ax_b.invert_yaxis()
ax_b.set_xlabel("mean Tajima's D")
for i, v in enumerate(mean_d.values):
    ha = "left" if v >= 0 else "right"
    ax_b.text(v, i, f" {v:.3f} ", va="center", ha=ha, fontsize=7)

fig.suptitle(f"Tajima's D by population — {ref_name}", fontweight="bold")
sns.despine(fig=fig)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
