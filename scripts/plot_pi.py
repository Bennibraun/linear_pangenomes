"""Per-population nucleotide diversity (π) summary.

Input is the long-format per-population windowed-π table (CHROM, BIN_START,
BIN_END, N_VARIANTS, PI, POP). We summarise diversity *across populations*
rather than along the genome (a manhattan mixes populations and buries the
comparison). Two panels:

  1. Violin + box of per-window π for each population, sorted by median, with a
     genome-wide median line. Answers "which populations are diverse vs
     bottlenecked."
  2. A bar of the genome-wide mean π per population (the headline number).

SV_* accessory contigs are excluded so the per-population summary reflects the
core genome (consistent across references).
"""
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


def _empty(msg):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)


if df.empty or "POP" not in df.columns:
    _empty(f"No per-population π data — {ref_name}")

df["PI"] = pd.to_numeric(df["PI"], errors="coerce")
df = df.dropna(subset=["PI"])
df = df[~df["CHROM"].astype(str).str.startswith("SV_")]
if df.empty:
    _empty(f"No core-genome π windows — {ref_name}")

# Order populations by median π (most diverse first).
med = df.groupby("POP")["PI"].median().sort_values(ascending=False)
order = med.index.tolist()
palette = dict(zip(order, sns.color_palette("husl", len(order))))

fig, (ax_v, ax_b) = plt.subplots(
    1, 2, figsize=(max(10, 1.2 * len(order) + 6), 5),
    gridspec_kw={"width_ratios": [2, 1]},
)

# Panel 1: per-window π distribution per population.
sns.violinplot(data=df, x="POP", y="PI", order=order, hue="POP",
               palette=palette, legend=False, cut=0, inner="box",
               linewidth=0.8, ax=ax_v)
ax_v.axhline(df["PI"].median(), color="grey", lw=0.8, linestyle="--",
             label="cohort median")
ax_v.set_xlabel("")
ax_v.set_ylabel("π (per window)")
ax_v.set_xticklabels(ax_v.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax_v.legend(frameon=False, fontsize=8)

# Panel 2: genome-wide mean π per population (headline number).
mean_pi = df.groupby("POP")["PI"].mean().reindex(order)
ax_b.barh(range(len(order)), mean_pi.values,
          color=[palette[p] for p in order], edgecolor="white")
ax_b.set_yticks(range(len(order)))
ax_b.set_yticklabels(order, fontsize=8)
ax_b.invert_yaxis()
ax_b.set_xlabel("mean π")
for i, v in enumerate(mean_pi.values):
    ax_b.text(v, i, f" {v:.4f}", va="center", fontsize=7)

fig.suptitle(f"Nucleotide diversity by population — {ref_name}",
             fontweight="bold")
sns.despine(fig=fig)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
