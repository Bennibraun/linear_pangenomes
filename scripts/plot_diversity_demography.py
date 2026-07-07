"""Demographic-quadrant plot: genome-wide π vs Tajima's D, one point per population.

This is the compact summary that ties diversity and the allele-frequency
spectrum together. Each population is a single point (mean π, mean Tajima's D):

  - low π,  D < 0  : bottleneck followed by expansion / strong purifying selection
  - high π, D < 0  : large, expanding population
  - low π,  D > 0  : contraction / recent bottleneck
  - high π, D > 0  : population structure or balancing selection

Both means are over the whole genome (core + accessory SV_* windows).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

pi_path     = snakemake.input.pi
tajima_path = snakemake.input.tajima
out_pdf     = snakemake.output.pdf
out_png     = snakemake.output.png
ref_name    = snakemake.params.ref


def _clean(df, col):
    # Accessory (SV_*) windows are kept — the whole genome, core + accessory.
    df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=[col])


def _empty(msg):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)


pi = pd.read_csv(pi_path, sep="\t")
td = pd.read_csv(tajima_path, sep="\t")
if pi.empty or td.empty or "POP" not in pi.columns or "POP" not in td.columns:
    _empty(f"No per-population diversity data — {ref_name}")

pi = _clean(pi, "PI")
td = _clean(td, "TajimaD")

pi_mean = pi.groupby("POP")["PI"].mean()
td_mean = td.groupby("POP")["TajimaD"].mean()
summary = pd.concat([pi_mean, td_mean], axis=1).dropna()
summary.columns = ["mean_pi", "mean_tajd"]
if summary.empty:
    _empty(f"No overlapping populations for π and Tajima's D — {ref_name}")

pops = summary.index.tolist()
palette = dict(zip(pops, sns.color_palette("husl", len(pops))))

fig, ax = plt.subplots(figsize=(8, 6))
ax.axhline(0, color="grey", lw=0.8, linestyle="--")
for pop, row in summary.iterrows():
    ax.scatter(row["mean_pi"], row["mean_tajd"],
               color=palette[pop], s=110, edgecolors="white", linewidths=0.8,
               zorder=3)
    ax.annotate(pop, (row["mean_pi"], row["mean_tajd"]),
                textcoords="offset points", xytext=(6, 4), fontsize=8)

ax.set_xlabel("mean π (nucleotide diversity)")
ax.set_ylabel("mean Tajima's D")
ax.set_title(f"Diversity vs demography by population — {ref_name}",
             fontweight="bold")
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
