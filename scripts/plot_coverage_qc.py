"""Per-sample coverage QC plot.

Reads the cohort alignment_metrics.tsv (one row per sample × ref) and emits
a per-sample mean-depth bar chart faceted by reference. Samples with
depth < min_depth_threshold are flagged in red so degraded samples surface
at a glance.

mc_graph mean_depth is reported as 'NA' because it comes from vg stats on
the GAM, so we drop it from the depth plot — only BAM-backed refs are
plotted.
"""

from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

metrics_path = snakemake.input.metrics
out_pdf = snakemake.output.pdf
out_png = snakemake.output.png
min_depth = float(snakemake.params.min_depth)

df = pd.read_csv(metrics_path, sep="\t")

# mean_depth is "NA" for the cactus / mc_graph rows (GAM-derived); drop those.
df = df[df["mean_depth"].astype(str) != "NA"].copy()
df["mean_depth"] = pd.to_numeric(df["mean_depth"], errors="coerce")
df = df.dropna(subset=["mean_depth"])

if df.empty:
    # Emit empty placeholder rather than crashing — keeps `snakemake --touch`
    # workflows happy.
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, "no coverage data available",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    raise SystemExit(0)

refs = sorted(df["alignment_type"].unique())
n_refs = len(refs)
fig, axes = plt.subplots(n_refs, 1, figsize=(10, 2.5 * n_refs), sharex=True)
if n_refs == 1:
    axes = [axes]

for ax, ref in zip(axes, refs):
    sub = df[df["alignment_type"] == ref].sort_values("sample")
    colors = ["#c0392b" if d < min_depth else "#2c7fb8" for d in sub["mean_depth"]]
    ax.bar(sub["sample"], sub["mean_depth"], color=colors, edgecolor="white")
    ax.axhline(min_depth, color="black", linestyle="--", linewidth=0.8,
               label=f"min depth threshold = {min_depth:g}×")
    ax.set_title(f"Mean depth — {ref}", fontweight="bold")
    ax.set_ylabel("× coverage")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    sns.despine(ax=ax)

axes[-1].set_xlabel("sample")
plt.setp(axes[-1].get_xticklabels(), rotation=45, ha="right")
fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
