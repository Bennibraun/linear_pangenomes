"""Per-population nucleotide diversity (π) summary — accessory sequence INCLUDED.

Input: long-format per-population windowed-π table (CHROM, BIN_START, BIN_END,
N_VARIANTS, PI, POP). Accessory (SV_*) contigs are the whole point of the
augmented reference, so they are kept. Panels:

  1. Per-population violin/box of π (all windows, core + accessory).
  2. Genome-wide mean π per population (bar).
  3. Core vs accessory π distribution (only when the ref has accessory sequence),
     with a Mann-Whitney p and Cliff's delta — the augref-specific result.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def region_class(chrom_series):
    is_acc = chrom_series.astype(str).str.startswith(("SV_", "UNMAP_"))
    return np.where(is_acc, "accessory", "core")


def has_accessory(df, chrom_col="CHROM"):
    return bool(df[chrom_col].astype(str).str.startswith(("SV_", "UNMAP_")).any())


def core_accessory_stats(values_core, values_acc):
    """(mannwhitneyu_p, cliffs_delta_acc_vs_core), or (nan, nan) if a side is empty."""
    core = np.asarray(values_core)
    acc = np.asarray(values_acc)
    if len(core) == 0 or len(acc) == 0:
        return np.nan, np.nan
    u, p = mannwhitneyu(acc, core, alternative="two-sided")
    delta = 2.0 * u / (len(acc) * len(core)) - 1.0
    return p, delta

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
if df.empty:
    _empty(f"No π windows — {ref_name}")

df["region"] = region_class(df["CHROM"])
accessory = has_accessory(df)

med = df.groupby("POP")["PI"].median().sort_values(ascending=False)
order = med.index.tolist()
colors = sns.color_palette("husl", len(order))
palette = dict(zip(order, colors))

n_panels = 3 if accessory else 2
fig, axes = plt.subplots(
    1, n_panels, figsize=(max(10, 1.2 * len(order) + 6) + (4 if accessory else 0), 5),
    gridspec_kw={"width_ratios": [2, 1] + ([1.2] if accessory else [])},
)
ax_v, ax_b = axes[0], axes[1]

# Panel 1: per-window π per population (core + accessory together).
sns.violinplot(data=df, x="POP", y="PI", order=order,
               palette=colors, cut=0, inner="box", linewidth=0.8, ax=ax_v)
ax_v.axhline(df["PI"].median(), color="grey", lw=0.8, linestyle="--",
             label="cohort median")
ax_v.set_xlabel("")
ax_v.set_ylabel("π (per window)")
ax_v.set_xticklabels(ax_v.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax_v.legend(frameon=False, fontsize=8)

# Panel 2: genome-wide mean π per population.
mean_pi = df.groupby("POP")["PI"].mean().reindex(order)
ax_b.barh(range(len(order)), mean_pi.values,
          color=[palette[p] for p in order], edgecolor="white")
ax_b.set_yticks(range(len(order)))
ax_b.set_yticklabels(order, fontsize=8)
ax_b.invert_yaxis()
ax_b.set_xlabel("mean π")
for i, v in enumerate(mean_pi.values):
    ax_b.text(v, i, f" {v:.4f}", va="center", fontsize=7)

# Panel 3: core vs accessory (augref only).
if accessory:
    ax_a = axes[2]
    core = df.loc[df["region"] == "core", "PI"]
    acc = df.loc[df["region"] == "accessory", "PI"]
    sns.violinplot(data=df, x="region", y="PI", order=["core", "accessory"],
                   palette=["#4C72B0", "#C44E52"], cut=0, inner="box",
                   linewidth=0.8, ax=ax_a)
    p, delta = core_accessory_stats(core, acc)
    p_txt = "p<1e-300" if p == 0 else f"p={p:.1e}"
    ax_a.set_title(f"core vs accessory\n{p_txt}, δ={delta:.2f}", fontsize=9)
    ax_a.set_xlabel("")
    ax_a.set_ylabel("π (per window)")

fig.suptitle(f"Nucleotide diversity by population — {ref_name}",
             fontweight="bold")
sns.despine(fig=fig)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
