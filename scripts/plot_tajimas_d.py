"""Per-population Tajima's D summary — accessory sequence INCLUDED.

Input: long-format per-population Tajima's D table (CHROM, BIN_START, N_SNPS,
TajimaD, POP). Accessory (SV_*) contigs are kept (they are the reason augref
exists). Panels:

  1. Per-population violin/box of Tajima's D (core + accessory), D = 0 line.
  2. Genome-wide mean Tajima's D per population (bar).
  3. Core vs accessory Tajima's D (only when the ref has accessory sequence),
     with Mann-Whitney p and Cliff's delta.
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
if df.empty:
    _empty(f"No Tajima's D windows — {ref_name}")

df["region"] = region_class(df["CHROM"])
accessory = has_accessory(df)

med = df.groupby("POP")["TajimaD"].median().sort_values(ascending=False)
order = med.index.tolist()
colors = sns.color_palette("husl", len(order))
palette = dict(zip(order, colors))

n_panels = 3 if accessory else 2
fig, axes = plt.subplots(
    1, n_panels, figsize=(max(10, 1.2 * len(order) + 6) + (4 if accessory else 0), 5),
    gridspec_kw={"width_ratios": [2, 1] + ([1.2] if accessory else [])},
)
ax_v, ax_b = axes[0], axes[1]

# Panel 1: per-window Tajima's D per population (core + accessory).
sns.violinplot(data=df, x="POP", y="TajimaD", order=order,
               palette=colors, cut=0, inner="box", linewidth=0.8, ax=ax_v)
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

# Panel 3: core vs accessory (augref only).
if accessory:
    ax_a = axes[2]
    core = df.loc[df["region"] == "core", "TajimaD"]
    acc = df.loc[df["region"] == "accessory", "TajimaD"]
    sns.violinplot(data=df, x="region", y="TajimaD", order=["core", "accessory"],
                   palette=["#4C72B0", "#C44E52"], cut=0, inner="box",
                   linewidth=0.8, ax=ax_a)
    ax_a.axhline(0, color="black", lw=0.8, linestyle="--")
    p, delta = core_accessory_stats(core, acc)
    p_txt = "p<1e-300" if p == 0 else f"p={p:.1e}"
    ax_a.set_title(f"core vs accessory\n{p_txt}, δ={delta:.2f}", fontsize=9)
    ax_a.set_xlabel("")
    ax_a.set_ylabel("Tajima's D (per window)")

fig.suptitle(f"Tajima's D by population — {ref_name}", fontweight="bold")
sns.despine(fig=fig)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
