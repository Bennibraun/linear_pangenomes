"""Depth/MAPQ vs Tajima's D, core vs accessory (augref) — the artifact-check
plot. Joins region_qc_summary.tsv (per-sample mean depth/MAPQ/missingness,
split core vs accessory) against the windowed Tajima's D table (aggregated to
the same core/accessory split via the SV_* CHROM prefix).

If accessory Tajima's D is elevated because depth/MAPQ there is genuinely
lower (samtools undercalling rare variants in poorly-mapped regions), the
per-sample points should show a negative depth/MAPQ vs Tajima's D relationship
concentrated in the accessory region. If accessory Tajima's D is elevated
despite comparable depth/MAPQ to core, the artifact explanation is much less
likely and the shared-variation interpretation gains support."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

qc      = pd.read_csv(snakemake.input.region_qc, sep="\t")
tajima  = pd.read_csv(snakemake.input.tajima, sep="\t")
out_pdf = snakemake.output.pdf
out_png = snakemake.output.png

tajima = tajima.dropna(subset=["TajimaD"])
tajima["region"] = tajima["CHROM"].astype(str).str.startswith("SV_").map(
    {True: "accessory", False: "core"}
)
tajd_by_region = (
    tajima.groupby(["POP", "region"])["TajimaD"]
    .median()
    .reset_index()
    .rename(columns={"TajimaD": "median_tajd"})
)

qc_mean = (
    qc.groupby("region")[["mean_depth", "mean_mapq", "missing_rate"]]
    .mean()
    .reset_index()
)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Panel 1: per-sample depth by region (paired lines show whether every
# sample individually has lower accessory depth, not just the mean).
ax = axes[0]
for sample, sub in qc.groupby("sample"):
    sub = sub.set_index("region").reindex(["core", "accessory"])
    ax.plot(["core", "accessory"], sub["mean_depth"], color="grey", alpha=0.4, lw=0.8, marker="o", ms=3)
ax.set_ylabel("Mean depth")
ax.set_title("Per-sample depth by region", fontsize=10)
sns.despine(ax=ax)

# Panel 2: per-sample MAPQ by region.
ax = axes[1]
for sample, sub in qc.groupby("sample"):
    sub = sub.set_index("region").reindex(["core", "accessory"])
    ax.plot(["core", "accessory"], sub["mean_mapq"], color="grey", alpha=0.4, lw=0.8, marker="o", ms=3)
ax.set_ylabel("Mean MAPQ")
ax.set_title("Per-sample MAPQ by region", fontsize=10)
sns.despine(ax=ax)

# Panel 3: region-level median Tajima's D (per population) vs region-level
# mean depth — the direct artifact-vs-biology comparison. If low-depth
# regions are also the high-Tajima's-D regions across populations, that's
# the calling-artifact signature.
ax = axes[2]
merged = tajd_by_region.merge(qc_mean, on="region")
palette = {"core": "#4C72B0", "accessory": "#C44E52"}
for region, sub in merged.groupby("region"):
    ax.scatter(sub["mean_depth"], sub["median_tajd"], color=palette[region],
               label=region, s=70, edgecolors="white", linewidths=0.5)
    for _, row in sub.iterrows():
        ax.annotate(row["POP"], (row["mean_depth"], row["median_tajd"]),
                    fontsize=6, ha="left", va="bottom", alpha=0.8)
ax.set_xlabel("Mean depth (region)")
ax.set_ylabel("Median Tajima's D (per population)")
ax.set_title("Depth vs Tajima's D by region", fontsize=10)
ax.axhline(0, color="black", lw=0.6, linestyle="--")
ax.legend(frameon=False, fontsize=8)
sns.despine(ax=ax)

fig.suptitle("Core vs accessory QC diagnostic — augref", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
