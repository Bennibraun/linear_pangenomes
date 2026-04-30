import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="ticks", font_scale=0.9)

df      = pd.read_csv(snakemake.input.metrics, sep="\t")
out_pdf = snakemake.output.pdf
out_png = snakemake.output.png

df["mapping_rate_pct"] = pd.to_numeric(df["mapping_rate"], errors="coerce") * 100

pivot = (
    df.pivot(index="short_sample", columns="long_sample", values="mapping_rate_pct")
    .sort_index(axis=0)
    .sort_index(axis=1)
)

n_rows, n_cols = pivot.shape
figw = max(4.5, 0.75 * n_cols + 2.0)
figh = max(3.5, 0.60 * n_rows + 1.5)

fig, ax = plt.subplots(figsize=(figw, figh), constrained_layout=True)

lo = max(0.0, pivot.min().min() - 2.0)

# Annotate cells only when the grid is small enough to read
annot = n_rows * n_cols <= 400

sns.heatmap(
    pivot,
    ax=ax,
    cmap="YlOrRd",
    vmin=lo,
    vmax=100,
    annot=annot,
    fmt=".1f",
    annot_kws={"size": 7},
    linewidths=0.35,
    linecolor="white",
    cbar_kws={"label": "Mapping rate (%)", "shrink": 0.75, "pad": 0.02},
)

ax.set_xlabel("Assembly (long-read sample)", fontsize=9, labelpad=6)
ax.set_ylabel("Short-read sample", fontsize=9, labelpad=6)
ax.set_title(
    "Short-read mapping rate to individual assemblies",
    fontsize=10, fontweight="bold", pad=8,
)
ax.tick_params(axis="x", rotation=35, labelsize=8)
ax.tick_params(axis="y", rotation=0,  labelsize=8)

# colour-bar label size
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.set_label("Mapping rate (%)", fontsize=8)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close(fig)
