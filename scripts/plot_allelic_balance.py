import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref
samples  = snakemake.params.samples
paths    = list(snakemake.input)

n     = len(samples)
ncols = min(4, n)
nrows = (n + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
fig.suptitle("Allelic balance \u2014 " + ref_name, fontsize=11, fontweight="bold")

for idx, (sample, path) in enumerate(zip(samples, paths)):
    ax   = axes[idx // ncols][idx % ncols]
    # File has two sections: summary rows first, then bin rows starting after "bin_low" header
    lines = open(path).read().splitlines()
    bin_start = next(i for i, l in enumerate(lines) if l.startswith("bin_low"))
    bins = pd.read_csv(path, sep="\t", skiprows=bin_start, nrows=len(lines) - bin_start - 1)
    bins["bin_low"] = bins["bin_low"].astype(float)
    bins["count"]   = bins["count"].astype(int)
    labels = [format(r["bin_low"], ".1f") for _, r in bins.iterrows()]
    ax.bar(range(len(bins)), bins["count"],
           color=sns.color_palette("muted")[0], edgecolor="white")
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_title(sample, fontsize=8)
    ax.set_xlabel("Ref ratio", fontsize=7)
    ax.set_ylabel("Sites", fontsize=7)
    sns.despine(ax=ax)

for idx in range(n, nrows * ncols):
    axes[idx // ncols][idx % ncols].set_visible(False)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
