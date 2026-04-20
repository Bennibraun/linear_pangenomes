import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

cov_file    = snakemake.input.cov
sample_file = snakemake.input.samples
out_pdf     = snakemake.output.pdf
out_png     = snakemake.output.png
ref_name    = snakemake.params.ref
sample_pop_str = snakemake.params.sample_pop

sample_to_pop = {}
for pair in sample_pop_str.split("|"):
    k, v = pair.split(",", 1)
    sample_to_pop[k] = v

cov     = np.loadtxt(cov_file)
samples = open(sample_file).read().splitlines()

eigenvalues, eigenvectors = np.linalg.eigh(cov)
idx          = np.argsort(eigenvalues)[::-1]
eigenvalues  = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]
pve          = eigenvalues / eigenvalues.sum() * 100

populations = [sample_to_pop.get(s, "unknown") for s in samples]
pop_names   = sorted(set(populations))
palette     = sns.color_palette("tab10", len(pop_names))
pop_color   = {p: palette[i] for i, p in enumerate(pop_names)}

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("PCA \u2014 " + ref_name, fontsize=13, fontweight="bold")

for ax, (px, py) in zip(axes, [(0, 1), (0, 2)]):
    for i, (s, pop) in enumerate(zip(samples, populations)):
        ax.scatter(eigenvectors[i, px], eigenvectors[i, py],
                   color=pop_color[pop], s=60, edgecolors="white", linewidths=0.4)
        ax.annotate(s, (eigenvectors[i, px], eigenvectors[i, py]),
                    fontsize=5, ha="left", va="bottom", alpha=0.7)
    ax.set_xlabel("PC" + str(px + 1) + " (" + format(pve[px], ".1f") + "% var)")
    ax.set_ylabel("PC" + str(py + 1) + " (" + format(pve[py], ".1f") + "% var)")
    ax.axhline(0, color="grey", lw=0.5, linestyle="--")
    ax.axvline(0, color="grey", lw=0.5, linestyle="--")
    sns.despine(ax=ax)

legend_handles = [mpatches.Patch(color=pop_color[p], label=p) for p in pop_names]
fig.legend(handles=legend_handles, title="Population",
           loc="lower center", ncol=len(pop_names),
           bbox_to_anchor=(0.5, -0.04), frameon=False)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
