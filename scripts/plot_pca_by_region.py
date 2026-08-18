"""Core-only vs accessory-only PCA, side by side, for augref. Directly tests
whether the three populations delineate clearly on core sites but blur
together on accessory sites — the expected signature if the accessory
sequence is capturing variation shared broadly across the cohort rather than
population-specific variation (low FST in accessory)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

cov_core   = snakemake.input.cov_core
cov_acc    = snakemake.input.cov_accessory
samples_core = snakemake.input.samples_core
samples_acc  = snakemake.input.samples_accessory
out_pdf    = snakemake.output.pdf
out_png    = snakemake.output.png
sample_pop_str = snakemake.params.sample_pop

is_pruned = open(snakemake.input.pruned_flag).read().strip() == "true"

sample_to_pop = {}
for pair in sample_pop_str.split("|"):
    k, v = pair.split(",", 1)
    sample_to_pop[k] = v

pop_names = sorted(set(sample_to_pop.values()))
palette   = sns.color_palette("tab10", len(pop_names))
pop_color = {p: palette[i] for i, p in enumerate(pop_names)}


def pc12(cov_file, sample_file):
    cov = np.loadtxt(cov_file)
    samples = open(sample_file).read().splitlines()
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    pve = eigenvalues / eigenvalues.sum() * 100
    return samples, eigenvectors, pve


fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (title, cov_file, sample_file) in zip(
    axes,
    [("core", cov_core, samples_core), ("accessory", cov_acc, samples_acc)],
):
    samples, eigenvectors, pve = pc12(cov_file, sample_file)
    populations = [sample_to_pop.get(s, "unknown") for s in samples]
    for i, (s, pop) in enumerate(zip(samples, populations)):
        ax.scatter(eigenvectors[i, 0], eigenvectors[i, 1],
                   color=pop_color[pop], s=60, edgecolors="white", linewidths=0.4)
        ax.annotate(s, (eigenvectors[i, 0], eigenvectors[i, 1]),
                    fontsize=5, ha="left", va="bottom", alpha=0.7)
    ax.set_xlabel("PC1 (" + format(pve[0], ".1f") + "% var)")
    ax.set_ylabel("PC2 (" + format(pve[1], ".1f") + "% var)")
    ax.axhline(0, color="grey", lw=0.5, linestyle="--")
    ax.axvline(0, color="grey", lw=0.5, linestyle="--")
    ax.set_title(title, fontsize=11)
    sns.despine(ax=ax)

title = "PCA: core vs accessory sites — augref"
if not is_pruned:
    title += "  [UNPRUNED: N<50 samples, LD pruning skipped]"
fig.suptitle(title, fontsize=13, fontweight="bold",
             color="crimson" if not is_pruned else "black")
legend_handles = [mpatches.Patch(color=pop_color[p], label=p) for p in pop_names]
fig.legend(handles=legend_handles, title="Population",
           loc="lower center", ncol=len(pop_names),
           bbox_to_anchor=(0.5, -0.04), frameon=False)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
