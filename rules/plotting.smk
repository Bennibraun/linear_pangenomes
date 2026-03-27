# ============================================================================
# Plotting Rules
# ============================================================================
# All rules use shell: + Python heredoc so the conda env is properly activated.
# Snakemake interpolates {input.*}, {output.*}, {params.*} before bash runs;
# Python dict literals and any other braces in the heredoc body use {{ }}.
# Points are colored/labeled by the 'grouping' column from the reads manifest
# via SAMPLE_TO_POP / POP_NAMES, passed as params.


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------
rule plot_pca:
    input:
        cov=PCANGSD_OUTDIR / "{ref}/pcangsd.cov",
        samples=PCANGSD_OUTDIR / "{ref}/samples.txt",
    output:
        pdf=PLOT_OUTDIR / "pca/{ref}_pca.pdf",
        png=PLOT_OUTDIR / "pca/{ref}_pca.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        sample_pop=lambda wildcards: "\n".join(
            s + "\t" + SAMPLE_TO_POP.get(s, "unknown")
            for s in SHORT_SAMPLES
        ),
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

cov_file    = "{input.cov}"
sample_file = "{input.samples}"
out_pdf     = "{output.pdf}"
out_png     = "{output.png}"
ref_name    = "{params.ref}"

# sample → population from manifest
sample_pop_raw = """{params.sample_pop}"""
sample_to_pop = {{}}
for line in sample_pop_raw.strip().splitlines():
    parts = line.split("\t")
    if len(parts) == 2:
        sample_to_pop[parts[0]] = parts[1]

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
pop_color   = {{p: palette[i] for i, p in enumerate(pop_names)}}

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
PYEOF
        """


# ---------------------------------------------------------------------------
# Selection scan (Manhattan plot of PCAngsd chi-squared scores)
# ---------------------------------------------------------------------------
rule plot_selection_scan:
    input:
        outliers=PCANGSD_OUTDIR / "{ref}/fst_outliers.tsv",
    output:
        pdf=PLOT_OUTDIR / "selection/{ref}_selection_scan.pdf",
        png=PLOT_OUTDIR / "selection/{ref}_selection_scan.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        fdr=0.05,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv("{input.outliers}", sep="\t")
out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"
fdr      = {params.fdr}

chroms        = sorted(df["chrom"].unique())
chrom_palette = sns.color_palette("tab20", len(chroms))
chrom_color   = {{c: chrom_palette[i % len(chrom_palette)] for i, c in enumerate(chroms)}}

offsets, centers = {{}}, {{}}
offset = 0
for chrom in chroms:
    sub = df[df["chrom"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub["pos"].max() - sub["pos"].min()) / 2
    offset += sub["pos"].max() - sub["pos"].min() + 1_000_000

df["x"] = df.apply(lambda r: offsets[r["chrom"]] + r["pos"], axis=1)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["chrom"] == chrom]
    ax.scatter(sub["x"], sub["chi2"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.6, linewidths=0)

outliers = df[df["outlier"]]
if len(outliers):
    ax.scatter(outliers["x"], outliers["chi2"],
               c="crimson", s=12, zorder=5, label="FDR < " + str(fdr))
    ax.axhline(outliers["chi2"].min(), color="crimson", lw=0.8, linestyle="--", alpha=0.7)
    ax.legend(frameon=False, fontsize=8)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Chi-squared score")
ax.set_title("Selection scan \u2014 " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
PYEOF
        """


# ---------------------------------------------------------------------------
# FST
# ---------------------------------------------------------------------------
rule plot_fst:
    input:
        fst=FST_OUTDIR / "fst/{ref}/{pop1}_vs_{pop2}.weir.fst",
    output:
        pdf=PLOT_OUTDIR / "fst/{ref}_{pop1}_vs_{pop2}_fst.pdf",
        png=PLOT_OUTDIR / "fst/{ref}_{pop1}_vs_{pop2}_fst.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        pop1=lambda wildcards: wildcards.pop1,
        pop2=lambda wildcards: wildcards.pop2,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv("{input.fst}", sep="\t")
out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"
pop1     = "{params.pop1}"
pop2     = "{params.pop2}"

is_windowed = "BIN_START" in df.columns
pos_col     = "BIN_START" if is_windowed else "POS"
fst_col     = "WEIGHTED_FST" if is_windowed else "WEIR_AND_COCKERHAM_FST"

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {{c: palette[i % len(palette)] for i, c in enumerate(chroms)}}

offsets, centers = {{}}, {{}}
offset = 0
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub[pos_col].max() - sub[pos_col].min()) / 2
    offset += sub[pos_col].max() - sub[pos_col].min() + 1_000_000

df["x"]  = df.apply(lambda r: offsets[r["CHROM"]] + r[pos_col], axis=1)
df_plot  = df[df[fst_col] >= 0]

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df_plot[df_plot["CHROM"] == chrom]
    ax.scatter(sub["x"], sub[fst_col],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("FST")
ax.set_ylim(bottom=0)
ax.set_title("FST \u2014 " + ref_name + ": " + pop1 + " vs " + pop2, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
PYEOF
        """


# ---------------------------------------------------------------------------
# Allele Frequency Spectrum
# ---------------------------------------------------------------------------
rule plot_afs:
    input:
        afs=FST_OUTDIR / "afs/{ref}.afs.tsv",
    output:
        pdf=PLOT_OUTDIR / "afs/{ref}_afs.pdf",
        png=PLOT_OUTDIR / "afs/{ref}_afs.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv("{input.afs}", sep="\t")
out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"

df = df[df["bin_low"].apply(lambda x: str(x).replace(".", "", 1).lstrip("-").isdigit())]
df["bin_low"]  = df["bin_low"].astype(float)
df["bin_high"] = df["bin_high"].astype(float)
df["count"]    = df["count"].astype(int)
labels = [format(r["bin_low"], ".2f") + "\u2013" + format(r["bin_high"], ".2f")
          for _, r in df.iterrows()]

fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(range(len(df)), df["count"],
       color=sns.color_palette("Blues_d", len(df)),
       edgecolor="white", linewidth=0.5)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.set_xlabel("Minor allele frequency")
ax.set_ylabel("Number of sites")
ax.set_title("Allele frequency spectrum \u2014 " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
PYEOF
        """


# ---------------------------------------------------------------------------
# Nucleotide diversity (π)
# ---------------------------------------------------------------------------
rule plot_pi:
    input:
        pi=PI_OUTDIR / "{ref}.windowed.pi",
    output:
        pdf=PLOT_OUTDIR / "pi/{ref}_pi.pdf",
        png=PLOT_OUTDIR / "pi/{ref}_pi.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv("{input.pi}", sep="\t")
out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"

chroms      = sorted(df["CHROM"].unique())
palette     = sns.color_palette("tab20", len(chroms))
chrom_color = {{c: palette[i % len(palette)] for i, c in enumerate(chroms)}}

offsets, centers = {{}}, {{}}
offset = 0
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    offsets[chrom] = offset
    centers[chrom] = offset + (sub["BIN_START"].max() - sub["BIN_START"].min()) / 2
    offset += sub["BIN_START"].max() - sub["BIN_START"].min() + 1_000_000

df["x"] = df.apply(lambda r: offsets[r["CHROM"]] + r["BIN_START"], axis=1)

fig, ax = plt.subplots(figsize=(16, 4))
for chrom in chroms:
    sub = df[df["CHROM"] == chrom]
    ax.scatter(sub["x"], sub["PI"],
               c=[chrom_color[chrom]] * len(sub), s=4, alpha=0.5, linewidths=0)

ax.set_xticks([centers[c] for c in chroms])
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=6)
ax.set_xlabel("Chromosome")
ax.set_ylabel("\u03c0")
ax.set_title("Nucleotide diversity \u2014 " + ref_name, fontweight="bold")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
PYEOF
        """


# ---------------------------------------------------------------------------
# Allelic balance
# ---------------------------------------------------------------------------
rule plot_allelic_balance:
    input:
        expand(
            AB_OUTDIR / "{{ref}}/{sample}.allelic_balance.tsv",
            sample=SHORT_SAMPLES,
        ),
    output:
        pdf=PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.pdf",
        png=PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        input_files=lambda wildcards, input: "\n".join(input),
        samples="\n".join(SHORT_SAMPLES),
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"
samples  = """{params.samples}""".strip().splitlines()
paths    = """{params.input_files}""".strip().splitlines()

n     = len(samples)
ncols = min(4, n)
nrows = (n + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
fig.suptitle("Allelic balance \u2014 " + ref_name, fontsize=11, fontweight="bold")

for idx, (sample, path) in enumerate(zip(samples, paths)):
    ax   = axes[idx // ncols][idx % ncols]
    full = pd.read_csv(path, sep="\t")
    bins = full[full["bin_low"].apply(
        lambda x: str(x).replace(".", "", 1).lstrip("-").isdigit()
    )].copy()
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
PYEOF
        """


# ---------------------------------------------------------------------------
# Runs of homozygosity (F_ROH)
# ---------------------------------------------------------------------------
rule plot_roh:
    input:
        summary=ROH_OUTDIR / "f_roh_summary.tsv",
    output:
        pdf=PLOT_OUTDIR / "roh/{ref}_f_roh.pdf",
        png=PLOT_OUTDIR / "roh/{ref}_f_roh.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        sample_pop=lambda wildcards: "\n".join(
            s + "\t" + SAMPLE_TO_POP.get(s, "unknown")
            for s in SHORT_SAMPLES
        ),
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pdf})
        python - << 'PYEOF'
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

df       = pd.read_csv("{input.summary}", sep="\t")
out_pdf  = "{output.pdf}"
out_png  = "{output.png}"
ref_name = "{params.ref}"

sample_pop_raw = """{params.sample_pop}"""
sample_to_pop = {{}}
for line in sample_pop_raw.strip().splitlines():
    parts = line.split("\t")
    if len(parts) == 2:
        sample_to_pop[parts[0]] = parts[1]

df = df[df["ref"] == ref_name].copy()
df["population"] = df["sample"].map(sample_to_pop).fillna("unknown")
df = df.sort_values(["population", "sample"])

pop_names = sorted(df["population"].unique())
palette   = sns.color_palette("tab10", len(pop_names))
pop_color = {{p: palette[i] for i, p in enumerate(pop_names)}}
colors    = [pop_color[p] for p in df["population"]]

fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.6), 4))
ax.bar(range(len(df)), df["f_roh"], color=colors, edgecolor="white", linewidth=0.5)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["sample"], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("F_ROH")
ax.set_ylim(0, max(df["f_roh"].max() * 1.15, 0.01))
ax.set_title("Runs of homozygosity \u2014 " + ref_name, fontweight="bold")

legend_handles = [mpatches.Patch(color=pop_color[p], label=p) for p in pop_names]
ax.legend(handles=legend_handles, title="Population",
          frameon=False, bbox_to_anchor=(1, 1), loc="upper left")
sns.despine(ax=ax)

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
PYEOF
        """
