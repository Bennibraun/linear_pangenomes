import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

df       = pd.read_csv(snakemake.input.summary, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref
sample_pop_str = snakemake.params.sample_pop

sample_to_pop = {}
for pair in sample_pop_str.split("|"):
    k, v = pair.split(",", 1)
    sample_to_pop[k] = v

df = df[df["ref"] == ref_name].copy()
df["population"] = df["sample"].map(sample_to_pop).fillna("unknown")
df = df.sort_values(["population", "sample"])

pop_names = sorted(df["population"].unique())
palette   = sns.color_palette("tab10", len(pop_names))
pop_color = {p: palette[i] for i, p in enumerate(pop_names)}
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
