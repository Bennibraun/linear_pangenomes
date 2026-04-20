import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df       = pd.read_csv(snakemake.input.afs, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
ref_name = snakemake.params.ref

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
