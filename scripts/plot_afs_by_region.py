"""Core vs accessory folded SFS overlay for augref. If the accessory sequence
carries variation shared broadly across the cohort (user's hypothesis), its
spectrum should be shifted toward common variants relative to core. If instead
low depth/mappability in accessory regions is undercalling rare variants
(the confound), the same rightward shift appears but should track the
depth/MAPQ gap in region_qc_summary.tsv rather than reflecting real
demography — this plot alone can't distinguish the two, see plot_region_qc.py."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

df      = pd.read_csv(snakemake.input.afs_by_region, sep="\t")
out_pdf = snakemake.output.pdf
out_png = snakemake.output.png

df["bin_low"]  = df["bin_low"].astype(float)
df["bin_high"] = df["bin_high"].astype(float)
df["count"]    = df["count"].astype(int)

for region in ("core", "accessory"):
    total = df.loc[df["region"] == region, "count"].sum()
    df.loc[df["region"] == region, "proportion"] = (
        df.loc[df["region"] == region, "count"] / total if total else 0.0
    )

labels = [format(r["bin_low"], ".2f") + "–" + format(r["bin_high"], ".2f")
          for _, r in df[df["region"] == "core"].iterrows()]
x = np.arange(len(labels))
width = 0.38

fig, ax = plt.subplots(figsize=(8, 4.5))
core = df[df["region"] == "core"].sort_values("bin_low")
acc  = df[df["region"] == "accessory"].sort_values("bin_low")
ax.bar(x - width / 2, core["proportion"], width, label="core", color="#4C72B0")
ax.bar(x + width / 2, acc["proportion"], width, label="accessory", color="#C44E52")
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.set_xlabel("Minor allele frequency")
ax.set_ylabel("Proportion of sites")
ax.set_title("Folded SFS: core vs accessory — augref", fontweight="bold")
ax.legend(frameon=False)
n_core = int(core["count"].sum())
n_acc  = int(acc["count"].sum())
ax.text(0.98, 0.95, f"n_core={n_core}\nn_accessory={n_acc}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
