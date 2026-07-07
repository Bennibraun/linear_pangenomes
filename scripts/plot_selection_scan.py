"""PCAngsd selection-scan summary: QQ plot + ranked outlier table.

The per-SNP chi-squared selection statistics (from pcangsd --selection, one test
per PC, aggregated) are better summarised by a QQ plot than a manhattan: the QQ
shows genome-wide inflation and exactly where the outlier tail departs from the
chi-squared null, in one compact panel. Alongside it we list the top outlier
SNPs by significance so the actual candidates are legible.

Input columns (fst_outliers.tsv): chrom, pos, id, chi2, best_pc, pval, qval,
outlier, plus per-PC chi2/pval. SV_* accessory contigs are excluded.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

df       = pd.read_csv(snakemake.input.outliers, sep="\t")
out_pdf  = snakemake.output.pdf
out_png  = snakemake.output.png
out_tsv  = snakemake.output.top_outliers
ref_name = snakemake.params.ref
fdr      = snakemake.params.fdr


def _empty(msg):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    pd.DataFrame(columns=["chrom", "pos", "chi2", "pval", "qval"]).to_csv(
        out_tsv, sep="\t", index=False
    )
    raise SystemExit(0)


df = df[~df["chrom"].astype(str).str.startswith("SV_")].copy()
df["chi2"] = pd.to_numeric(df["chi2"], errors="coerce")
df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
df = df.dropna(subset=["chi2", "pval"])
if df.empty:
    _empty(f"No selection-scan sites — {ref_name}")

# --- Ranked outlier table (top candidates by combined p, then chi2) ---------
top = df.sort_values(["pval", "chi2"], ascending=[True, False]).head(50)
top_cols = [c for c in ["chrom", "pos", "id", "chi2", "best_pc", "pval", "qval",
                        "outlier"] if c in top.columns]
top[top_cols].to_csv(out_tsv, sep="\t", index=False)

# --- QQ panel: observed vs expected -log10(p) under the null ----------------
p = np.clip(df["pval"].to_numpy(), 1e-300, 1.0)
p_sorted = np.sort(p)
n = p_sorted.size
expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
observed = -np.log10(p_sorted)

# Genomic inflation factor lambda from the median chi2 (df=1).
lam = np.median(df["chi2"].to_numpy()) / stats.chi2.ppf(0.5, df=1)

fig, (ax_qq, ax_hist) = plt.subplots(1, 2, figsize=(13, 5))

lim = max(expected.max(), observed.max()) * 1.05
ax_qq.plot([0, lim], [0, lim], color="grey", lw=1.0, linestyle="--")
n_out = int(df["outlier"].sum()) if "outlier" in df.columns else 0
is_out = df.sort_values("pval")["outlier"].to_numpy() if "outlier" in df.columns \
    else np.zeros(n, dtype=bool)
ax_qq.scatter(expected[~is_out], observed[~is_out], s=8, color="#4C72B0",
              alpha=0.6, linewidths=0, rasterized=True, label="sites")
if is_out.any():
    ax_qq.scatter(expected[is_out], observed[is_out], s=18, color="crimson",
                  zorder=5, label=f"FDR < {fdr} (n={n_out})")
ax_qq.set_xlabel(r"expected $-\log_{10}(p)$")
ax_qq.set_ylabel(r"observed $-\log_{10}(p)$")
ax_qq.set_title(f"QQ (λ = {lam:.2f})")
ax_qq.legend(frameon=False, fontsize=8)

# --- chi2 distribution vs null ----------------------------------------------
ax_hist.hist(df["chi2"], bins=80, density=True, color="#4C72B0", alpha=0.6,
             label="observed χ²")
xs = np.linspace(0, np.percentile(df["chi2"], 99.5), 200)
ax_hist.plot(xs, stats.chi2.pdf(xs, df=1), color="black", lw=1.2,
             label="χ²(df=1) null")
ax_hist.set_xlabel("selection χ² (max across PCs)")
ax_hist.set_ylabel("density")
ax_hist.legend(frameon=False, fontsize=8)

fig.suptitle(f"Selection scan (PCAngsd) — {ref_name}", fontweight="bold")
sns.despine(fig=fig)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=150)
plt.close(fig)
