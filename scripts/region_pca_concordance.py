"""E5: Procrustes + Mantel concordance between the core-only and accessory-only
sample configurations (augref region-split PCA).

High concordance => accessory variation co-segregates with core demography.
Also reports per-morph accessory depth so the coverage confound stays on record
(surface accessory depth is ~half cave's, so per-morph configuration differences
are coverage-driven, not necessarily biological).

Consumes the region PCAngsd covariances + region_qc_summary; writes one tidy
stats TSV. No heavy inputs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr

cov_core = np.loadtxt(snakemake.input.cov_core)
cov_acc = np.loadtxt(snakemake.input.cov_accessory)
samples_core = [s.strip() for s in Path(snakemake.input.samples_core).read_text().splitlines() if s.strip()]
samples_acc = [s.strip() for s in Path(snakemake.input.samples_accessory).read_text().splitlines() if s.strip()]
out_path = Path(snakemake.output.stats)

if samples_core != samples_acc:
    raise ValueError("core/accessory sample order differs; cannot pair configurations")


def pca_coords(cov, k):
    cov = (cov + cov.T) / 2
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    coords = vecs * np.sqrt(np.clip(vals, 0, None))
    return coords[:, :k]


def mantel(A, B, n=9999, seed=0):
    """Mantel test: correlation of the two pairwise-distance matrices, with a
    permutation p-value. Returns (observed r, p)."""
    a, b = pdist(A), pdist(B)
    r_obs = pearsonr(a, b)[0]
    Bsq = squareform(b)
    rng = np.random.default_rng(seed)
    ns = Bsq.shape[0]
    cnt = 0
    for _ in range(n):
        perm = rng.permutation(ns)
        b_perm = squareform(Bsq[np.ix_(perm, perm)])
        if abs(pearsonr(a, b_perm)[0]) >= abs(r_obs):
            cnt += 1
    return r_obs, (cnt + 1) / (n + 1)


rows = []
for k in (2, 4, 6):
    Ccore = pca_coords(cov_core, k)
    Cacc = pca_coords(cov_acc, k)
    _, _, disparity = procrustes(Ccore, Cacc)
    r, p = mantel(Ccore, Cacc)
    rows.append({
        "n_pcs": k,
        "procrustes_disparity_m2": round(float(disparity), 6),
        "mantel_r": round(float(r), 6),
        "mantel_p": round(float(p), 6),
    })

stats = pd.DataFrame(rows)

# --- coverage confound record: per-morph accessory depth --------------------
SURFACE = {"Riochoy", "Mante", "Rascon", "Choy", "Rio"}


def pop_of(sample_id):
    return sample_id.split("_")[0]


def morph_of(sample_id):
    p = pop_of(sample_id).lower()
    surface_prefixes = {"choy", "rio", "riochoy", "mante", "rascon"}
    return "surface" if p in surface_prefixes else "cave"


rq = pd.read_csv(snakemake.input.region_qc, sep="\t")
acc = rq[rq["region"] == "accessory"].copy()
acc["morph"] = acc["sample"].map(morph_of)
depth_by_morph = acc.groupby("morph")["mean_depth"].mean()

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as fh:
    stats.to_csv(fh, sep="\t", index=False)
    fh.write("\n# coverage confound (accessory mean depth by morph):\n")
    for morph, d in depth_by_morph.items():
        fh.write(f"# {morph}\t{d:.3f}x\n")
    fh.write("# NOTE: per-morph configuration differences are coverage-confounded.\n")

print(stats.to_string(index=False))
print("accessory depth by morph:\n", depth_by_morph.round(3).to_string())
