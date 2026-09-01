"""Merge per-sample PAV presence columns into a contig x sample matrix, then
compute per-population presence frequencies and pairwise Hudson PAV FST with a
label-permutation null.

The permutation null is the negative control the PAV analysis needs: shuffle the
population labels across samples, recompute the genome-wide Hudson FST, repeat
`n_perm` times, and report where the observed FST falls in that distribution
(empirical p and a z-score). If observed PAV FST is not beyond the permuted
null, the "PAV separates populations" claim is not supported -- so this is what
turns a bare FST number into a defensible result.

(A hetspec-style mismatched-reference null does not apply to PAV: accessory
contigs exist only in the augref, so there is no accessory presence/absence to
score against hetspec. The permutation null is the appropriate substitute.)
"""
import numpy as np
import pandas as pd

col_files    = snakemake.input.cols
out_matrix   = snakemake.output.matrix
out_pop_freq = snakemake.output.pop_freq
out_fst      = snakemake.output.fst
pop_samples  = dict(snakemake.params.pop_samples)
pop_pairs    = [tuple(p) for p in snakemake.params.pop_pairs]
n_perm       = int(snakemake.params.n_perm)
seed         = int(snakemake.params.seed)

from pathlib import Path

# --- merge presence columns into one contig x sample 0/1 matrix ---
mat = None
for f in col_files:
    col = pd.read_csv(f, sep="\t")
    mat = col if mat is None else mat.merge(col, on="contig", how="outer")
mat = mat.fillna(0)
sample_cols = [c for c in mat.columns if c != "contig"]
mat[sample_cols] = mat[sample_cols].astype(int)
Path(out_matrix).parent.mkdir(parents=True, exist_ok=True)
mat.to_csv(out_matrix, sep="\t", index=False)

# presence matrix as a float array (contigs x samples), columns aligned to sample_cols
P = mat[sample_cols].to_numpy(dtype=float)
col_index = {s: i for i, s in enumerate(sample_cols)}

# --- per-population presence frequency per contig ---
pops = {p: [s for s in members if s in col_index] for p, members in pop_samples.items()}
freq = pd.DataFrame({"contig": mat["contig"]})
for p, members in pops.items():
    if members:
        idx = [col_index[s] for s in members]
        freq[p] = P[:, idx].mean(axis=1)
freq.to_csv(out_pop_freq, sep="\t", index=False)


def hudson_fst(p1, n1, p2, n2):
    """Genome-wide Hudson FST (ratio of averages) from per-contig presence freqs.
    p1, p2 are per-contig frequency vectors; n1, n2 the sample counts."""
    if n1 < 2 or n2 < 2:
        return float("nan")
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    den_sum = den.sum()
    return float(num.sum() / den_sum) if den_sum else float("nan")


def freqs(idx_a, idx_b):
    p1 = P[:, idx_a].mean(axis=1) if idx_a else np.zeros(P.shape[0])
    p2 = P[:, idx_b].mean(axis=1) if idx_b else np.zeros(P.shape[0])
    return p1, p2


rng = np.random.default_rng(seed)
rows = []
for (pa, pb) in pop_pairs:
    a = [col_index[s] for s in pops.get(pa, [])]
    b = [col_index[s] for s in pops.get(pb, [])]
    na, nb = len(a), len(b)
    p1, p2 = freqs(a, b)
    obs = hudson_fst(p1, na, p2, nb)

    # Permutation null: shuffle labels across the pooled a+b samples.
    pooled = a + b
    null = np.empty(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(pooled)
        pa_idx, pb_idx = list(perm[:na]), list(perm[na:])
        q1, q2 = freqs(pa_idx, pb_idx)
        null[k] = hudson_fst(q1, na, q2, nb)
    null = null[~np.isnan(null)]
    if null.size and not np.isnan(obs):
        p_emp = (1 + np.sum(null >= obs)) / (1 + null.size)
        z = (obs - null.mean()) / null.std() if null.std() > 0 else float("nan")
    else:
        p_emp, z = float("nan"), float("nan")
    rows.append({
        "pop_a": pa, "pop_b": pb, "pav_fst": obs,
        "null_mean": float(null.mean()) if null.size else float("nan"),
        "null_sd": float(null.std()) if null.size else float("nan"),
        "perm_p": p_emp, "z_score": z, "n_perm": int(null.size),
    })

pd.DataFrame(rows).to_csv(out_fst, sep="\t", index=False)
print(f"PAV FST + {n_perm}-permutation null written for {len(rows)} pairs")
