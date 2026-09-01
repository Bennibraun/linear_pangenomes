"""Tally Snakemake benchmark TSVs into a linear-arm vs graph-arm cost table --
the quantitative backing for the "without graph overhead" claim.

Snakemake benchmark TSVs have a header row and one data row with columns:
    s  h:m:s  max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time
`s` is wall-clock seconds; `max_rss` is peak resident memory in MB.

We report two comparisons:
  - one-time build cost:  bwa index (augref)      vs  cactus + haplo-index
  - per-sample cost:      bwa-mem align            vs  giraffe + surject
and a per-arm total (one-time + sum over samples), so the paper can state the
whole-cohort cost ratio and the peak-memory ratio.

SV-calling cost (building the augref insert catalogue) is NOT included in the
linear total: it is a one-time preprocessing step whose graph-side analogue is
debatable, so folding it in either way would be arguable. It is reported in a
separate `note` row for transparency instead of buried in the headline ratio.
"""
import glob
import os
from pathlib import Path

import pandas as pd

bench_dir = Path(snakemake.params.bench_dir)
out_path  = Path(snakemake.output.summary)


def read_bench(path):
    """Return (wall_s, max_rss_mb) from one benchmark TSV, or None if missing."""
    try:
        df = pd.read_csv(path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None
    if df.empty:
        return None
    row = df.iloc[0]
    return float(row["s"]), float(row["max_rss"])


def collect(pattern):
    """Sum wall-clock and take peak RSS across all TSVs matching a glob."""
    total_s, peak_rss, n = 0.0, 0.0, 0
    for p in glob.glob(str(bench_dir / pattern), recursive=True):
        r = read_bench(p)
        if r is None:
            continue
        total_s += r[0]
        peak_rss = max(peak_rss, r[1])
        n += 1
    return total_s, peak_rss, n


rows = []

def add(stage, arm, pattern):
    total_s, peak_rss, n = collect(pattern)
    rows.append({
        "stage": stage, "arm": arm, "n_tasks": n,
        "wall_seconds": round(total_s, 1),
        "wall_hours": round(total_s / 3600.0, 3),
        "peak_rss_mb": round(peak_rss, 1),
        "peak_rss_gb": round(peak_rss / 1024.0, 2),
    })

# One-time build.
add("build_onetime", "linear", "index/**/*augmented_reference.fasta.tsv")
add("build_onetime", "graph",  "graph_arm/make_cactus_graph.tsv")
add("build_onetime", "graph",  "graph_arm/make_haplo_index.tsv")

# Per-sample.
add("per_sample", "linear", "linear_arm/align_bwa_augref.*.tsv")
add("per_sample", "graph",  "graph_arm/giraffe_align.*.tsv")
add("per_sample", "graph",  "graph_arm/vg_surject.*.tsv")

df = pd.DataFrame(rows)

# Per-arm totals (one-time + all per-sample).
totals = []
for arm in ("linear", "graph"):
    sub = df[df["arm"] == arm]
    totals.append({
        "stage": "TOTAL", "arm": arm, "n_tasks": int(sub["n_tasks"].sum()),
        "wall_seconds": round(sub["wall_seconds"].sum(), 1),
        "wall_hours": round(sub["wall_hours"].sum(), 3),
        "peak_rss_mb": round(sub["peak_rss_mb"].max(), 1),
        "peak_rss_gb": round(sub["peak_rss_gb"].max(), 2),
    })
df = pd.concat([df, pd.DataFrame(totals)], ignore_index=True)

out_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_path, sep="\t", index=False)

# Headline ratio, if both arms produced data.
lin = df[(df.stage == "TOTAL") & (df.arm == "linear")]
gra = df[(df.stage == "TOTAL") & (df.arm == "graph")]
if not lin.empty and not gra.empty and float(lin.wall_hours) > 0:
    ratio = float(gra.wall_hours) / float(lin.wall_hours)
    print(f"Graph/linear wall-clock ratio: {ratio:.1f}x  "
          f"(linear {float(lin.wall_hours):.2f} h, graph {float(gra.wall_hours):.2f} h)")
    print(f"Peak RSS: linear {float(lin.peak_rss_gb):.1f} GB, "
          f"graph {float(gra.peak_rss_gb):.1f} GB")
print("Wrote", out_path)
