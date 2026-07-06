"""Per-sample QC table.

Aggregates per-sample metrics across all references into a single TSV:
  sample, ref, mean_depth, mapping_rate, n_het, n_missing, het_rate,
  qc_status, qc_flags

Each metric is checked against a configurable threshold from
config['sample_qc']. qc_status is "PASS" if all checks pass for that
(sample, ref), else "FAIL". qc_flags lists the specific failed checks.

This is report-only: failed samples are still consumed by downstream
analyses. The user decides what to exclude based on this table.

mc_graph rows pull mean_depth from the heterozygosity table's averageDepth
column (since BAM-flagstat mean_depth is GAM-derived/NA for that ref).
"""

from pathlib import Path

import pandas as pd

metrics_path = snakemake.input.metrics
het_paths = list(snakemake.input.heterozygosity)
out_path = Path(snakemake.output.qc)

min_depth = float(snakemake.params.min_depth)
min_mapping_rate = float(snakemake.params.min_mapping_rate)
min_het_rate = float(snakemake.params.min_het_rate)
max_het_rate = float(snakemake.params.max_het_rate)
max_missing_rate = float(snakemake.params.max_missing_rate)


# --- Load alignment metrics (mean_depth, mapping_rate) ---------------------
metrics = pd.read_csv(metrics_path, sep="\t")
metrics = metrics.rename(columns={"alignment_type": "ref"})
metrics["mean_depth"] = pd.to_numeric(metrics["mean_depth"], errors="coerce")
metrics["mapping_rate"] = pd.to_numeric(metrics["mapping_rate"], errors="coerce")
# Drop the duplicate "cactus" rows (mc_graph has its own row that we keep).
metrics = metrics[metrics["ref"] != "cactus"][
    ["sample", "ref", "mean_depth", "mapping_rate"]
]


# --- Load heterozygosity (n_het, n_missing, het_rate, mean_depth backup) ---
het_frames = []
for path in het_paths:
    h = pd.read_csv(path, sep="\t")
    het_frames.append(h)
het = pd.concat(het_frames, ignore_index=True) if het_frames else pd.DataFrame(
    columns=["sample", "ref", "n_het", "n_missing", "mean_depth", "het_rate"]
)
het = het.rename(columns={"mean_depth": "het_mean_depth"})
het = het[["sample", "ref", "n_ref_hom", "n_nonref_hom", "n_het", "n_missing", "het_mean_depth", "het_rate"]]


# --- Merge -----------------------------------------------------------------
df = metrics.merge(het, on=["sample", "ref"], how="outer")

# For mc_graph mean_depth is NA in alignment_metrics (GAM-derived); fall back
# to the per-sample average depth bcftools stats computed from the surjected
# BAM-backed VCF.
df["mean_depth"] = df["mean_depth"].fillna(df["het_mean_depth"])

# Missingness rate: n_missing / (n_called_or_missing). bcftools stats's PSC
# block only gives counts, not a denominator; we approximate the denom as
# (called + missing) per sample.
total_sites = (df["n_ref_hom"].fillna(0) + df["n_nonref_hom"].fillna(0)
               + df["n_het"].fillna(0) + df["n_missing"].fillna(0))
df["missing_rate"] = (df["n_missing"].fillna(0) / total_sites).where(
    total_sites > 0, 0.0
)

df = df.drop(columns=["het_mean_depth", "n_ref_hom", "n_nonref_hom"])


# --- Apply thresholds ------------------------------------------------------
def evaluate(row):
    flags = []
    if pd.notna(row["mean_depth"]) and row["mean_depth"] < min_depth:
        flags.append(f"low_depth({row['mean_depth']:.1f}<{min_depth:g})")
    if pd.notna(row["mapping_rate"]) and row["mapping_rate"] < min_mapping_rate:
        flags.append(f"low_mapping_rate({row['mapping_rate']:.3f}<{min_mapping_rate:g})")
    if pd.notna(row["het_rate"]):
        if row["het_rate"] < min_het_rate:
            flags.append(f"low_het({row['het_rate']:.4f}<{min_het_rate:g})")
        elif row["het_rate"] > max_het_rate:
            flags.append(f"high_het({row['het_rate']:.4f}>{max_het_rate:g})")
    if pd.notna(row["missing_rate"]) and row["missing_rate"] > max_missing_rate:
        flags.append(f"high_missing({row['missing_rate']:.3f}>{max_missing_rate:g})")
    return ("PASS" if not flags else "FAIL", ";".join(flags) if flags else "")


df[["qc_status", "qc_flags"]] = df.apply(
    lambda r: pd.Series(evaluate(r)), axis=1
)

# --- Order columns and write ------------------------------------------------
cols = [
    "sample", "ref", "mean_depth", "mapping_rate",
    "n_het", "n_missing", "het_rate", "missing_rate",
    "qc_status", "qc_flags",
]
df = df[cols].sort_values(["sample", "ref"])

out_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
