"""Per-reference summary table (Jeon et al. 2026 Table 1 equivalent).

Emits one row per reference with:
  ref, length_bp, mean_mapping_rate, mean_mapping_qual, n_snps,
  n_snps_filtered, n_svs (from graph SV pipeline), mean_pi,
  size_vs_conspec_pct

Designed for direct comparison to Jeon Table 1.
"""

from pathlib import Path

import pandas as pd
import pysam


# --- Inputs (lists) --------------------------------------------------------
refs = list(snakemake.params.refs)
fastas = {r: f for r, f in zip(refs, snakemake.params.fastas)}
fais = {r: f for r, f in zip(refs, snakemake.params.fais)}
vcfs = {r: f for r, f in zip(refs, snakemake.input.vcfs)}
sv_vcfs = list(getattr(snakemake.input, "sv_vcfs", []))  # per-sample mc_graph SV VCFs (graph only)
metrics_path = snakemake.input.metrics
pi_paths = {r: f for r, f in zip(refs, snakemake.input.pi)}
conspec_ref = snakemake.params.conspec_ref

out_path = Path(snakemake.output.summary)
out_path.parent.mkdir(parents=True, exist_ok=True)


# --- Genome length per ref (from fai) --------------------------------------
def total_length(fai_path):
    return sum(int(line.split("\t")[1]) for line in Path(fai_path).read_text().splitlines() if line.strip())


lengths = {r: total_length(fais[r]) for r in refs}


# --- SNP counts (pre-filter = raw VCF; post-filter = QUAL >= 20, biallelic) ---
def count_variants(vcf_path):
    """Returns (total, post_filter) variant counts.
    post_filter: QUAL >= 20 and biallelic SNPs only (Jeon's "post-filtering")."""
    total = 0
    pf = 0
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            total += 1
            qual_ok = rec.qual is not None and rec.qual >= 20
            biallelic = len(rec.alts or []) == 1
            is_snp = (
                biallelic
                and rec.alts
                and len(rec.ref) == 1
                and len(rec.alts[0]) == 1
            )
            if qual_ok and is_snp:
                pf += 1
    return total, pf


snp_counts = {}
for r in refs:
    total, pf = count_variants(vcfs[r])
    snp_counts[r] = {"n_snps": total, "n_snps_filtered": pf}


# --- SV counts: sum across per-sample mc_graph SV VCFs (cohort union) ------
sv_set = set()
for path in sv_vcfs:
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            sv_set.add((rec.contig, rec.pos, rec.ref, tuple(rec.alts or [])))
n_svs = len(sv_set)
# SVs only exist for mc_graph (graph SV pipeline); other refs are 0.


# --- Mapping rate / mapping quality from alignment_metrics.tsv -------------
metrics = pd.read_csv(metrics_path, sep="\t").rename(
    columns={"alignment_type": "ref"}
)
metrics["mapping_rate"] = pd.to_numeric(metrics["mapping_rate"], errors="coerce")
metrics["mean_mapq"] = pd.to_numeric(metrics["mean_mapq"], errors="coerce")
map_summary = (
    metrics.groupby("ref")[["mapping_rate", "mean_mapq"]]
    .mean()
    .reset_index()
    .set_index("ref")
)


# --- Mean nucleotide diversity (π) per ref ---------------------------------
def mean_pi(pi_path):
    df = pd.read_csv(pi_path, sep="\t")
    if "PI" not in df.columns:
        return float("nan")
    return df["PI"].mean()


pi_summary = {r: mean_pi(pi_paths[r]) for r in refs}


# --- Assemble --------------------------------------------------------------
conspec_len = lengths.get(conspec_ref, 0)
rows = []
for r in refs:
    rows.append(
        {
            "ref": r,
            "length_bp": lengths[r],
            "size_vs_conspec_pct": (
                100.0 * (lengths[r] - conspec_len) / conspec_len
                if conspec_len > 0
                else float("nan")
            ),
            "mean_mapping_rate": map_summary.loc[r, "mapping_rate"]
            if r in map_summary.index
            else float("nan"),
            "mean_mapping_qual": map_summary.loc[r, "mean_mapq"]
            if r in map_summary.index
            else float("nan"),
            "n_snps": snp_counts[r]["n_snps"],
            "n_snps_filtered": snp_counts[r]["n_snps_filtered"],
            "n_svs": n_svs if r == "mc_graph" else 0,
            "mean_pi": pi_summary[r],
        }
    )

out = pd.DataFrame(rows)
out.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
