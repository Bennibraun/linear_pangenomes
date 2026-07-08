"""Per-reference summary table.

Emits one row per reference with:
  ref, length_bp, size_vs_conspec_pct, graph_total_bp (graph only),
  mean_mapping_rate, mean_mapping_qual, n_snps, n_snps_filtered, n_svs,
  mean_pi, mean_tajimas_d, mean_fst

SV counts come from the ref-appropriate source (these are INDEPENDENT discovery
paths, not the same set):
  - augref  : long-read pan-sample SV catalog (sniffles/cuteSV -> SURVIVOR)
  - mc_graph: its own vg-call VCFs
  - conspec/hetspec: 0 (plain references, no SVs added)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pysam


# --- Inputs (lists) --------------------------------------------------------
refs = list(snakemake.params.refs)
fastas = {r: f for r, f in zip(refs, snakemake.params.fastas)}
fais = {r: f for r, f in zip(refs, snakemake.params.fais)}
vcfs = {r: f for r, f in zip(refs, snakemake.input.vcfs)}
sv_vcfs = list(getattr(snakemake.input, "sv_vcfs", []))  # per-sample mc_graph SV VCFs (graph only)
sv_catalog = snakemake.input.sv_catalog                  # augref long-read SV catalog
metrics_path = snakemake.input.metrics
pi_paths = {r: f for r, f in zip(refs, snakemake.input.pi)}
tajima_paths = {r: f for r, f in zip(refs, snakemake.input.tajima)}
conspec_ref = snakemake.params.conspec_ref
pop_pairs = [tuple(p) for p in snakemake.params.pop_pairs]
include_graph = bool(snakemake.params.include_graph)

# FST inputs are flattened (ref x pair). Map each to its ref by path — the FST
# files live under fst/<ref>/<pair>.weir.fst — so we don't depend on Snakemake
# preserving input order.
fst_input = list(getattr(snakemake.input, "fst", []))
fst_paths = {r: [] for r in refs}
for path in fst_input:
    parts = Path(path).parts
    # .../fst/<ref>/<pair>.weir.fst  -> ref is the parent directory name.
    ref = Path(path).parent.name
    if ref in fst_paths:
        fst_paths[ref].append(path)

# Graph total length (precomputed by graph_total_length rule; graph only).
graph_length_files = list(getattr(snakemake.input, "graph_length", []))
graph_total_bp = np.nan
if graph_length_files:
    try:
        graph_total_bp = int(Path(graph_length_files[0]).read_text().split()[0])
    except (ValueError, IndexError):
        graph_total_bp = np.nan

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


# --- SV counts (ref-appropriate, independent discovery paths) --------------
# Both are size-filtered to >= min_sv_size so counts are comparable. The long-
# read catalog (augref) is already filtered upstream (sniffles/cuteSV
# --minsvlen), but vg call for mc_graph filters by SNARL traversal length, not
# by the emitted record's allele-length change, so its raw VCF still contains
# SNP- and small-indel-scale records. We filter by allele-length difference here.
min_sv_size = int(snakemake.params.min_sv_size)


def record_sv_len(rec):
    """Size of a VCF record = max allele-length difference vs REF.
    Prefer INFO/SVLEN when the header declares it, else compute from REF/ALT
    lengths. pysam raises on rec.info access for undeclared INFO keys (the vg
    call VCFs don't declare SVLEN), so we check the header first."""
    if "SVLEN" in rec.header.info:
        svlen = rec.info.get("SVLEN")
        if svlen is not None:
            vals = svlen if isinstance(svlen, (tuple, list)) else [svlen]
            nums = [abs(int(v)) for v in vals if v is not None]
            if nums:
                return max(nums)
    ref_len = len(rec.ref) if rec.ref else 0
    return max((abs(len(a) - ref_len) for a in (rec.alts or [])), default=0)


# mc_graph: union across per-sample vg-call VCFs, SVs >= min_sv_size only.
mc_graph_sv_set = set()
for path in sv_vcfs:
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            if record_sv_len(rec) >= min_sv_size:
                mc_graph_sv_set.add((rec.contig, rec.pos, rec.ref, tuple(rec.alts or [])))
n_svs_mc_graph = len(mc_graph_sv_set)


# augref: count records in the long-read pan-sample SV catalog (SURVIVOR VCF).
# Already >= min_sv_size upstream, but we apply the same filter for consistency.
def count_svs(path):
    n = 0
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            if record_sv_len(rec) >= min_sv_size:
                n += 1
    return n


n_svs_augref = count_svs(sv_catalog)


def n_svs_for(ref):
    if ref == "augref":
        return n_svs_augref
    if ref == "mc_graph":
        return n_svs_mc_graph
    return 0  # conspec / hetspec: plain references


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


# --- Mean Tajima's D per ref -----------------------------------------------
def mean_tajd(path):
    df = pd.read_csv(path, sep="\t")
    if "TajimaD" not in df.columns:
        return float("nan")
    return pd.to_numeric(df["TajimaD"], errors="coerce").mean()


tajd_summary = {r: mean_tajd(tajima_paths[r]) for r in refs}


# --- Mean FST per ref (averaged over pairs and windows) --------------------
_FST_NA = ["-nan", "nan", "-NaN", "NaN", "-inf", "inf"]


def mean_fst(paths):
    vals = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", na_values=_FST_NA)
        col = "WEIGHTED_FST" if "WEIGHTED_FST" in df.columns \
            else "WEIR_AND_COCKERHAM_FST"
        if col not in df.columns:
            continue
        v = pd.to_numeric(df[col], errors="coerce").dropna().clip(lower=0)
        vals.extend(v.tolist())
    return float(np.mean(vals)) if vals else float("nan")


fst_summary = {r: mean_fst(fst_paths[r]) for r in refs}


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
            # Total pangenome graph sequence (mc_graph only); distinct from
            # length_bp, which is the conspec coordinate space used for calling.
            "graph_total_bp": graph_total_bp if r == "mc_graph" else np.nan,
            "mean_mapping_rate": map_summary.loc[r, "mapping_rate"]
            if r in map_summary.index
            else float("nan"),
            "mean_mapping_qual": map_summary.loc[r, "mean_mapq"]
            if r in map_summary.index
            else float("nan"),
            "n_snps": snp_counts[r]["n_snps"],
            "n_snps_filtered": snp_counts[r]["n_snps_filtered"],
            "n_svs": n_svs_for(r),
            "mean_pi": pi_summary[r],
            "mean_tajimas_d": tajd_summary[r],
            "mean_fst": fst_summary[r],
        }
    )

out = pd.DataFrame(rows)
out.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
