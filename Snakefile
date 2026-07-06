import os

# Cap implicit BLAS / numerical-library thread fans at 1 BEFORE pandas/numpy
# import below — set in the snakemake parent process so every subprocess
# (shell, script, run) inherits it. Without this, anything that touches
# numpy/scipy/sklearn/matplotlib (vg, pcangsd, plotting scripts) tries to
# spawn one OpenBLAS thread per host core and trips RLIMIT_NPROC inside the
# slurm cgroup. shell.prefix below is belt-and-suspenders for shell rules.
for _v in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_v, "1")

from pathlib import Path

_scratch_tmpdir = str(Path("results/tmp").resolve())
Path(_scratch_tmpdir).mkdir(parents=True, exist_ok=True)
for _v in ("TMPDIR", "TEMPDIR", "TMP", "TEMP"):
    os.environ[_v] = _scratch_tmpdir
import pandas as pd

configfile: "config/config.yaml"

# ============================================================================
# Parse reads manifest (TSV with columns: sample_id, seq_type, platform, grouping, fastq_r1, fastq_r2)
# ============================================================================
INPUTS_CFG = config["inputs"]
MANIFEST_FILE = INPUTS_CFG["reads_manifest"]

# Read manifest
manifest_df = pd.read_csv(MANIFEST_FILE, sep="\t")

# Extract sample lists by seq_type
LONG_SAMPLES = manifest_df[manifest_df["seq_type"] == "long"]["sample_id"].tolist()
SHORT_SAMPLES = manifest_df[manifest_df["seq_type"] == "short"]["sample_id"].tolist()

# Build lookup tables for fastq paths
LONG_READS = {row["sample_id"]: row["fastq_r1"] for _, row in manifest_df[manifest_df["seq_type"] == "long"].iterrows()}
PLATFORM_MAP = {row["sample_id"]: row["platform"].upper() for _, row in manifest_df[manifest_df["seq_type"] == "long"].iterrows()}
SHORT_READS_R1 = {row["sample_id"]: row["fastq_r1"] for _, row in manifest_df[manifest_df["seq_type"] == "short"].iterrows()}
SHORT_READS_R2 = {row["sample_id"]: row["fastq_r2"] for _, row in manifest_df[manifest_df["seq_type"] == "short"].iterrows()}

# Population/grouping structure
SAMPLE_TO_POP = dict(zip(manifest_df["sample_id"], manifest_df["grouping"]))
POP_NAMES = sorted(manifest_df["grouping"].unique().tolist())
POP_SAMPLES = {pop: manifest_df[manifest_df["grouping"] == pop]["sample_id"].tolist() for pop in POP_NAMES}

# ============================================================================
# Reference genome configurations
# ============================================================================
REFS_NESTED = config["references"]
REFERENCE_NAMES = list(REFS_NESTED.keys())
ALIGN_AUGREF = REFS_NESTED["augref"]["fasta"]
ALIGN_CONSPEC = REFS_NESTED["conspec"]["fasta"]
ALIGN_HETSPEC = REFS_NESTED["hetspec"]["fasta"]

# ============================================================================
# Optional-stage toggles (default OFF)
# ============================================================================
# INCLUDE_GRAPH gates the entire Minigraph-Cactus graph path: the cactus build,
# haplotype-sampling index, giraffe alignment, surjection, vg SV calls, and
# mc_graph as a reference in the popgen/plot rules. When off, mc_graph is not
# added to REFERENCE_NAMES, so every {ref}-wildcarded rule stops expanding it
# automatically. SV calling / augref are independent and stay on.
# INCLUDE_ASSEMBLY_ALIGN gates the short-reads-to-assembly QC metrics + plot
# (NOT the assembly-to-reference PAF that SV calling needs, and NOT the
# assemblies themselves, which SV/augref/graph depend on).
INCLUDE_GRAPH = config.get("include_graph", False)
INCLUDE_ASSEMBLY_ALIGN = config.get("include_assembly_align", False)

# mc_graph is the Minigraph-Cactus pangenome reference. SNPs come from
# bcftools on surjected BAMs (conspec coordinates), written to the same
# bcftools/{ref}/combined/merged.vcf.gz path as the linear refs so
# {ref}-wildcarded popgen rules (FST, AFS, π, allelic balance, PCAngsd) pick it
# up automatically. SVs are called separately via vg call on the cohort-
# augmented graph and stored under sv/vg/ for benchmarking. Only added when the
# graph is enabled.
if INCLUDE_GRAPH:
    REFS_NESTED["mc_graph"] = dict(REFS_NESTED["conspec"])
    REFERENCE_NAMES = REFERENCE_NAMES + ["mc_graph"]

def get_ref_fasta(name):
    return REFS_NESTED[name]["fasta"]

# ============================================================================
# Population pairs for FST/AFS analysis
# ============================================================================
POP_PAIRS = config["population_pairs"]
POP_PAIR_TUPLES = [(p[0], p[1]) for p in POP_PAIRS]

# ============================================================================
# Stage-specific configurations
# ============================================================================
ASSEMBLY_CFG = config["assembly"]
ASSEMBLY_OUTDIR = Path(ASSEMBLY_CFG.get("outdir", "results/assemblies"))
ASSEMBLY_LINEAGE = ASSEMBLY_CFG.get("lineage", "hymenoptera_odb10")
ASSEMBLY_THREADS = ASSEMBLY_CFG.get("threads", 1)

SV_CFG = config["sv_calling"]
SV_OUTDIR = Path(SV_CFG.get("outdir", "results/sv_calls"))
SV_ASSEMBLY_DIR = SV_CFG.get("assembly_dir", str(ASSEMBLY_OUTDIR / "hifiasm"))
SV_THREADS = SV_CFG.get("threads", 8)
SV_SURVIVOR = SV_CFG.get("survivor_exec", "SURVIVOR")
SV_MIN_SIZE = SV_CFG.get("min_sv_size", 50)
SV_MIN_SUPPORT = SV_CFG.get("min_read_support", 3)
SV_BREAKPOINT_SLOP = SV_CFG.get("breakpoint_slop", 1000)
SV_JASMINE_SLOP = SV_CFG.get("jasmine_slop", 500)
SV_FLANK = SV_CFG.get("flank", 200)

CACTUS_CFG = config["cactus"]
CACTUS_IMAGE = CACTUS_CFG["image"]
CACTUS_SEQFILE = Path("results") / "cactus_work" / "seqfile.txt"  # Auto-generated, not manual
CACTUS_OUTDIR = Path(CACTUS_CFG.get("outdir", "results/cactus"))
CACTUS_JOBSTORE = CACTUS_OUTDIR.parent / "cactus_jobstore"
CACTUS_OUTNAME = CACTUS_CFG.get("outname", "cactus_graph")
CACTUS_MAX_CORES = CACTUS_CFG.get("max_cores", 8)
CACTUS_REF_CONTIGS = CACTUS_CFG.get("ref_contigs", "")
CACTUS_EXTRA_ARGS = CACTUS_CFG.get("extra_args", "")
# Reference sequence names to drop from the graph reference only (not the linear
# pipeline). Used to exclude circular/organellar contigs (e.g. the mitochondrion)
# that break `vg haplotypes` haplotype sampling with "top-level chain is a loop".
CACTUS_EXCLUDE_REF_SEQS = CACTUS_CFG.get("exclude_ref_seqs", []) or []
if isinstance(CACTUS_EXCLUDE_REF_SEQS, str):
    CACTUS_EXCLUDE_REF_SEQS = CACTUS_EXCLUDE_REF_SEQS.split()

VG_IMAGE = "docker://quay.io/vgteam/vg:v1.74.0"

ALIGN_CFG = config["align_wgs"]
ALIGN_OUTDIR = Path(ALIGN_CFG.get("outdir", "results/wgs_alignments"))
ALIGN_CACTUS_GBZ = ALIGN_CFG.get("cactus_gbz", str(CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz"))
ALIGN_THREADS = ALIGN_CFG.get("threads", 4)
ALIGN_DOWNSAMPLE_BASES = ALIGN_CFG.get("downsample_bases", "10G")

METRICS_CFG = config["align_metrics"]
METRICS_OUTDIR = Path(METRICS_CFG.get("outdir", "results/align_metrics"))
METRICS_REF_TYPES = REFERENCE_NAMES
METRICS_GAM = METRICS_CFG.get("include_gam", True)
METRICS_MIN_MAPQ = METRICS_CFG.get("min_mapq", 0)
METRICS_MIN_BASEQ = METRICS_CFG.get("min_baseq", 0)
METRICS_THREADS = METRICS_CFG.get("threads", 2)

VC_CFG = config["variant_calling"]
BCFTOOLS_CFG = VC_CFG["bcftools"]
VG_CFG = VC_CFG["vg"]
VC_OUTDIR = Path(VC_CFG.get("outdir", "results/variants"))
BCFTOOLS_THREADS = BCFTOOLS_CFG.get("threads", 4)
VG_GBZ = VG_CFG.get("graph_gbz", str(CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz"))
VG_THREADS = VG_CFG.get("threads", 4)


FST_CFG = config["fst_afs"]
FST_OUTDIR = Path(FST_CFG.get("outdir", "results/fst_afs"))
FST_WINDOW_SIZE = FST_CFG.get("window_size", 0)
FST_WINDOW_STEP = FST_CFG.get("window_step", 0)
FST_WINDOW_ARGS = ""
if FST_WINDOW_SIZE and FST_WINDOW_STEP:
    FST_WINDOW_ARGS = f"--fst-window-size {FST_WINDOW_SIZE} --fst-window-step {FST_WINDOW_STEP}"
elif FST_WINDOW_SIZE:
    FST_WINDOW_ARGS = f"--fst-window-size {FST_WINDOW_SIZE}"
AFS_BINS = FST_CFG.get("afs_bins", [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])

PI_CFG = config["pi"]
PI_OUTDIR = Path(PI_CFG.get("outdir", "results/pi"))
PI_WINDOW_SIZE = PI_CFG.get("window_size", 10000)
PI_WINDOW_STEP = PI_CFG.get("window_step", 5000)

AB_CFG = config["allelic_balance"]
AB_OUTDIR = Path(AB_CFG.get("outdir", "results/allelic_balance"))
AB_BINS = AB_CFG.get("bins", [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
AB_MIN_DEPTH = AB_CFG.get("min_depth", 10)

# ROH (runs of homozygosity) removed: with SV_* contigs excluded from ROH
# calling, augref ROH is identical to conspec by construction (no added
# information), so the whole ROH stage was cut. Rules, targets, config, and the
# f_roh columns in the summary tables are all disabled. To revive, uncomment
# this block, the ROH rules in rules/popgen.smk + rules/plotting.smk, the target
# lines below, the roh: config block, and restore the froh inputs/columns in
# scripts/per_ref_summary.py and scripts/sample_qc.py.
# ROH_CFG = config["roh"]
# ROH_OUTDIR = Path(ROH_CFG.get("outdir", "results/roh"))
# ROH_GENOME_LENGTHS = ROH_CFG.get("genome_lengths", {})
# ROH_MIN_LENGTH = ROH_CFG.get("min_roh_length", 100000)
# ROH_RECOMB_RATE = ROH_CFG.get("recomb_rate", 30)
# ROH_BCFTOOLS_ARGS = ROH_CFG.get("bcftools_args", "")

RELATEDNESS_CFG = config.get("relatedness", {})
RELATEDNESS_OUTDIR = Path(RELATEDNESS_CFG.get("outdir", "results/relatedness"))

HET_CFG = config.get("heterozygosity", {})
HET_OUTDIR = Path(HET_CFG.get("outdir", "results/heterozygosity"))

QC_CFG = config.get("sample_qc", {})
QC_OUTDIR = Path(QC_CFG.get("outdir", "results/sample_qc"))
QC_MIN_DEPTH = QC_CFG.get("min_depth", 5.0)
QC_MIN_MAPPING_RATE = QC_CFG.get("min_mapping_rate", 0.8)
QC_MIN_HET_RATE = QC_CFG.get("min_het_rate", 0.0001)
QC_MAX_HET_RATE = QC_CFG.get("max_het_rate", 0.05)
QC_MAX_MISSING_RATE = QC_CFG.get("max_missing_rate", 0.2)

# Reference pairs for F1 concordance. Only refs that share the conspec
# coordinate system can be intersected directly (hetspec is on its own
# coordinates so it's excluded). Override via config['f1_concordance'][
# 'compatible_refs'] for non-default reference setups.
_F1_CFG = config.get("f1_concordance", {})
_F1_COMPATIBLE = _F1_CFG.get(
    "compatible_refs",
    [r for r in REFERENCE_NAMES if r in ("conspec", "augref", "mc_graph")],
)
# Drop any refs that aren't active (e.g. mc_graph when the graph is disabled)
# so an explicit compatible_refs list can't request unbuildable outputs.
_F1_COMPATIBLE = [r for r in _F1_COMPATIBLE if r in REFERENCE_NAMES]
F1_OUTDIR = Path(_F1_CFG.get("outdir", "results/f1_concordance"))
F1_PAIRS = [
    (a, b)
    for i, a in enumerate(_F1_COMPATIBLE)
    for b in _F1_COMPATIBLE[i + 1 :]
]

# Derive canonical chromosome names from the conspec .fai so we never need
# to maintain a manual list. These are the chromosomes present in all three
# references (augref just appends SV_* contigs on top of conspec).
_conspec_fai = Path(ALIGN_CONSPEC).with_suffix(".fai")
if not _conspec_fai.exists():
    _conspec_fai = Path(f"{ALIGN_CONSPEC}.fai")
if _conspec_fai.exists():
    ROH_CANONICAL_CHROMS = [
        line.split("\t")[0]
        for line in _conspec_fai.read_text().splitlines()
        if line.strip()
    ]
else:
    ROH_CANONICAL_CHROMS = []  # fai not yet built; --regions omitted, min_roh_length guards SV contigs

PCANGSD_CFG = config.get("pcangsd", {})
PCANGSD_OUTDIR = Path(PCANGSD_CFG.get("outdir", "results/pcangsd"))
PCANGSD_LD_WINDOW = PCANGSD_CFG.get("ld_window", 50)
PCANGSD_LD_STEP = PCANGSD_CFG.get("ld_step", 10)
PCANGSD_LD_R2 = PCANGSD_CFG.get("ld_r2", 0.3)
PCANGSD_MAF = PCANGSD_CFG.get("maf", 0.05)
PCANGSD_N_PCS = PCANGSD_CFG.get("n_pcs", 10)

PLOT_OUTDIR = Path(config.get("plot", {}).get("outdir", "results/plots"))


shell.executable("bash")

# Cap implicit BLAS / numerical-library thread fans at 1 for EVERY shell rule.
# Without this, anything that imports numpy/scipy/sklearn/matplotlib (vg, pcangsd,
# any python plotting rule) tries to spin up one OpenBLAS thread per host core
# and trips RLIMIT_NPROC inside the slurm cgroup. Rules that genuinely need
# parallel BLAS can still override these in their shell body. Applies globally
# so we stop playing whack-a-mole rule-by-rule.
shell.prefix(
    "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "
    "NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ; "
)

# Global wildcard constraints: prevent greedy matching of {ref} across underscores
# in paths like {ref}_{pop1}_vs_{pop2}_fst.png where both can contain underscores.
# sample is constrained globally to known sample IDs so wildcards in paths
# like results/variants/vg/chunks/{sample}/{contig}.vcf.gz don't greedily
# match sub-paths (e.g. {sample}=chunks/ERR.../NC_...).
_ALL_SAMPLES = sorted(set(SHORT_SAMPLES) | set(LONG_SAMPLES))
wildcard_constraints:
    ref="|".join(REFERENCE_NAMES),
    sample="|".join(_ALL_SAMPLES) if _ALL_SAMPLES else "x^",
    contig="[A-Za-z0-9._-]+"

# All references now have a real BAM — mc_graph uses a surjected BAM from
# the giraffe GAM so its SNPs go through bcftools like every other ref.
_BCFTOOLS_REFS = REFERENCE_NAMES

# Linear references only (everything except mc_graph). Used by rules whose
# inputs/behaviour assume a true bwa-mem BAM (e.g. ANGSD genotype likelihoods,
# samtools flagstat-based mapping metrics — surjection loss biases these for
# mc_graph so it gets its own GAM-derived or VCF-derived path instead).
_LINEAR_REFS = [r for r in REFERENCE_NAMES if r != "mc_graph"]

rule all:
    input:
        # Assembly and QC (long reads only)
        expand(ASSEMBLY_OUTDIR / "nanostat/{sample}_nanostat.txt", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "nanoplot/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "quast/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "busco/{sample}", sample=LONG_SAMPLES),
        ASSEMBLY_OUTDIR / "assembly_qc_summary.tsv",
        # SV calling
        SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.jasmine.vcf",
        SV_OUTDIR / "pan_sample_catalog/catalog_stats.txt",
        SV_OUTDIR / "pan_sample_catalog/sv_support_matrix.txt",
        SV_OUTDIR / "pan_sample_catalog/novel_sequence_summary.tsv",
        SV_OUTDIR / "pan_sample_catalog/sv_sharing_summary.tsv",
        SV_OUTDIR / "augref/augmented_reference.fasta",
        # Cactus graph (optional; INCLUDE_GRAPH)
        *([
            CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz",
            CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gfa.gz",
            CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vcf.gz",
        ] if INCLUDE_GRAPH else []),
        # WGS alignment to linear references (short reads only)
        expand(ALIGN_OUTDIR / "{sample}/{sample}.augref.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam.bai", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam.bai", sample=SHORT_SAMPLES),
        # Graph alignment: giraffe GAM + surjected mc_graph BAM (optional)
        *([
            *expand(ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam", sample=SHORT_SAMPLES),
            *expand(ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam", sample=SHORT_SAMPLES),
            *expand(ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam.bai", sample=SHORT_SAMPLES),
            *expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES),
            *expand(METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv", sample=SHORT_SAMPLES),
        ] if INCLUDE_GRAPH else []),
        # Alignment metrics (linear refs)
        expand(METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv", sample=SHORT_SAMPLES, ref=_BCFTOOLS_REFS),
        METRICS_OUTDIR / "alignment_metrics.tsv",
        # Short reads → individual assemblies QC (optional; INCLUDE_ASSEMBLY_ALIGN)
        *([
            *expand(METRICS_OUTDIR / "short_to_assembly/{short_sample}__{long_sample}.metrics.tsv",
                    short_sample=SHORT_SAMPLES, long_sample=LONG_SAMPLES),
            METRICS_OUTDIR / "short_to_assembly_metrics.tsv",
        ] if INCLUDE_ASSEMBLY_ALIGN else []),
        # SNP calling (bcftools on BAMs for all refs including mc_graph)
        expand(VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz.tbi", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz", ref=REFERENCE_NAMES),
        expand(VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi", ref=REFERENCE_NAMES),
        # SV calling (vg call on augmented graph, for benchmarking; optional)
        *([
            *expand(VC_OUTDIR / "sv/vg/{sample}.vcf.gz", sample=SHORT_SAMPLES),
            *expand(VC_OUTDIR / "sv/vg/{sample}.vcf.gz.tbi", sample=SHORT_SAMPLES),
        ] if INCLUDE_GRAPH else []),
        # Population genetics
        expand(FST_OUTDIR / "afs/{ref}.afs.tsv", ref=REFERENCE_NAMES),
        [FST_OUTDIR / f"fst/{ref}/{pop1}_vs_{pop2}.weir.fst" for ref in REFERENCE_NAMES for (pop1, pop2) in POP_PAIR_TUPLES],
        expand(PI_OUTDIR / "{ref}.windowed.pi", ref=REFERENCE_NAMES),
        expand(PI_OUTDIR / "{ref}.Tajima.D", ref=REFERENCE_NAMES),
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.raw.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        AB_OUTDIR / "allelic_balance_summary.tsv",
        # ROH removed (see note near ROH_CFG above):
        # expand(ROH_OUTDIR / "{ref}/{sample}.roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        # expand(ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        # ROH_OUTDIR / "f_roh_summary.tsv",
        # PCAngsd PCA + selection scan
        expand(PCANGSD_OUTDIR / "{ref}/pcangsd.cov", ref=REFERENCE_NAMES),
        expand(PCANGSD_OUTDIR / "{ref}/fst_outliers.tsv", ref=REFERENCE_NAMES),
        # Plots
        expand(PLOT_OUTDIR / "pca/{ref}_pca.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "selection/{ref}_selection_scan.png", ref=REFERENCE_NAMES),
        [PLOT_OUTDIR / f"fst/{ref}_{pop1}_vs_{pop2}_fst.png" for ref in REFERENCE_NAMES for (pop1, pop2) in POP_PAIR_TUPLES],
        expand(PLOT_OUTDIR / "afs/{ref}_afs.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "pi/{ref}_pi.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "tajimas_d/{ref}_tajimas_d.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.png", ref=REFERENCE_NAMES),
        # ROH removed: expand(PLOT_OUTDIR / "roh/{ref}_f_roh.png", ref=REFERENCE_NAMES),
        PLOT_OUTDIR / "alignment/alignment_rates.png",
        *([PLOT_OUTDIR / "alignment/assembly_alignment.png"] if INCLUDE_ASSEMBLY_ALIGN else []),
        PLOT_OUTDIR / "qc/coverage_qc.png",
        # Per-sample heterozygosity
        expand(HET_OUTDIR / "{ref}/per_sample_heterozygosity.tsv", ref=REFERENCE_NAMES),
        # Relatedness (per-ref + cohort summary)
        expand(RELATEDNESS_OUTDIR / "{ref}/relatedness.tsv", ref=REFERENCE_NAMES),
        RELATEDNESS_OUTDIR / "relatedness_summary.tsv",
        # Cohort QC and summary tables
        QC_OUTDIR / "sample_qc.tsv",
        QC_OUTDIR / "per_ref_summary.tsv",
        # F1 concordance between SNP call sets across refs
        F1_OUTDIR / "f1_summary.tsv"

rule nanostat:
    input:
        fastq=lambda wildcards: LONG_READS[wildcards.sample],
    output:
        stats=ASSEMBLY_OUTDIR / "nanostat/{sample}_nanostat.txt",
    conda:
        "envs/assembly.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=4000,
        cpus=4,
    shell:
        r"""
        set -euo pipefail
        NanoStat \
          --fastq {input.fastq} \
          --name {wildcards.sample}_nanostat.txt \
          --outdir {ASSEMBLY_OUTDIR}/nanostat \
          --threads {threads}
        """


rule nanoplot:
    input:
        fastq=lambda wildcards: LONG_READS[wildcards.sample],
    output:
        outdir=directory(ASSEMBLY_OUTDIR / "nanoplot/{sample}"),
    conda:
        "envs/assembly.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=4000,
        cpus=4,
    shell:
        r"""
        set -euo pipefail
        NanoPlot \
          --fastq {input.fastq} \
          --prefix {wildcards.sample} \
          --outdir {output.outdir} \
          --threads {threads}
        """


rule hifiasm_assemble:
    input:
        fastq=lambda wildcards: LONG_READS[wildcards.sample],
    output:
        assembly=ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta",
    conda:
        "envs/assembly.yaml"
    threads: ASSEMBLY_THREADS
    resources:
        slurm_partition="long",
        runtime=1440,
        mem_mb=32000,
        cpus=ASSEMBLY_THREADS,
    params:
        prefix=lambda wildcards: ASSEMBLY_OUTDIR / f"hifiasm/{wildcards.sample}/{wildcards.sample}",
        ont_flag=lambda wildcards: "--ont" if PLATFORM_MAP.get(wildcards.sample, "") == "ONT" else "",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        hifiasm \
          -o {params.prefix} \
          -t {threads} \
          --primary \
          -f0 \
          {params.ont_flag} \
          {input.fastq}
        awk '/^S/{{print ">"$2; print $3}}' {params.prefix}.p_ctg.gfa > {output.assembly}
        """


rule quast_assembly:
    input:
        assembly=ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta",
    output:
        outdir=directory(ASSEMBLY_OUTDIR / "quast/{sample}"),
    conda:
        "envs/assembly.yaml"
    threads: 8
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=8,
    shell:
        r"""
        set -euo pipefail
        quast \
          {input.assembly} \
          -t {threads} \
          -o {output.outdir}
        """


rule busco_assembly:
    input:
        assembly=ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta",
    output:
        outdir=directory(ASSEMBLY_OUTDIR / "busco/{sample}"),
    conda:
        "envs/assembly.yaml"
    threads: ASSEMBLY_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=32000,
        cpus=ASSEMBLY_THREADS,
    params:
        lineage=ASSEMBLY_LINEAGE,
    shell:
        r"""
        set -euo pipefail
        # BUSCO calls BBTools (stats.sh), which needs an explicit JVM heap on large assemblies.
        export JAVA_TOOL_OPTIONS="-Xmx28g"
        export _JAVA_OPTIONS="-Xmx28g"
        busco \
          -i {input.assembly} \
          -m genome \
          --lineage_dataset {params.lineage} \
          -c {threads} \
          -o {wildcards.sample} \
          --out_path {ASSEMBLY_OUTDIR}/busco \
          -f
        """

rule assembly_qc_summary:
    """Aggregate per-sample QUAST + BUSCO into a single TSV that's easy to
    drop into a paper supplement (N50, total length, BUSCO C/D/F/M, ...)."""
    input:
        quast=expand(ASSEMBLY_OUTDIR / "quast/{sample}", sample=LONG_SAMPLES),
        busco=expand(ASSEMBLY_OUTDIR / "busco/{sample}", sample=LONG_SAMPLES)
    output:
        summary=ASSEMBLY_OUTDIR / "assembly_qc_summary.tsv"
    conda:
        "envs/plotting.yaml"
    params:
        samples=LONG_SAMPLES,
        quast_dirs=[str(ASSEMBLY_OUTDIR / f"quast/{s}") for s in LONG_SAMPLES],
        busco_dirs=[str(ASSEMBLY_OUTDIR / f"busco/{s}") for s in LONG_SAMPLES],
        lineage=ASSEMBLY_LINEAGE
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1
    script:
        "scripts/assembly_qc_summary.py"


include: "rules/sv_calling.smk"

# ---------------------------------------------------------------------------
# Post-processing: Build augmented reference from pan-sample SV catalog
# ---------------------------------------------------------------------------
rule sv_extract_sequences:
    input:
        vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        reference=ALIGN_CONSPEC,
    output:
        fasta=SV_OUTDIR / "augref/extracted_flanked_sv_seqs.fasta",
    conda:
        "envs/sv_calling.yaml"
    params:
        flank=SV_FLANK,
        min_ins=SV_MIN_SIZE,
        samples=" ".join(LONG_SAMPLES),
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    script:
        "scripts/sv_extract_sequences.py"

rule sv_dedup_sequences:
    input:
        fasta=SV_OUTDIR / "augref/extracted_flanked_sv_seqs.fasta",
    output:
        fasta=SV_OUTDIR / "augref/extracted_flanked_sv_seqs.dedup.fasta",
    conda:
        "envs/sv_calling.yaml"
    shell:
        r"""
        set -euo pipefail
        cd-hit-est -i {input.fasta} -o {output.fasta} -c 0.95 -n 10
        """

rule sv_build_augmented_reference:
    input:
        reference=ALIGN_CONSPEC,
        sv_seqs=SV_OUTDIR / "augref/extracted_flanked_sv_seqs.dedup.fasta",
    output:
        augref=SV_OUTDIR / "augref/augmented_reference.fasta",
    shell:
        r"""
        set -euo pipefail
        cat {input.reference} {input.sv_seqs} > {output.augref}
        """

CACTUS_STAGED_REF = Path("results") / "cactus_work" / "reference.fasta"

rule stage_cactus_reference:
    """Stage the reference FASTA used to build the pangenome graph.

    This is the ONLY place the graph reference diverges from the linear pipeline:
    the linear tools keep using ALIGN_CONSPEC (the full reference) while the graph
    uses this filtered copy. Sequences named in cactus.exclude_ref_seqs are dropped
    here so circular/organellar contigs (e.g. the mitochondrion NC_001566.1) don't
    end up in the graph, where they form top-level chain loops that break
    `vg haplotypes` haplotype sampling. With no exclusions this is a plain copy."""
    input:
        ALIGN_CONSPEC
    output:
        CACTUS_STAGED_REF
    params:
        exclude=" ".join(CACTUS_EXCLUDE_REF_SEQS)
    container:
        VG_IMAGE
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.exclude}" ]; then
            cp {input} {output}
        else
            # Names actually present in the FASTA (first token after '>').
            present=$(grep '^>' {input} | sed 's/^>//; s/[[:space:]].*//')
            # Fail loudly if any requested exclusion isn't in the FASTA, rather
            # than silently producing an unfiltered graph reference.
            for s in {params.exclude}; do
                if ! printf '%s\n' "$present" | grep -qxF "$s"; then
                    echo "ERROR: cactus.exclude_ref_seqs name '$s' not found in {input}" >&2
                    echo "Available sequence names:" >&2
                    printf '%s\n' "$present" >&2
                    exit 1
                fi
            done
            # Build the keep-list (everything not excluded), preserving FASTA order.
            keep=$(printf '%s\n' "$present" | grep -vxF -f <(printf '%s\n' {params.exclude}) || true)
            if [ -z "$keep" ]; then
                echo "ERROR: cactus.exclude_ref_seqs would remove every sequence from {input}" >&2
                exit 1
            fi
            samtools faidx {input} $keep > {output}
        fi
        """

rule generate_cactus_seqfile:
    input:
        reference=CACTUS_STAGED_REF,
        assemblies=expand(ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta", sample=LONG_SAMPLES)
    output:
        seqfile=CACTUS_SEQFILE
    run:
        with open(output.seqfile, "w") as f:
            f.write(f"reference {input.reference}\n")
            for sample, assembly in zip(LONG_SAMPLES, input.assemblies):
                f.write(f"{sample} {assembly}\n")

rule make_cactus_graph:
    input:
        seqfile=CACTUS_SEQFILE,
        reference=ALIGN_CONSPEC,
        assemblies=expand(ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta", sample=LONG_SAMPLES)
    output:
        gbz=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz",
        gfa=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gfa.gz",
        vcf=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vcf.gz"
    container:
        CACTUS_IMAGE
    threads:
        CACTUS_MAX_CORES
    resources:
        slurm_partition="highmem",
        runtime=2880,
        mem_mb=128000,
        cpus=CACTUS_MAX_CORES
    params:
        jobstore=CACTUS_JOBSTORE,
        outdir=CACTUS_OUTDIR,
        outname=CACTUS_OUTNAME,
        ref_contigs=CACTUS_REF_CONTIGS,
        extra_args=CACTUS_EXTRA_ARGS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        rm -rf {params.jobstore}
        cactus-pangenome \
          {params.jobstore} \
          {input.seqfile} \
          --reference "reference" \
          --outDir {params.outdir} \
          --outName "{params.outname}" \
          --vcf \
          --gfa \
          --gbz \
          --maxCores {threads} \
          {params.ref_contigs} \
          {params.extra_args}
        if [ -f "{params.outdir}/{params.outname}.gfa" ] && [ ! -f "{params.outdir}/{params.outname}.gfa.gz" ]; then
            gzip -f "{params.outdir}/{params.outname}.gfa"
        fi
        """

rule make_haplo_index:
    """Build the haplotype-sampling indexes (.dist, .ri, .hapl) used by
    `vg giraffe -H`. These are normally produced by `cactus-pangenome --haplo`,
    but the cactus container's `vg haplotypes` aborts with "cannot run N threads
    in parallel on this system" on some nodes. We reproduce the index sequence
    here in the standalone VG_IMAGE, running `vg haplotypes` single-threaded to
    avoid that parallelism check. The .gbz is unchanged from cactus.

    The distance index is built with `-P "reference"` so the reference sample is
    used as the backbone when orienting top-level chains. Without it, vg
    haplotypes fails with "top-level chain N is a loop; haplotype sampling cannot
    be used with this graph". "reference" matches `--reference` in
    make_cactus_graph (the PanSN sample name of the reference assembly).

    We also rebuild the short-read minimizer/zipcodes index here. `vg giraffe`
    auto-discovers sibling indexes by name (.dist, .shortread.withzip.min,
    .shortread.zipcodes) and refuses to run if .dist is newer than the minimizer
    index derived from it. Since we regenerate .dist, the minimizer index must be
    rebuilt from the same .dist or giraffe aborts with a staleness error. Building
    them together here keeps the whole giraffe index set mutually consistent.

    Note the .dist is a FULL distance index (no --no-nested-distance). vg
    haplotypes alone would tolerate the reduced top-level-only index that cactus
    uses, but vg minimizer needs nested-snarl distances, so we build one full
    index and share it across both steps."""
    input:
        gbz=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz"
    output:
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.dist",
        ri=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.ri",
        hapl=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.hapl",
        min=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.shortread.withzip.min",
        zipcodes=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.shortread.zipcodes"
    container:
        VG_IMAGE
    threads:
        CACTUS_MAX_CORES
    resources:
        slurm_partition="highmem",
        runtime=480,
        mem_mb=32000,
        cpus=CACTUS_MAX_CORES
    shell:
        r"""
        set -euo pipefail
        # Full distance index (NO --no-nested-distance). vg haplotypes only needs
        # top-level chain distances, but vg minimizer's zipcode/payload caching
        # requires nested-snarl distances ("is_regular_snarl requires distances in
        # the distance index"). A full index satisfies both, so we build one .dist
        # and share it; the reduced (--no-nested-distance) index would crash vg
        # minimizer.
        vg index -t {threads} -j {output.dist} -P "reference" {input.gbz}
        vg gbwt -p --num-threads {threads} -r {output.ri} -Z {input.gbz}
        vg haplotypes -v 2 -t 1 -H {output.hapl} -d {output.dist} -r {output.ri} {input.gbz}
        vg minimizer -t {threads} -d {output.dist} \
          -o {output.min} -z {output.zipcodes} {input.gbz}
        """

rule bwa_index:
    input:
        fasta="{fasta}"
    output:
        bwt="{fasta}.bwt",
        amb="{fasta}.amb",
        ann="{fasta}.ann",
        pac="{fasta}.pac",
        sa="{fasta}.sa"
    conda:
        "envs/align_wgs.yaml"
    shell:
        r"""
        set -euo pipefail
        bwa index {input.fasta}
        """

rule downsample_short_reads:
    """Subsample paired-end FASTQs to ALIGN_DOWNSAMPLE_BASES (default 10 G)
    of uncompressed sequence so that oversized libraries don't bottleneck
    alignment and variant calling.  rasusa passes all reads through unchanged
    when the input is already smaller than the target."""
    input:
        r1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        r2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
    output:
        r1=ALIGN_OUTDIR / "{sample}/{sample}.R1.ds.fastq.gz",
        r2=ALIGN_OUTDIR / "{sample}/{sample}.R2.ds.fastq.gz",
    conda:
        "envs/align_wgs.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        rasusa reads \
            -b {ALIGN_DOWNSAMPLE_BASES} \
            -o {output.r1} -o {output.r2} \
            {input.r1} {input.r2}
        """

rule make_augref_alt:
    """Build the BWA .alt sidecar that tells bwa-mem the SV_* contigs are
    alternative representations of conspec breakpoints. Generated by
    aligning the deduplicated SV insertion sequences back to conspec — the
    raw bwa-mem SAM output IS the .alt file format. With this present next
    to the augref index, bwa-mem stops collapsing the shared 200 bp flanks
    to MAPQ 0 and assigns proper MAPQ across SV breakpoints, which propagates
    through bcftools and the popgen stack with no other changes."""
    input:
        conspec=ALIGN_CONSPEC,
        conspec_bwt=f"{ALIGN_CONSPEC}.bwt",
        conspec_amb=f"{ALIGN_CONSPEC}.amb",
        conspec_ann=f"{ALIGN_CONSPEC}.ann",
        conspec_pac=f"{ALIGN_CONSPEC}.pac",
        conspec_sa=f"{ALIGN_CONSPEC}.sa",
        sv_seqs=SV_OUTDIR / "augref/extracted_flanked_sv_seqs.dedup.fasta",
        # Augref must exist (and be indexed) before alignments consume the
        # .alt file — it ends up next to the augref index by naming convention.
        augref=ALIGN_AUGREF,
        augref_bwt=f"{ALIGN_AUGREF}.bwt"
    output:
        alt=f"{ALIGN_AUGREF}.alt"
    conda:
        "envs/align_wgs.yaml"
    threads: 8
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=16000,
        cpus=8
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} {input.conspec} {input.sv_seqs} > {output.alt}
        """


rule align_bwa_augref:
    input:
        fq1=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R1.ds.fastq.gz"),
        fq2=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R2.ds.fastq.gz"),
        ref=ALIGN_AUGREF,
        bwt=f"{ALIGN_AUGREF}.bwt",
        alt=ancient(f"{ALIGN_AUGREF}.alt")
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai"
    conda:
        "envs/align_wgs.yaml"
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=960,
        mem_mb=32000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        RG="@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA"
        bwa mem -k 17 -O 5,5 -E 2,2 -Y -R "$RG" -t {threads} {input.ref} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule align_bwa_conspec:
    input:
        fq1=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R1.ds.fastq.gz"),
        fq2=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R2.ds.fastq.gz"),
        ref=ALIGN_CONSPEC,
        bwt=f"{ALIGN_CONSPEC}.bwt"
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam.bai"
    conda:
        "envs/align_wgs.yaml"
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=960,
        mem_mb=32000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        RG="@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA"
        bwa mem -k 17 -O 5,5 -E 2,2 -Y -R "$RG" -t {threads} {input.ref} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule align_bwa_hetspec:
    input:
        fq1=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R1.ds.fastq.gz"),
        fq2=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R2.ds.fastq.gz"),
        ref=ALIGN_HETSPEC,
        bwt=f"{ALIGN_HETSPEC}.bwt"
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam.bai"
    conda:
        "envs/align_wgs.yaml"
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=960,
        mem_mb=32000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        RG="@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA"
        bwa mem -k 17 -O 5,5 -E 2,2 -Y -R "$RG" -t {threads} {input.ref} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule giraffe_align:
    input:
        fq1=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R1.ds.fastq.gz"),
        fq2=ancient(ALIGN_OUTDIR / "{sample}/{sample}.R2.ds.fastq.gz"),
        gbz=ALIGN_CACTUS_GBZ,
        hapl=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.hapl",
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.dist",
        min=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.shortread.withzip.min",
        zipcodes=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.shortread.zipcodes"
    output:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    container:
        VG_IMAGE
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=960,
        mem_mb=32000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        export OPENBLAS_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        vg giraffe -Z {input.gbz} -H {input.hapl} \
          -t {threads} -f {input.fq1} -f {input.fq2} -p \
          --rescue-attempts 5 --sample {wildcards.sample} > {output.gam}
        """

rule vg_surject:
    """Project the giraffe GAM onto the conspec reference paths to produce a
    coordinate-sorted BAM in conspec coordinates. Used for bcftools SNP calling
    on mc_graph so SNP calls are on the same coordinates as the linear refs.

    -F: restrict surjection to PanSN paths whose final component matches a
        conspec contig (so surjection targets the linear reference paths only,
        not assembly haplotypes or other paths in the GBZ).
    --prune-low-cplx: drop low-complexity alignments that inflate false SNPs.

    After surjection, CHROM names are still PanSN-prefixed (e.g.
    "reference#0#NC_001234.1"). samtools reheader strips the prefix so the BAM
    matches the conspec FASTA used by bcftools mpileup downstream.

    NOTE: A small fraction of reads in the GAM align to non-reference graph
    paths only (e.g. SV alt paths absent from conspec). vg surject drops those;
    they appear in the SV pipeline but not in the SNP pipeline. This is the
    surjection loss the user warned about — it's expected and unavoidable for
    SNP calling on linear coordinates."""
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam",
        gbz=VG_GBZ,
        conspec_fai=f"{ALIGN_CONSPEC}.fai",
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam.bai",
    container:
        VG_IMAGE
    threads: ALIGN_THREADS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=16000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        export OMP_NUM_THREADS={threads}
        mkdir -p $(dirname {output.bam})
        _workdir=$(dirname {output.bam})

        # Resolve PanSN path names in the GBZ that correspond to conspec
        # contigs. The .fai gives us bare contig names; vg paths -L gives us
        # whatever's in the GBZ. Match the last "#"-delimited field, or fall
        # back to exact-match for non-PanSN graphs.
        contigs_tmp=$(mktemp -p "$_workdir")
        paths_tmp=$(mktemp -p "$_workdir")
        path_list=$(mktemp -p "$_workdir")
        rename_map=$(mktemp -p "$_workdir")
        trap "rm -f $contigs_tmp $paths_tmp $path_list $rename_map" EXIT

        awk -F'\t' '{{print $1}}' {input.conspec_fai} > "$contigs_tmp"
        vg paths -x {input.gbz} -L > "$paths_tmp"

        awk -v contigs="$contigs_tmp" '
            BEGIN {{
                while ((getline c < contigs) > 0) want[c] = 1
            }}
            {{
                n = split($0, parts, "#")
                last = (n >= 3) ? parts[n] : $0
                if (last in want) print $0
            }}
        ' "$paths_tmp" > "$path_list"

        if [ ! -s "$path_list" ]; then
            echo "ERROR: no GBZ paths matched any conspec contig" >&2
            exit 1
        fi

        # Build a CHROM rename map (PanSN path → bare contig name) for
        # samtools reheader. Same matching logic as path_list.
        awk 'BEGIN{{FS="#"}} {{
            if (NF >= 3) print $0 "\t" $NF
            else print $0 "\t" $0
        }}' "$path_list" > "$rename_map"

        vg surject -x {input.gbz} {input.gam} \
            -b --prune-low-cplx \
            -F "$path_list" \
            -t {threads} \
          | samtools sort -n -@ {threads} \
          | samtools fixmate -m -@ {threads} - - \
          | samtools sort -@ {threads} -O bam -o "$_workdir/{wildcards.sample}.surject.tmp.bam"

        # Strip PanSN prefixes from @SQ headers so CHROM matches conspec.
        # samtools reheader only rewrites the header, so we edit @SQ SN: in
        # place using the rename map.
        samtools view -H "$_workdir/{wildcards.sample}.surject.tmp.bam" \
          | awk -v map="$rename_map" '
                BEGIN {{
                    while ((getline line < map) > 0) {{
                        split(line, a, "\t")
                        rn[a[1]] = a[2]
                    }}
                }}
                /^@SQ/ {{
                    n = split($0, fields, "\t")
                    for (i = 1; i <= n; i++) {{
                        if (fields[i] ~ /^SN:/) {{
                            old = substr(fields[i], 4)
                            if (old in rn) fields[i] = "SN:" rn[old]
                        }}
                    }}
                    line = fields[1]
                    for (i = 2; i <= n; i++) line = line "\t" fields[i]
                    print line
                    next
                }}
                {{ print }}
            ' > "$_workdir/{wildcards.sample}.reheader.sam"

        samtools reheader "$_workdir/{wildcards.sample}.reheader.sam" \
            "$_workdir/{wildcards.sample}.surject.tmp.bam" > {output.bam}
        samtools index -@ {threads} {output.bam}

        rm -f "$_workdir/{wildcards.sample}.surject.tmp.bam" \
              "$_workdir/{wildcards.sample}.reheader.sam"
        """


rule count_short_reads:
    """Derive the canonical per-sample read count from the conspec BAM.
    bwa-mem writes every input read to the BAM (mapped or not), so
    `samtools view -c` on the conspec BAM equals the downsampled-fastq
    read count exactly — without touching the temp fastq files, which
    would trigger a full mtime cascade if they needed to be recreated."""
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam.bai",
    output:
        txt=ALIGN_OUTDIR / "{sample}/{sample}.read_count.txt"
    conda:
        "envs/align_metrics.yaml"
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1
    wildcard_constraints:
        sample="|".join(SHORT_SAMPLES) if SHORT_SAMPLES else "x^"
    shell:
        r"""
        set -euo pipefail
        samtools view -c {input.bam} > {output.txt}
        """


rule align_metrics_per_bam:
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam.bai",
        read_count=ALIGN_OUTDIR / "{sample}/{sample}.read_count.txt"
    output:
        METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv"
    conda:
        "envs/align_metrics.yaml"
    threads:
        METRICS_THREADS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=METRICS_THREADS
    params:
        min_mapq=METRICS_MIN_MAPQ,
        min_baseq=METRICS_MIN_BASEQ
    wildcard_constraints:
        # mc_graph BAM metrics come from align_metrics_per_gam (vg stats on
        # the GAM) rather than samtools flagstat — surjection loss would make
        # mapping rate look artificially low.
        ref="|".join(_LINEAR_REFS) if _LINEAR_REFS else "x^"
    shell:
        r"""
        set -euo pipefail
        # vg pulls in OpenBLAS which tries to spawn one thread per core and
        # blows past RLIMIT_NPROC on the slurm short partition. Pin to 1.
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        mkdir -p {METRICS_OUTDIR}/{wildcards.sample}

        # Denominator is the fastq read count, not the BAM total, so the
        # mapping rate denominator is consistent with how mc_graph reports
        # it (vg stats Total alignments).
        total=$(cat {input.read_count})
        mapped=$(samtools view -c -F 0x904 {input.bam})
        mean_mapq=$(samtools view -F 0x904 {input.bam} | awk '{{sum+=$5; n++}} END {{if(n>0) printf "%.6f", sum/n; else print "0"}}')
        mean_depth=$(samtools depth -a -q {params.min_mapq} -Q {params.min_baseq} {input.bam} | \
            awk '{{sum+=$3; n++}} END {{if (n>0) printf "%.6f", sum/n; else print "0"}}')
        map_rate=$(awk -v m="$mapped" -v t="$total" 'BEGIN {{if (t>0) printf "%.6f", m/t; else print "0"}}')

        cat > {output} << EOF
sample	alignment_type	total_reads	aligned_reads	mapping_rate	mean_mapq	mean_depth
{wildcards.sample}	{wildcards.ref}	$total	$mapped	$map_rate	$mean_mapq	$mean_depth
EOF
        """

rule align_metrics_per_gam:
    """Emit GAM-derived metrics under both 'cactus' (raw giraffe) and
    'mc_graph' (the canonical graph-reference label used downstream).
    Uses vg stats on the GAM directly — surjected BAM metrics are not used
    here since surjection loss would make mc_graph look artificially worse."""
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    output:
        cactus=METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv",
        mc_graph=METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv"
    container:
        VG_IMAGE
    threads:
        METRICS_THREADS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=METRICS_THREADS
    shell:
        r"""
        set -euo pipefail
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        mkdir -p {METRICS_OUTDIR}/{wildcards.sample}

        stats=$(vg stats -a {input.gam})
        total_reads=$(echo "$stats" | awk '/Total alignments:/ {{print $3}}')
        aligned_reads=$(echo "$stats" | awk '/Total aligned:/ {{print $3}}')
        mean_mapq=$(echo "$stats" | awk '/Mapping quality:/ {{print $5}}' | tr -d ',')
        map_rate=$(awk -v m="$aligned_reads" -v t="$total_reads" \
            'BEGIN {{if (t>0) printf "%.6f", m/t; else print "0"}}')

        cat > {output.cactus} << EOF
sample	alignment_type	total_reads	aligned_reads	mapping_rate	mean_mapq	mean_depth
{wildcards.sample}	cactus	$total_reads	$aligned_reads	$map_rate	$mean_mapq	NA
EOF

        cat > {output.mc_graph} << EOF
sample	alignment_type	total_reads	aligned_reads	mapping_rate	mean_mapq	mean_depth
{wildcards.sample}	mc_graph	$total_reads	$aligned_reads	$map_rate	$mean_mapq	NA
EOF
        """

rule align_short_to_assembly:
    """Map each short-read sample to each long-read assembly and emit a one-line
    metrics TSV. We stream minimap2 → samtools flagstat with no on-disk BAM,
    since these alignments are only used to score how well an assembly
    represents the short-read population. This intentionally avoids the
    bwa-based heavy alignment used for variant calling."""
    input:
        fq1=ancient(ALIGN_OUTDIR / "{short_sample}/{short_sample}.R1.ds.fastq.gz"),
        fq2=ancient(ALIGN_OUTDIR / "{short_sample}/{short_sample}.R2.ds.fastq.gz"),
        assembly=ASSEMBLY_OUTDIR / "hifiasm/{long_sample}/assembly.fasta"
    output:
        metrics=METRICS_OUTDIR / "short_to_assembly/{short_sample}__{long_sample}.metrics.tsv"
    conda:
        "envs/align_wgs.yaml"
    threads: ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=16000,
        cpus=ALIGN_THREADS
    wildcard_constraints:
        short_sample="|".join(SHORT_SAMPLES) if SHORT_SAMPLES else "x^",
        long_sample="|".join(LONG_SAMPLES) if LONG_SAMPLES else "x^"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {METRICS_OUTDIR}/short_to_assembly
        flagstat=$(minimap2 -ax sr -t {threads} {input.assembly} {input.fq1} {input.fq2} 2>/dev/null \
                   | samtools flagstat -)
        total=$(echo "$flagstat" | awk '/in total/ {{print $1; exit}}')
        primary=$(echo "$flagstat" | awk '/primary$/ {{print $1; exit}}')
        mapped=$(echo "$flagstat" | awk '/primary mapped/ {{print $1; exit}}')
        # Some samtools versions report "mapped" without "primary mapped"; fall back.
        if [ -z "$mapped" ]; then
            mapped=$(echo "$flagstat" | awk '/[0-9]+ mapped \(/ {{print $1; exit}}')
        fi
        denom=${{primary:-$total}}
        rate=$(awk -v m="$mapped" -v t="$denom" 'BEGIN {{if (t>0) printf "%.6f", m/t; else print "0"}}')

        cat > {output.metrics} << EOF
short_sample	long_sample	total_reads	mapped_reads	mapping_rate
{wildcards.short_sample}	{wildcards.long_sample}	$denom	$mapped	$rate
EOF
        """


rule align_short_to_assembly_summary:
    input:
        expand(METRICS_OUTDIR / "short_to_assembly/{short_sample}__{long_sample}.metrics.tsv",
               short_sample=SHORT_SAMPLES, long_sample=LONG_SAMPLES)
    output:
        metrics=METRICS_OUTDIR / "short_to_assembly_metrics.tsv"
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1
    run:
        from pathlib import Path
        out_path = Path(output.metrics)
        rows, header = [], None
        for path in input:
            text = Path(path).read_text().strip().splitlines()
            if not text:
                continue
            if header is None:
                header = text[0]
            if len(text) > 1:
                rows.append(text[1])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join([header] + rows) + "\n")


rule align_metrics_summary:
    input:
        expand(METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv",
               sample=SHORT_SAMPLES, ref=_BCFTOOLS_REFS),
        # Graph alignment metrics only when the graph is enabled.
        *([
            *expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES),
            *expand(METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv", sample=SHORT_SAMPLES),
        ] if INCLUDE_GRAPH else []),
    output:
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv"
    conda:
        "envs/align_metrics.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    run:
        from pathlib import Path
        
        out_path = Path(output.metrics)
        rows = []
        header = None
        for path in input:
            text = Path(path).read_text().strip().splitlines()
            if not text:
                continue
            if header is None:
                header = text[0]
            if len(text) > 1:
                rows.append(text[1])
        
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join([header] + rows) + "\n")

def _ref_fasta(wildcards):
    if wildcards.ref not in REFS_NESTED:
        raise ValueError(f"Unknown reference: {wildcards.ref}")
    return REFS_NESTED[wildcards.ref]["fasta"]

def _ref_fai(wildcards):
    return f"{_ref_fasta(wildcards)}.fai"

def _ref_dict(wildcards):
    return str(Path(_ref_fasta(wildcards)).with_suffix(".dict"))

rule index_augref:
    input:
        fasta=ALIGN_AUGREF
    output:
        fai=f"{ALIGN_AUGREF}.fai",
        dict=str(Path(ALIGN_AUGREF).with_suffix(".dict")),
        bwt=f"{ALIGN_AUGREF}.bwt",
        amb=f"{ALIGN_AUGREF}.amb",
        ann=f"{ALIGN_AUGREF}.ann",
        pac=f"{ALIGN_AUGREF}.pac",
        sa=f"{ALIGN_AUGREF}.sa"
    conda:
        "envs/index.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        samtools faidx {input.fasta}
        samtools dict {input.fasta} -o {output.dict}
        bwa index {input.fasta}
        """

rule index_conspec:
    input:
        fasta=ALIGN_CONSPEC
    output:
        fai=f"{ALIGN_CONSPEC}.fai",
        dict=str(Path(ALIGN_CONSPEC).with_suffix(".dict")),
        bwt=f"{ALIGN_CONSPEC}.bwt",
        amb=f"{ALIGN_CONSPEC}.amb",
        ann=f"{ALIGN_CONSPEC}.ann",
        pac=f"{ALIGN_CONSPEC}.pac",
        sa=f"{ALIGN_CONSPEC}.sa"
    conda:
        "envs/index.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        samtools faidx {input.fasta}
        samtools dict {input.fasta} -o {output.dict}
        bwa index {input.fasta}
        """

rule index_hetspec:
    input:
        fasta=ALIGN_HETSPEC
    output:
        fai=f"{ALIGN_HETSPEC}.fai",
        dict=str(Path(ALIGN_HETSPEC).with_suffix(".dict")),
        bwt=f"{ALIGN_HETSPEC}.bwt",
        amb=f"{ALIGN_HETSPEC}.amb",
        ann=f"{ALIGN_HETSPEC}.ann",
        pac=f"{ALIGN_HETSPEC}.pac",
        sa=f"{ALIGN_HETSPEC}.sa"
    conda:
        "envs/index.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        samtools faidx {input.fasta}
        samtools dict {input.fasta} -o {output.dict}
        bwa index {input.fasta}
        """

def _bams_for_ref(wildcards):
    return [str(ALIGN_OUTDIR / sample / f"{sample}.{wildcards.ref}.bam") for sample in SHORT_SAMPLES]

def _bais_for_ref(wildcards):
    return [str(ALIGN_OUTDIR / sample / f"{sample}.{wildcards.ref}.bam.bai") for sample in SHORT_SAMPLES]

def _ref_fasta_for_bcftools(wildcards):
    """Return the FASTA for bcftools mpileup. mc_graph surjects onto conspec,
    so it shares the conspec FASTA for pileup."""
    if wildcards.ref == "mc_graph":
        return ALIGN_CONSPEC
    return _ref_fasta(wildcards)

def _ref_fai_for_bcftools(wildcards):
    if wildcards.ref == "mc_graph":
        return f"{ALIGN_CONSPEC}.fai"
    return _ref_fai(wildcards)


rule bcftools_joint_call:
    """Joint SNP calling across all samples for one reference. mc_graph uses
    surjected BAMs (conspec coordinates) so it runs through the same pipeline
    as the linear refs and produces directly comparable SNP calls."""
    input:
        bams=_bams_for_ref,
        bais=_bais_for_ref,
        fasta=lambda wildcards: _ref_fasta_for_bcftools(wildcards),
        fai=lambda wildcards: _ref_fai_for_bcftools(wildcards),
    output:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    threads:
        BCFTOOLS_THREADS
    params:
        mpileup_extra=lambda wc: "-A" if wc.ref == "mc_graph" else "",
    resources:
        slurm_partition="long",
        runtime=2880,
        mem_mb=8000,
        cpus=BCFTOOLS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/bcftools/{wildcards.ref}/combined
        bam_list=$(mktemp -p {VC_OUTDIR}/bcftools/{wildcards.ref}/combined)
        tmp_vcf=$(mktemp -p {VC_OUTDIR}/bcftools/{wildcards.ref}/combined --suffix=.vcf.gz)
        trap "rm -f $bam_list $tmp_vcf" EXIT
        printf '%s\n' {input.bams} > "$bam_list"
        bcftools mpileup \
          {params.mpileup_extra} \
          -f {input.fasta} \
          -b "$bam_list" \
          -a AD,DP \
          -q 20 -Q 20 \
          --threads {threads} \
          2> {output.vcf}.mpileup.log \
        | bcftools call \
          -m \
          -v \
          -a GQ \
          --threads {threads} \
          -Oz \
          -o "$tmp_vcf"

        # Strip paths from sample names: "results/.../SRR123.mc_graph.bam" → "SRR123"
        bcftools query -l "$tmp_vcf" \
          | sed 's|.*/||; s|\..*||' \
          > "$tmp_vcf.samples"
        bcftools reheader -s "$tmp_vcf.samples" "$tmp_vcf" \
          -o {output.vcf}
        rm -f "$tmp_vcf.samples"
        bcftools index -t {output.vcf}
        """

# ============================================================================
# Graph SV genotyping (mc_graph only, for benchmarking/comparison)
# Genotypes SVs already present in the Cactus pangenome graph by packing
# per-sample GAMs directly against the GBZ and calling per-sample snarls.
# vg wiki confirms de novo SV calling via augment does not work; genotyping
# existing graph structure is the supported workflow.
# Output goes to a separate sv/vg/ path — not consumed by the popgen stack.
# ============================================================================

rule vg_sv_snarls:
    """Compute snarls on the pangenome graph once, shared across all samples."""
    input:
        gbz=VG_GBZ,
    output:
        snarls=VC_OUTDIR / "sv/vg/graph.snarls",
    container:
        VG_IMAGE
    threads: VG_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=64000,
        cpus=VG_THREADS
    shell:
        r"""
        set -euo pipefail
        export OMP_NUM_THREADS={threads}
        mkdir -p $(dirname {output.snarls})
        vg snarls -t {threads} {input.gbz} > {output.snarls}
        """


rule vg_sv_pack:
    """Compute per-sample read support against the pangenome graph."""
    input:
        gbz=VG_GBZ,
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam",
    output:
        pack=VC_OUTDIR / "sv/vg/{sample}.pack",
    container:
        VG_IMAGE
    threads: VG_THREADS
    resources:
        slurm_partition="long",
        runtime=720,
        mem_mb=32000,
        cpus=VG_THREADS
    shell:
        r"""
        set -euo pipefail
        export OMP_NUM_THREADS={threads}
        mkdir -p $(dirname {output.pack})
        vg pack -t {threads} -Q 5 -x {input.gbz} -g {input.gam} -o {output.pack}
        """


rule vg_sv_call_raw:
    """Genotype SVs per sample from the pangenome graph.
    -c 50 / -C 100000: only snarls with traversals between 50 bp and 100 kb
    (SVs only; skips SNP-scale snarls and the pathologically large snarls
    that cause multi-day runtimes).
    -A: call all snarls including nested (Jeon et al. 2026). We don't pass -a
    (which would emit 0/0 reference calls at every snarl); for the variant
    counts and SV-based popgen we only need variant sites."""
    input:
        gbz=VG_GBZ,
        snarls=VC_OUTDIR / "sv/vg/graph.snarls",
        pack=VC_OUTDIR / "sv/vg/{sample}.pack",
    output:
        vcf=temp(VC_OUTDIR / "sv/vg/{sample}.raw.vcf"),
    container:
        VG_IMAGE
    threads: VG_THREADS
    resources:
        slurm_partition="long",
        runtime=2880,
        mem_mb=16000,
        cpus=VG_THREADS
    shell:
        r"""
        set -euo pipefail
        export OMP_NUM_THREADS={threads}
        mkdir -p $(dirname {output.vcf})
        vg call -t {threads} \
            -k {input.pack} \
            -r {input.snarls} \
            -s {wildcards.sample} \
            -A \
            -c 50 -C 100000 \
            {input.gbz} \
          > {output.vcf}
        """


rule vg_sv_call_postprocess:
    """Strip PanSN prefixes from CHROM, compress, and index."""
    input:
        vcf=VC_OUTDIR / "sv/vg/{sample}.raw.vcf",
    output:
        vcf=VC_OUTDIR / "sv/vg/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "sv/vg/{sample}.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        _outdir=$(dirname {output.vcf})

        # Strip PanSN prefixes (sample#hap#contig → contig) in headers and data.
        rename_map=$(mktemp -p "$_outdir")
        trap "rm -f $rename_map" EXIT
        grep '^##contig=<ID=' {input.vcf} \
            | sed 's/##contig=<ID=//;s/,.*//' \
            | while read -r chrom; do
                n=$(echo "$chrom" | awk -F'#' '{{print NF}}')
                if [ "$n" -ge 3 ]; then
                    plain=$(echo "$chrom" | cut -d'#' -f3)
                    echo "$chrom $plain"
                fi
            done > "$rename_map"

        if [ -s "$rename_map" ]; then
            bcftools annotate --rename-chrs "$rename_map" -Oz -o {output.vcf} {input.vcf}
        else
            bgzip -c {input.vcf} > {output.vcf}
        fi
        tabix -p vcf {output.vcf}
        """

rule split_vcf_per_sample:
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        vcf=VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    wildcard_constraints:
        sample="|".join(SHORT_SAMPLES) if SHORT_SAMPLES else "x^"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/bcftools/{wildcards.ref}/per_sample
        bcftools view -s {wildcards.sample} -v snps -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """

include: "rules/popgen.smk"
include: "rules/plotting.smk"
