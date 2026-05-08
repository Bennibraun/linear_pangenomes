from pathlib import Path
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

# mc_graph is the Minigraph-Cactus pangenome reference. Variants come from
# vg call on the graph (chunked + per-sample-merged in vg_merge_samples),
# written to the same bcftools/{ref}/combined/merged.vcf.gz path as the
# linear refs so {ref}-wildcarded popgen rules (FST, AFS, π, allelic
# balance, ROH, PCAngsd) pick it up automatically. It shares the conspec
# FASTA only for indexing (the graph variants live in conspec coordinates
# because vg call linearizes against the reference path).
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

ROH_CFG = config["roh"]
ROH_OUTDIR = Path(ROH_CFG.get("outdir", "results/roh"))
ROH_GENOME_LENGTHS = ROH_CFG.get("genome_lengths", {})
ROH_MIN_LENGTH = ROH_CFG.get("min_roh_length", 100000)
ROH_BCFTOOLS_ARGS = ROH_CFG.get("bcftools_args", "")

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

# Dynamically bind any absolute reference paths into singularity
_bind_paths = set()
for ref in REFS_NESTED.values():
    p = Path(ref["fasta"])
    if p.is_absolute():
        _bind_paths.add(str(p.parent))

singularity_args = ("--bind " + ",".join(_bind_paths)) if _bind_paths else ""


shell.executable("bash")

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

# References that have a real BAM (everything except mc_graph, which now
# routes through vg call instead of a surjected BAM). Used by metrics and
# bcftools rules that must not match mc_graph.
_BCFTOOLS_REFS = [r for r in REFERENCE_NAMES if r != "mc_graph"]

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
        # Cactus graph
        CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz",
        CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gfa.gz",
        CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vcf.gz",
        # WGS alignment (short reads only)
        expand(ALIGN_OUTDIR / "{sample}/{sample}.augref.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam.bai", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam.bai", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam", sample=SHORT_SAMPLES),
        # Alignment metrics — mc_graph uses GAM-derived stats, not a BAM
        expand(METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv", sample=SHORT_SAMPLES, ref=_BCFTOOLS_REFS),
        expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES),
        expand(METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv", sample=SHORT_SAMPLES),
        METRICS_OUTDIR / "alignment_metrics.tsv",
        # Short reads → individual assemblies (one row per short × long pair)
        expand(METRICS_OUTDIR / "short_to_assembly/{short_sample}__{long_sample}.metrics.tsv",
               short_sample=SHORT_SAMPLES, long_sample=LONG_SAMPLES),
        METRICS_OUTDIR / "short_to_assembly_metrics.tsv",
        # Variant calling
        expand(VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz.tbi", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "vg/{sample}.vcf.gz", sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "vg/{sample}.vcf.gz.tbi", sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz", ref=REFERENCE_NAMES),
        expand(VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi", ref=REFERENCE_NAMES),
        # Population genetics
        expand(FST_OUTDIR / "afs/{ref}.afs.tsv", ref=REFERENCE_NAMES),
        [FST_OUTDIR / f"fst/{ref}/{pop1}_vs_{pop2}.weir.fst" for ref in REFERENCE_NAMES for (pop1, pop2) in POP_PAIR_TUPLES],
        expand(PI_OUTDIR / "{ref}.windowed.pi", ref=REFERENCE_NAMES),
        expand(PI_OUTDIR / "{ref}.Tajima.D", ref=REFERENCE_NAMES),
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.raw.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        AB_OUTDIR / "allelic_balance_summary.tsv",
        expand(ROH_OUTDIR / "{ref}/{sample}.roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        ROH_OUTDIR / "f_roh_summary.tsv",
        # PCAngsd PCA + selection scan
        expand(PCANGSD_OUTDIR / "{ref}/pcangsd.cov", ref=REFERENCE_NAMES),
        expand(PCANGSD_OUTDIR / "{ref}/fst_outliers.tsv", ref=REFERENCE_NAMES),
        # Plots
        expand(PLOT_OUTDIR / "pca/{ref}_pca.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "selection/{ref}_selection_scan.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "fst/{ref}_{pop1}_vs_{pop2}_fst.png",
               ref=REFERENCE_NAMES,
               pop1=[p[0] for p in POP_PAIR_TUPLES],
               pop2=[p[1] for p in POP_PAIR_TUPLES]),
        expand(PLOT_OUTDIR / "afs/{ref}_afs.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "pi/{ref}_pi.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "tajimas_d/{ref}_tajimas_d.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "roh/{ref}_f_roh.png", ref=REFERENCE_NAMES),
        PLOT_OUTDIR / "alignment/alignment_rates.png",
        PLOT_OUTDIR / "alignment/assembly_alignment.png"

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
    input:
        ALIGN_CONSPEC
    output:
        CACTUS_STAGED_REF
    shell:
        "cp {input} {output}"

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
        vcf=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vcf.gz",
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.d2.dist"
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
          --collapse \
          --outDir {params.outdir} \
          --outName "{params.outname}" \
          --giraffe \
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

rule vg_index:
    input:
        gbz=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz"
    output:
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.dist"
    conda:
        "envs/align_wgs.yaml"
    threads: 8
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=32000,
        cpus=8
    shell:
        r"""
        set -euo pipefail
        vg index -t {threads} --dist-name {output.dist} {input.gbz}
        """

rule vg_minimizer:
    input:
        gbz=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz",
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.dist"
    output:
        min_idx=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.withzip.min",
        zipcodes=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.zipcodes"
    conda:
        "envs/align_wgs.yaml"
    threads: 8
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=32000,
        cpus=8
    shell:
        r"""
        set -euo pipefail
        vg minimizer -t {threads} -d {input.dist} -z {output.zipcodes} -o {output.min_idx} {input.gbz}
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
        dist=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.dist",
        min_idx=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.withzip.min",
        zipcodes=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.vg.zipcodes"
    output:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
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
        export OPENBLAS_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        vg giraffe -Z {input.gbz} --dist-name {input.dist} \
          -m {input.min_idx} \
          -t {threads} -f {input.fq1} -f {input.fq2} -p \
          --rescue-attempts 0 --sample {wildcards.sample} > {output.gam}
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
        # mc_graph has no BAM (vg call works directly on the GAM). Its metrics
        # are emitted by align_metrics_per_gam from vg stats instead.
        ref="|".join(_BCFTOOLS_REFS) if _BCFTOOLS_REFS else "x^"
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
    'mc_graph' (the canonical graph-reference label used downstream). With
    surject removed, mc_graph mapping rate is the true graph mapping rate
    from giraffe, computed straight from vg stats — no surject loss to
    account for."""
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    output:
        cactus=METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv",
        mc_graph=METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv"
    conda:
        "envs/align_metrics.yaml"
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
        expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES),
        expand(METRICS_OUTDIR / "{sample}/{sample}.mc_graph.metrics.tsv", sample=SHORT_SAMPLES)
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


rule bcftools_joint_call:
    input:
        bams=_bams_for_ref,
        bais=_bais_for_ref,
        fasta=lambda wildcards: _ref_fasta(wildcards),
        fai=lambda wildcards: _ref_fai(wildcards),
    output:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    threads:
        BCFTOOLS_THREADS
    wildcard_constraints:
        # mc_graph 'reference' draws its variants from vg_merge_samples (graph-
        # native vg call output), not from a surjected BAM, so this rule must
        # not match it. Without this constraint Snakemake sees two rules
        # producing the same output path.
        ref="|".join(_BCFTOOLS_REFS) if _BCFTOOLS_REFS else "x^"
    resources:
        slurm_partition="long",
        runtime=2880,
        mem_mb=8000,
        cpus=BCFTOOLS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/bcftools/{wildcards.ref}/combined
        bam_list=$(mktemp)
        trap "rm -f $bam_list" EXIT
        printf '%s\n' {input.bams} > "$bam_list"
        bcftools mpileup \
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
          --write-index=tbi \
          -o {output.vcf}
        """

# vg call chunks: contigs are bin-packed by length into bundles so the DAG
# stays tractable (one job per contig × sample blew the SLURM driver's
# RLIMIT_NPROC on the bee genome's ~177-contig fai). Each bundle's outputs
# hold every reference path in that bundle in a single .pg / .gam / .vcf, and
# vg call iterates them in one invocation.
VG_CHUNK_DIR = CACTUS_OUTDIR / "chunks"
VG_BUNDLE_TARGET_BP = 20_000_000  # ceiling per bundle; tune to trade off DAG size vs per-job runtime


def _bin_pack_bundles(fai_path, target_bp):
    """Single-pass greedy bin-pack. Walk contigs in descending length,
    accumulating into the current bundle until its total exceeds target_bp,
    then start a new one. Contigs longer than target_bp end up alone (each
    closes its bundle on the first append)."""
    fai = Path(fai_path)
    if not fai.exists():
        return {}
    pairs = []
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        pairs.append((cols[0], int(cols[1])))
    pairs.sort(key=lambda x: -x[1])
    bundles, current, current_bp, idx = {}, [], 0, 0
    for contig, length in pairs:
        current.append(contig)
        current_bp += length
        if current_bp >= target_bp:
            bundles[f"bundle{idx:03d}"] = current
            idx += 1
            current, current_bp = [], 0
    if current:
        bundles[f"bundle{idx:03d}"] = current
    return bundles


VG_BUNDLES = _bin_pack_bundles(_conspec_fai, VG_BUNDLE_TARGET_BP)
VG_BUNDLE_NAMES = list(VG_BUNDLES.keys())


def _gbz_paths_for_bundle(gbz, bundle_name, outpath):
    """Resolve PanSN reference paths in the GBZ for every contig in this
    bundle, one per line. vg chunk's -P consumes the multi-line file and
    emits one chunk per path in a single invocation."""
    cmds = [f": > {outpath}"]
    for contig in VG_BUNDLES[bundle_name]:
        cmds.append(
            f"vg paths -x {gbz} -L "
            f"| awk -v c={contig} '$0 ~ \"#\" c \"$\" || $0 == c' "
            f"| head -n1 >> {outpath}"
        )
    return " && ".join(cmds)


rule vg_chunk_graph_bundle:
    """Extract one PackedGraph subgraph per bundle. The .pg holds every
    reference path assigned to the bundle plus 2000 nodes of context (enough
    to capture attached snarls), combined into one graph via vg combine."""
    input:
        gbz=VG_GBZ
    output:
        graph=VG_CHUNK_DIR / "vgchunk_{bundle}.pg"
    conda:
        "envs/vg_call.yaml"
    threads: 4
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=32000,
        cpus=4
    params:
        outdir=VG_CHUNK_DIR,
        path_picker=lambda wc, input, output: _gbz_paths_for_bundle(
            input.gbz, wc.bundle, f"{output.graph}.paths.txt"
        ),
    wildcard_constraints:
        bundle="|".join(VG_BUNDLE_NAMES) if VG_BUNDLE_NAMES else "x^"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        export OPENBLAS_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}

        {params.path_picker}
        if [ ! -s {output.graph}.paths.txt ]; then
            echo "ERROR: no GBZ paths matched bundle {wildcards.bundle}" >&2
            exit 1
        fi

        tmpdir=$(mktemp -d -p {params.outdir} chunk_{wildcards.bundle}_XXXX)
        trap "rm -rf $tmpdir" EXIT

        # -c 2000: 2000 nodes of context around each path. vg chunk requires
        # context (or -S snarls) when chunking on paths; this size pulls in
        # the attached variant bubbles without bloating the graph.
        vg chunk -x {input.gbz} -O pg -t {threads} -c 2000 \
            -P {output.graph}.paths.txt \
            -b "$tmpdir/chunk"

        produced=( "$tmpdir"/chunk_*.pg )
        if [ ${{#produced[@]}} -eq 1 ]; then
            mv "${{produced[0]}}" {output.graph}
        else
            vg combine "${{produced[@]}}" > {output.graph}
        fi
        rm -f {output.graph}.paths.txt
        """


rule vg_chunk_gam_bundle:
    """Per-sample GAM chunk for one bundle. Uses the same multi-path file as
    the graph chunk so GAM and graph cover identical reference regions."""
    input:
        gbz=VG_GBZ,
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    output:
        gam=ALIGN_OUTDIR / "{sample}/vg_chunks/{sample}_{bundle}.gam"
    conda:
        "envs/vg_call.yaml"
    threads: 2
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=16000,
        cpus=2
    params:
        outdir=lambda wc: ALIGN_OUTDIR / wc.sample / "vg_chunks",
        path_picker=lambda wc, input, output: _gbz_paths_for_bundle(
            input.gbz, wc.bundle, f"{output.gam}.paths.txt"
        ),
    wildcard_constraints:
        bundle="|".join(VG_BUNDLE_NAMES) if VG_BUNDLE_NAMES else "x^"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        export OPENBLAS_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}

        {params.path_picker}
        if [ ! -s {output.gam}.paths.txt ]; then
            echo "ERROR: no GBZ paths matched bundle {wildcards.bundle}" >&2
            exit 1
        fi

        tmpdir=$(mktemp -d -p {params.outdir} chunk_{wildcards.bundle}_XXXX)
        trap "rm -rf $tmpdir" EXIT

        vg chunk -x {input.gbz} -t {threads} -c 2000 \
            -P {output.gam}.paths.txt \
            -a {input.gam} -g \
            -b "$tmpdir/chunk"

        # GAM is a length-prefixed protobuf stream; cat-concatenation is a
        # valid GAM stream of the union.
        cat "$tmpdir"/chunk_*.gam > {output.gam}
        rm -f {output.gam}.paths.txt
        """


rule vg_call_chunk:
    """Pack + call on a single (sample, bundle). The bundle .pg contains
    every reference path; vg call iterates them in one invocation."""
    input:
        graph=VG_CHUNK_DIR / "vgchunk_{bundle}.pg",
        gam=ALIGN_OUTDIR / "{sample}/vg_chunks/{sample}_{bundle}.gam"
    output:
        vcf=VC_OUTDIR / "vg/chunks/{sample}/{bundle}.vcf.gz",
        tbi=VC_OUTDIR / "vg/chunks/{sample}/{bundle}.vcf.gz.tbi"
    conda:
        "envs/vg_call.yaml"
    threads: VG_THREADS
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=16000,
        cpus=VG_THREADS
    wildcard_constraints:
        bundle="|".join(VG_BUNDLE_NAMES) if VG_BUNDLE_NAMES else "x^"
    shell:
        r"""
        set -euo pipefail
        export OPENBLAS_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        mkdir -p $(dirname {output.vcf})
        pack_tmp=$(mktemp --suffix=.pack)
        trap "rm -f $pack_tmp" EXIT

        # vg call must run on the SAME graph object passed to vg pack — that's
        # why we pass the .pg bundle to both, not the GBZ.
        vg pack -t {threads} -e -Q 5 \
            -x {input.graph} -g {input.gam} -o "$pack_tmp"
        vg call -t {threads} -k "$pack_tmp" -s {wildcards.sample} \
            {input.graph} \
          | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vg_call:
    """Concatenate per-bundle vg call VCFs into the final per-sample VCF.
    Bundles are bin-packed by size, not genome order, so we sort after
    concat to restore canonical chromosome ordering for downstream rules."""
    input:
        vcfs=lambda wc: [VC_OUTDIR / "vg/chunks" / wc.sample / f"{b}.vcf.gz" for b in VG_BUNDLE_NAMES],
        tbis=lambda wc: [VC_OUTDIR / "vg/chunks" / wc.sample / f"{b}.vcf.gz.tbi" for b in VG_BUNDLE_NAMES]
    output:
        vcf=VC_OUTDIR / "vg/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "vg/{sample}.vcf.gz.tbi"
    conda:
        "envs/vg_call.yaml"
    threads: 2
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=2
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})
        tmp_concat=$(mktemp --suffix=.vcf.gz)
        trap "rm -f $tmp_concat" EXIT
        bcftools concat -a -Oz --threads {threads} -o "$tmp_concat" {input.vcfs}
        bcftools sort -Oz -o {output.vcf} "$tmp_concat"
        tabix -p vcf {output.vcf}
        """


rule vg_merge_samples:
    """Merge per-sample vg call VCFs into a cohort VCF that the popgen stack
    consumes the same way as a bcftools-derived merged VCF. Written to the
    bcftools/{ref}/combined path the {ref}-wildcarded popgen rules expect, so
    the mc_graph 'reference' is now a synthetic entry whose variants come
    entirely from vg call rather than from a surjected BAM + bcftools."""
    input:
        vcfs=expand(VC_OUTDIR / "vg/{sample}.vcf.gz", sample=SHORT_SAMPLES),
        tbis=expand(VC_OUTDIR / "vg/{sample}.vcf.gz.tbi", sample=SHORT_SAMPLES),
        fasta=ALIGN_CONSPEC,
        fai=f"{ALIGN_CONSPEC}.fai",
    output:
        vcf=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=8000,
        cpus=4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})

        # Normalize each per-sample VCF first: vg call snarl decomposition can
        # emit different ALT orderings per sample at multiallelic sites, which
        # breaks bcftools merge. norm -m -any splits multiallelics, -f
        # left-aligns against the conspec FASTA so all samples share REF/ALT
        # representations at each locus.
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        normed=()
        for vcf in {input.vcfs}; do
            base=$(basename "$vcf" .vcf.gz)
            out="$tmpdir/$base.norm.vcf.gz"
            bcftools norm -f {input.fasta} -m -any --threads {threads} \
                -Oz -o "$out" "$vcf"
            bcftools index -t "$out"
            normed+=("$out")
        done

        bcftools merge --threads {threads} -Oz -o {output.vcf} "${{normed[@]}}"
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
