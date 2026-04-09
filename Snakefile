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

# mc_graph is the Minigraph-Cactus pangenome surjected into conspec linear
# coordinates via vg surject. It shares the conspec FASTA and all its indices
# for downstream variant calling, so all {ref}-wildcarded rules (GATK, merge,
# FST, AFS, π, allelic balance, ROH) automatically pick it up from here.
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

METRICS_CFG = config["align_metrics"]
METRICS_OUTDIR = Path(METRICS_CFG.get("outdir", "results/align_metrics"))
METRICS_REF_TYPES = REFERENCE_NAMES
METRICS_GAM = METRICS_CFG.get("include_gam", True)
METRICS_MIN_MAPQ = METRICS_CFG.get("min_mapq", 0)
METRICS_MIN_BASEQ = METRICS_CFG.get("min_baseq", 0)
METRICS_THREADS = METRICS_CFG.get("threads", 2)

VC_CFG = config["variant_calling"]
GATK_CFG = VC_CFG["gatk"]
VG_CFG = VC_CFG["vg"]
VC_OUTDIR = Path(VC_CFG.get("outdir", "results/variants"))
GATK_THREADS = GATK_CFG.get("threads", 4)
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
wildcard_constraints:
    ref="|".join(REFERENCE_NAMES)

rule all:
    input:
        # Assembly and QC (long reads only)
        expand(ASSEMBLY_OUTDIR / "nanostat/{sample}_nanostat.txt", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "nanoplot/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "hifiasm/{sample}/assembly.fasta", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "quast/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "busco/{sample}", sample=LONG_SAMPLES),
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
        expand(ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam", sample=SHORT_SAMPLES),
        expand(ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam.bai", sample=SHORT_SAMPLES),
        # Alignment metrics
        expand(METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv", sample=SHORT_SAMPLES, ref=METRICS_REF_TYPES),
        expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES),
        METRICS_OUTDIR / "alignment_metrics.tsv",
        # Variant calling
        expand(VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "vg/{sample}.vcf.gz", sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "vg/{sample}.vcf.gz.tbi", sample=SHORT_SAMPLES),
        expand(VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz", ref=REFERENCE_NAMES),
        expand(VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi", ref=REFERENCE_NAMES),
        # Population genetics
        expand(FST_OUTDIR / "afs/{ref}.afs.tsv", ref=REFERENCE_NAMES),
        [FST_OUTDIR / f"fst/{ref}/{pop1}_vs_{pop2}.weir.fst" for ref in REFERENCE_NAMES for (pop1, pop2) in POP_PAIR_TUPLES],
        expand(PI_OUTDIR / "{ref}.windowed.pi", ref=REFERENCE_NAMES),
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
        expand(PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.png", ref=REFERENCE_NAMES),
        expand(PLOT_OUTDIR / "roh/{ref}_f_roh.png", ref=REFERENCE_NAMES)

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
    shell:
        r"""
        set -euo pipefail
        python - << 'PYEOF'
import pysam

vcf_path   = "{input.vcf}"
ref_path   = "{input.reference}"
out_fasta  = "{output.fasta}"
flank_size = {params.flank}
min_ins    = {params.min_ins}
samples    = "{params.samples}".split()

vcf       = pysam.VariantFile(vcf_path)
reference = pysam.FastaFile(ref_path)

count = 0
seen_headers = {{}}
with open(out_fasta, "w") as out:
    for record in vcf:
        alt_seq    = record.alts[0]
        ref_allele = record.ref
        if len(alt_seq) <= len(ref_allele) or alt_seq.startswith("<"):
            continue
        ins_len = len(alt_seq) - len(ref_allele)
        if ins_len < min_ins:
            continue
        novel_part = alt_seq[len(ref_allele):]
        chrom = record.chrom
        pos   = record.pos
        try:
            left_flank  = reference.fetch(chrom, max(0, pos - flank_size), pos)
            right_flank = reference.fetch(chrom, pos, pos + flank_size)
        except KeyError:
            continue
        supp_vec    = record.info.get("SUPP_VEC", "")
        sample_name = "unknown"
        for i, val in enumerate(supp_vec):
            if val == "1" and i < len(samples):
                sample_name = samples[i]
                break
        full_seq  = left_flank + novel_part + right_flank
        base_name = "SV_" + chrom + "_" + str(pos) + "_" + sample_name + "_ins" + str(ins_len)
        if base_name in seen_headers:
            seen_headers[base_name] += 1
            header = ">" + base_name + "_dup" + str(seen_headers[base_name])
        else:
            seen_headers[base_name] = 0
            header = ">" + base_name
        out.write(header + "\n" + full_seq + "\n")
        count += 1

print("Extracted", count, "flanked sequences to", out_fasta)
PYEOF
        """

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

rule align_bwa_augref:
    input:
        fq1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        fq2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
        ref=ALIGN_AUGREF,
        bwt=f"{ALIGN_AUGREF}.bwt"
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai"
    conda:
        "envs/align_wgs.yaml"
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
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
        fq1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        fq2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
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
        runtime=480,
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
        fq1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        fq2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
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
        runtime=480,
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
        fq1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        fq2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
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
        runtime=480,
        mem_mb=32000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}
        vg giraffe -Z {input.gbz} --dist-name {input.dist} \
          -m {input.min_idx} \
          -t {threads} -f {input.fq1} -f {input.fq2} -p \
          --rescue-attempts 0 --sample {wildcards.sample} > {output.gam}
        """

rule vg_surject:
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam",
        gbz=ALIGN_CACTUS_GBZ,
    output:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.mc_graph.bam.bai",
    conda:
        "envs/align_wgs.yaml"
    threads: ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=16000,
        cpus=ALIGN_THREADS,
    shell:
        r"""
        set -euo pipefail
        vg surject \
            -x {input.gbz} \
            -b \
            -t {threads} \
            {input.gam} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools view -H {output.bam} \
            | sed 's/SN:reference#[0-9]*#/SN:/g' \
            > {output.bam}.header
        samtools reheader {output.bam}.header {output.bam} > {output.bam}.tmp
        mv {output.bam}.tmp {output.bam}
        rm {output.bam}.header
        samtools index {output.bam}
        """

rule align_metrics_per_bam:
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam.bai"
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
    shell:
        r"""
        set -euo pipefail
        mkdir -p {METRICS_OUTDIR}/{wildcards.sample}

        total=$(samtools view -c -F 0x900 {input.bam})
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
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    output:
        METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv"
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
        mkdir -p {METRICS_OUTDIR}/{wildcards.sample}

        stats=$(vg stats -a {input.gam})
        total_reads=$(echo "$stats" | awk '/Total alignments:/ {{print $3}}')
        aligned_reads=$(echo "$stats" | awk '/Total aligned:/ {{print $3}}')
        mean_mapq=$(echo "$stats" | awk '/Mapping quality:/ {{print $NF}}')

        cat > {output} << EOF
sample	alignment_type	total_reads	aligned_reads	mapping_rate	mean_mapq	mean_depth
{wildcards.sample}	cactus	$total_reads	$aligned_reads	NA	$mean_mapq	NA
EOF
        """

rule align_metrics_summary:
    input:
        expand(METRICS_OUTDIR / "{sample}/{sample}.{ref}.metrics.tsv", sample=SHORT_SAMPLES, ref=METRICS_REF_TYPES),
        expand(METRICS_OUTDIR / "{sample}/{sample}.cactus.metrics.tsv", sample=SHORT_SAMPLES)
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

def _gatk_ref_fasta(wildcards):
    if wildcards.ref not in REFS_NESTED:
        raise ValueError(f"Unknown reference: {wildcards.ref}")
    return REFS_NESTED[wildcards.ref]["fasta"]

def _gatk_ref_fai(wildcards):
    return f"{_gatk_ref_fasta(wildcards)}.fai"

def _gatk_ref_dict(wildcards):
    return str(Path(_gatk_ref_fasta(wildcards)).with_suffix(".dict"))

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
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
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
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
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
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
        bwa index {input.fasta}
        """

rule gatk_haplotypecaller:
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.{ref}.bam.bai",
        fasta=lambda wildcards: _gatk_ref_fasta(wildcards),
        fai=lambda wildcards: _gatk_ref_fai(wildcards),
        dict=lambda wildcards: _gatk_ref_dict(wildcards)
    output:
        vcf=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi"
    conda:
        "envs/gatk.yaml"
    threads:
        GATK_THREADS
    resources:
        slurm_partition="long",
        runtime=2880,
        mem_mb=16000,
        cpus=GATK_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/gatk/{wildcards.ref}/merged
        gatk HaplotypeCaller \
          -R {input.fasta} \
          -I {input.bam} \
          -O {output.vcf} \
          --native-pair-hmm-threads {threads}
        """

rule vg_call:
    input:
        gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam",
        gbz=VG_GBZ
    output:
        vcf=VC_OUTDIR / "vg/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "vg/{sample}.vcf.gz.tbi"
    conda:
        "envs/vg_call.yaml"
    threads:
        VG_THREADS
    resources:
        slurm_partition="long",
        runtime=480,
        mem_mb=16000,
        cpus=VG_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/vg
        pack_tmp=$(mktemp --suffix=.pack)
        vg pack -t {threads} -x {input.gbz} -g {input.gam} -o "$pack_tmp"
        vg call -t {threads} -k "$pack_tmp" -s {wildcards.sample} {input.gbz} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        rm -f "$pack_tmp"
        """

def _gatk_vcfs_for_ref(wildcards):
    return [f"{VC_OUTDIR}/gatk/{wildcards.ref}/merged/{sample}.vcf.gz" for sample in SHORT_SAMPLES]

rule merge_gatk_vcfs:
    input:
        vcfs=_gatk_vcfs_for_ref,
        tbis=lambda wildcards: [f"{VC_OUTDIR}/gatk/{wildcards.ref}/merged/{sample}.vcf.gz.tbi" for sample in SHORT_SAMPLES]
    output:
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
    conda:
        "envs/fst_afs.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {VC_OUTDIR}/gatk/{wildcards.ref}/combined
        bcftools merge -m all -Oz -o {output.vcf} {input.vcfs}
        tabix -p vcf {output.vcf}
        """

include: "rules/popgen.smk"
include: "rules/plotting.smk"
