from pathlib import Path
import pandas as pd

configfile: "config/config.yaml"

SCRIPTS = config["scripts"]

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
ASSEMBLY_GENOME_SIZE = ASSEMBLY_CFG.get("genome_size", "225m")
ASSEMBLY_LINEAGE = ASSEMBLY_CFG.get("lineage", "hymenoptera_odb10")
ASSEMBLY_THREADS = ASSEMBLY_CFG.get("threads", 1)

SV_CFG = config["sv_calling"]
SV_OUTDIR = Path(SV_CFG.get("outdir", "results/sv_calls"))
SV_ASSEMBLY_DIR = SV_CFG.get("assembly_dir", str(ASSEMBLY_OUTDIR / "flye"))
SV_THREADS = SV_CFG.get("threads", 8)
SV_SURVIVOR = SV_CFG.get("survivor_exec", "SURVIVOR")
SV_MIN_SIZE = SV_CFG.get("min_sv_size", 50)
SV_MIN_SUPPORT = SV_CFG.get("min_read_support", 3)
SV_BREAKPOINT_SLOP = SV_CFG.get("breakpoint_slop", 1000)
SV_JASMINE_SLOP = SV_CFG.get("jasmine_slop", 500)
SV_FLANK = SV_CFG.get("flank", 200)

CACTUS_CFG = config["cactus"]
CACTUS_IMAGE = CACTUS_CFG["image"]
CACTUS_JOBSTORE = CACTUS_CFG.get("jobstore", "results/cactus_work")
CACTUS_SEQFILE = Path("results") / "cactus_work" / "seqfile.txt"  # Auto-generated, not manual
CACTUS_OUTDIR = Path(CACTUS_CFG.get("outdir", "results/cactus"))
CACTUS_OUTNAME = CACTUS_CFG.get("outname", "cactus_graph")
CACTUS_REFERENCE = CACTUS_CFG.get("reference", LONG_SAMPLES[0] if LONG_SAMPLES else "")
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
ROH_AUTOSOMES = ROH_CFG.get("autosomes", {})
ROH_BCFTOOLS_ARGS = ROH_CFG.get("bcftools_args", "")

shell.executable("bash")

rule all:
    input:
        # Assembly and QC (long reads only)
        expand(ASSEMBLY_OUTDIR / "nanostat/{sample}_nanostat.txt", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "nanoplot/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "flye/{sample}/assembly.fasta", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "quast/{sample}", sample=LONG_SAMPLES),
        expand(ASSEMBLY_OUTDIR / "busco/{sample}", sample=LONG_SAMPLES),
        # SV calling
        SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.jasmine.vcf",
        SV_OUTDIR / "pan_sample_catalog/catalog_stats.txt",
        SV_OUTDIR / "pan_sample_catalog/sv_support_matrix.txt",
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
        ROH_OUTDIR / "f_roh_summary.tsv"

rule assemble_and_qc:
    input:
        fastq=lambda wildcards: LONG_READS[wildcards.sample]
    output:
        nanostat=ASSEMBLY_OUTDIR / "nanostat/{sample}_nanostat.txt",
        nanoplot=directory(ASSEMBLY_OUTDIR / "nanoplot/{sample}"),
        assembly=ASSEMBLY_OUTDIR / "flye/{sample}/assembly.fasta",
        quast=directory(ASSEMBLY_OUTDIR / "quast/{sample}"),
        busco=directory(ASSEMBLY_OUTDIR / "busco/{sample}")
    conda:
        "envs/assembly.yaml"
    threads:
        ASSEMBLY_THREADS
    resources:
        slurm_partition="long",
        runtime=1440,
        mem_mb=32000,
        cpus=ASSEMBLY_THREADS
    params:
        genome_size=ASSEMBLY_GENOME_SIZE,
        lineage=ASSEMBLY_LINEAGE
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ASSEMBLY_OUTDIR}/nanostat {ASSEMBLY_OUTDIR}/nanoplot {ASSEMBLY_OUTDIR}/flye {ASSEMBLY_OUTDIR}/quast {ASSEMBLY_OUTDIR}/busco

        NanoStat \
          --fastq {input.fastq} \
          --name {wildcards.sample}_nanostat.txt \
          --outdir {ASSEMBLY_OUTDIR}/nanostat \
          --threads {threads}

        NanoPlot \
          --fastq {input.fastq} \
          --prefix {wildcards.sample} \
          --outdir {output.nanoplot} \
          --threads {threads}

        flye \
          --nano-hq {input.fastq} \
          --out-dir {ASSEMBLY_OUTDIR}/flye/{wildcards.sample} \
          --genome-size {params.genome_size} \
          --threads {threads}

        quast \
          {output.assembly} \
          -t {threads} \
          -o {output.quast}

        busco \
          -i {output.assembly} \
          -m genome \
          --lineage_dataset {params.lineage} \
          -c {threads} \
          -o {wildcards.sample} \
          --out_path {ASSEMBLY_OUTDIR}/busco \
          -f
        """

rule call_svs:
    input:
        reads=[LONG_READS[sample] for sample in LONG_SAMPLES],
        assemblies=[f"{SV_ASSEMBLY_DIR}/{sample}/assembly.fasta" for sample in LONG_SAMPLES]
    output:
        survivor=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        jasmine=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.jasmine.vcf",
        stats=SV_OUTDIR / "pan_sample_catalog/catalog_stats.txt",
        support=SV_OUTDIR / "pan_sample_catalog/sv_support_matrix.txt",
        augref=SV_OUTDIR / "augref/augmented_reference.fasta"
    conda:
        "envs/sv_calling.yaml"
    threads:
        SV_THREADS
    resources:
        slurm_partition="long",
        runtime=1440,
        mem_mb=32000,
        cpus=SV_THREADS
    params:
        script=SCRIPTS["call_svs"],
        reference=ALIGN_CONSPEC,
        assembly_dir=SV_ASSEMBLY_DIR,
        outdir=SV_OUTDIR,
        samples=",".join(LONG_SAMPLES),
        survivor=SV_SURVIVOR,
        min_size=SV_MIN_SIZE,
        min_support=SV_MIN_SUPPORT,
        breakpoint_slop=SV_BREAKPOINT_SLOP,
        jasmine_slop=SV_JASMINE_SLOP,
        flank=SV_FLANK
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SV_OUTDIR}
        
        REFERENCE="{params.reference}" \
        ASSEMBLY_DIR="{params.assembly_dir}" \
        OUTPUT_DIR="{params.outdir}" \
        SAMPLES="{params.samples}" \
        SURVIVOR_EXEC="{params.survivor}" \
        MIN_SV_SIZE="{params.min_size}" \
        MIN_READ_SUPPORT="{params.min_support}" \
        BREAKPOINT_SLOP="{params.breakpoint_slop}" \
        JASMINE_SLOP="{params.jasmine_slop}" \
        FLANK="{params.flank}" \
        THREADS="{threads}" \
        bash {params.script}
        """

rule generate_cactus_seqfile:
    input:
        reference=ALIGN_CONSPEC,
        assemblies=expand(ASSEMBLY_OUTDIR / "flye/{sample}/assembly.fasta", sample=LONG_SAMPLES)
    output:
        seqfile=CACTUS_SEQFILE
    params:
        ref_name=CACTUS_REFERENCE
    run:
        # Generate seqfile in cactus format with reference as backbone
        # First line: reference_name /path/to/reference.fasta
        # Following lines: sample_name /path/to/assembly.fasta
        with open(output.seqfile, "w") as f:
            # Write reference first (backbone)
            f.write(f"{params.ref_name} {input.reference}\n")
            # Write assemblies
            for sample, assembly in zip(LONG_SAMPLES, input.assemblies):
                f.write(f"{sample} {assembly}\n")

rule make_cactus_graph:
    input:
        seqfile=CACTUS_SEQFILE
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
        reference=CACTUS_REFERENCE,
        ref_contigs=CACTUS_REF_CONTIGS,
        extra_args=CACTUS_EXTRA_ARGS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        cactus-pangenome \
          {params.jobstore} \
          {input.seqfile} \
          --reference "{params.reference}" \
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

rule align_wgs:
    input:
        fq1=lambda wildcards: SHORT_READS_R1[wildcards.sample],
        fq2=lambda wildcards: SHORT_READS_R2[wildcards.sample],
        augref=ALIGN_AUGREF,
        augref_bwt=f"{ALIGN_AUGREF}.bwt",
        conspec=ALIGN_CONSPEC,
        conspec_bwt=f"{ALIGN_CONSPEC}.bwt",
        hetspec=ALIGN_HETSPEC,
        hetspec_bwt=f"{ALIGN_HETSPEC}.bwt",
        cactus=ALIGN_CACTUS_GBZ
    output:
        augref_bam=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam",
        augref_bai=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai",
        conspec_bam=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam",
        conspec_bai=ALIGN_OUTDIR / "{sample}/{sample}.conspec.bam.bai",
        hetspec_bam=ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam",
        hetspec_bai=ALIGN_OUTDIR / "{sample}/{sample}.hetspec.bam.bai",
        cactus_gam=ALIGN_OUTDIR / "{sample}/{sample}.cactus.gam"
    conda:
        "envs/align_wgs.yaml"
    threads:
        ALIGN_THREADS
    resources:
        slurm_partition="long",
        runtime=960,
        mem_mb=16000,
        cpus=ALIGN_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ALIGN_OUTDIR}/{wildcards.sample}

        bwa mem -k 17 -O 5,5 -E 2,2 -Y -t {threads} {input.augref} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.augref_bam}
        samtools index {output.augref_bam}

        bwa mem -k 17 -O 5,5 -E 2,2 -Y -t {threads} {input.conspec} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.conspec_bam}
        samtools index {output.conspec_bam}

        bwa mem -k 17 -O 5,5 -E 2,2 -Y -t {threads} {input.hetspec} {input.fq1} {input.fq2} | \
          samtools sort -@ {threads} -o {output.hetspec_bam}
        samtools index {output.hetspec_bam}

        vg giraffe -Z {input.cactus} -t {threads} -f {input.fq1} -f {input.fq2} -p --rescue-attempts 0 --sample {wildcards.sample} > {output.cactus_gam}
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
        mean_mapq=$(samtools view -F 4 {input.bam} | awk '{{sum+=$5; n++}} END {{if(n>0) printf "%.6f", sum/n; else print "0"}}')
        mean_depth=$(samtools depth -a -q {params.min_mapq} -Q {params.min_baseq} {input.bam} | \
            awk '{{sum+=$3; n++}} END {{if (n>0) printf "%.6f", sum/n; else print "0"}}')
        map_rate=$(awk -v m="$mapped" -v t="$total" 'BEGIN {{if (t>0) printf "%.6f", m/t; else print "0"}}')

        cat > {output} << 'EOF'
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
        aligned_reads=$(echo "$stats" | awk '/Total aligned:/ {{print $3}}')
        mean_mapq=$(echo "$stats" | awk -F'mean ' '/Mapping quality:/ {{print $2}}' | awk '{{print $1}}')

        cat > {output} << 'EOF'
sample	alignment_type	total_reads	aligned_reads	mapping_rate	mean_mapq	mean_depth
{wildcards.sample}	cactus	NA	$aligned_reads	NA	$mean_mapq	NA
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

def _ref_lengths_path(wildcards):
    if wildcards.ref in ROH_GENOME_LENGTHS:
        return ROH_GENOME_LENGTHS[wildcards.ref]
    return str(ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv")

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
        slurm_partition="short",
        runtime=240,
        mem_mb=8000,
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
        vg pack -t {threads} -Z {input.gbz} -g {input.gam} -o "$pack_tmp"
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

rule afs_per_ref:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
    output:
        afs=FST_OUTDIR / "afs/{ref}.afs.tsv"
    conda:
        "envs/fst_afs.yaml"
    params:
        bins=AFS_BINS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    run:
        from pathlib import Path
        import subprocess
        
        out_path = Path(output.afs)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        freq_prefix = Path(output.afs).parent / "freq_tmp"
        subprocess.run(
            f"vcftools --gzvcf {input.vcf} --freq2 --out {freq_prefix}",
            shell=True, check=True, capture_output=True
        )
        
        freq_path = freq_prefix.with_suffix(".frq")
        bins = params.bins
        counts = [0 for _ in range(len(bins) - 1)]
        
        for line in freq_path.read_text().strip().splitlines()[1:]:
            parts = line.split()
            if len(parts) < 6:
                continue
            freqs = []
            for item in parts[5:]:
                if ":" not in item:
                    continue
                try:
                    freqs.append(float(item.split(":")[1]))
                except ValueError:
                    continue
            if not freqs:
                continue
            maf = min(freqs)
            for i in range(len(bins) - 1):
                if bins[i] <= maf < bins[i + 1]:
                    counts[i] += 1
                    break
        
        out_path.write_text("bin_low\tbin_high\tcount\n")
        with out_path.open("a") as handle:
            for i in range(len(counts)):
                handle.write(f"{bins[i]}\t{bins[i+1]}\t{counts[i]}\n")
        
        freq_path.unlink(missing_ok=True)
        freq_prefix.with_suffix(".log").unlink(missing_ok=True)

rule fst_per_ref_pair:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
    output:
        FST_OUTDIR / "fst/{ref}/{pop1}_vs_{pop2}.weir.fst"
    conda:
        "envs/fst_afs.yaml"
    params:
        window_args=FST_WINDOW_ARGS,
        pop1_samples=lambda wildcards: "\n".join(POP_SAMPLES[wildcards.pop1]),
        pop2_samples=lambda wildcards: "\n".join(POP_SAMPLES[wildcards.pop2])
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FST_OUTDIR}/fst/{wildcards.ref}
        prefix={FST_OUTDIR}/fst/{wildcards.ref}/{wildcards.pop1}_vs_{wildcards.pop2}
        
        # Create temporary population files
        echo "{params.pop1_samples}" > /tmp/pop1_samples.txt
        echo "{params.pop2_samples}" > /tmp/pop2_samples.txt
        
        vcftools --gzvcf {input.vcf} \
          --weir-fst-pop /tmp/pop1_samples.txt \
          --weir-fst-pop /tmp/pop2_samples.txt \
          {params.window_args} \
          --out $prefix > /dev/null
        
        rm /tmp/pop1_samples.txt /tmp/pop2_samples.txt
        """

rule pi_per_ref:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
    output:
        PI_OUTDIR / "{ref}.windowed.pi"
    conda:
        "envs/pi.yaml"
    params:
        window_size=PI_WINDOW_SIZE,
        window_step=PI_WINDOW_STEP
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {PI_OUTDIR}
        vcftools --gzvcf {input.vcf} \
          --window-pi {params.window_size} \
          --window-pi-step {params.window_step} \
          --out {PI_OUTDIR}/{wildcards.ref} > /dev/null
        """

rule allelic_balance_per_sample:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi"
    output:
        summary=AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv",
        raw=AB_OUTDIR / "{ref}/{sample}.allelic_balance.raw.tsv"
    conda:
        "envs/allelic_balance.yaml"
    params:
        bins=AB_BINS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    run:
        import statistics
        from pathlib import Path
        import pysam
        
        bins = params.bins
        vcf_path = input.vcf
        out_path = Path(output.summary)
        raw_path = Path(output.raw)
        
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        counts = [0 for _ in range(len(bins) - 1)]
        ratios = []
        
        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        if len(samples) != 1:
            raise ValueError(f"Expected 1 sample in {vcf_path}, found {len(samples)}")
        sample = samples[0]
        
        raw_data = []
        for rec in vcf.fetch():
            if "AD" not in rec.format or "GT" not in rec.format:
                continue
            gt = rec.samples[sample].get("GT")
            if gt is None or len(gt) < 2:
                continue
            if gt[0] == gt[1]:
                continue
            ad = rec.samples[sample].get("AD")
            if ad is None or len(ad) < 2:
                continue
            ref_depth = ad[0]
            alt_depth = ad[1]
            if ref_depth is None or alt_depth is None:
                continue
            total = ref_depth + alt_depth
            if total == 0:
                continue
            ratio = ref_depth / total
            ratios.append(ratio)
            site = f"{rec.contig}:{rec.pos}"
            raw_data.append((site, ref_depth, alt_depth, ratio))
            for i in range(len(bins) - 1):
                if bins[i] <= ratio < bins[i + 1]:
                    counts[i] += 1
                    break
        
        mean_ratio = statistics.mean(ratios) if ratios else 0.0
        median_ratio = statistics.median(ratios) if ratios else 0.0
        total_hets = len(ratios)
        
        with out_path.open("w") as handle:
            handle.write("sample\tref\ttotal_hets\tmean_ref_ratio\tmedian_ref_ratio\n")
            handle.write(f"{sample}\t{wildcards.ref}\t{total_hets}\t{mean_ratio:.6f}\t{median_ratio:.6f}\n")
            handle.write("bin_low\tbin_high\tcount\n")
            for i in range(len(counts)):
                handle.write(f"{bins[i]}\t{bins[i+1]}\t{counts[i]}\n")
        
        with raw_path.open("w") as handle:
            handle.write("site\tref_depth\talt_depth\tref_ratio\n")
            for site, ref_depth, alt_depth, ratio in raw_data:
                handle.write(f"{site}\t{ref_depth}\t{alt_depth}\t{ratio:.6f}\n")

rule allelic_balance_summary:
    input:
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES)
    output:
        summary=AB_OUTDIR / "allelic_balance_summary.tsv"
    conda:
        "envs/allelic_balance.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    run:
        from pathlib import Path
        
        out_path = Path(output.summary)
        rows = ["sample\tref\ttotal_hets\tmean_ref_ratio\tmedian_ref_ratio"]
        
        for path in input:
            text = Path(path).read_text().splitlines()
            if len(text) < 2:
                continue
            rows.append(text[1])
        
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(rows) + "\n")

rule roh_per_sample:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi",
        lengths=lambda wildcards: ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv"
    output:
        roh=ROH_OUTDIR / "{ref}/{sample}.roh.tsv",
        froh=ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv"
    conda:
        "envs/roh.yaml"
    params:
        autosomes=ROH_AUTOSOMES,
        bcftools_args=ROH_BCFTOOLS_ARGS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    run:
        import subprocess
        from pathlib import Path
        import tempfile
        
        out_path = Path(output.roh)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Run bcftools roh
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            roh_tmp = tmp.name
        
        subprocess.run(
            f"bcftools roh {params.bcftools_args} -o {roh_tmp} {input.vcf}",
            shell=True, check=True, capture_output=True
        )
        
        roh_path = Path(roh_tmp)
        lengths_path = Path(input.lengths)
        out_roh = Path(output.roh)
        out_froh = Path(output.froh)
        
        # Read lengths
        lengths = {}
        for line in lengths_path.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            contig, size = parts[0], parts[1]
            try:
                lengths[contig] = int(size)
            except ValueError:
                continue
        
        # Get autosomes for this ref
        autosome_map = params.autosomes
        autosomes = autosome_map.get(wildcards.ref, [])
        if autosomes:
            genome_len = sum(lengths.get(c, 0) for c in autosomes)
        else:
            genome_len = sum(lengths.values())
        
        # Parse ROH results
        roh_rows = []
        roh_total = 0
        for line in roh_path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            chrom = parts[2] if len(parts) >= 5 else parts[1]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            length = end - start + 1
            if autosomes and chrom not in autosomes:
                continue
            roh_total += length
            roh_rows.append((chrom, start, end, length))
        
        # Write outputs
        out_roh.write_text("chrom\tstart\tend\tlength\n")
        with out_roh.open("a") as handle:
            for chrom, start, end, length in roh_rows:
                handle.write(f"{chrom}\t{start}\t{end}\t{length}\n")
        
        f_roh = (roh_total / genome_len) if genome_len > 0 else 0.0
        out_froh.write_text("sample\tref\troh_bp\tgenome_bp\tf_roh\n")
        with out_froh.open("a") as handle:
            handle.write(f"{wildcards.sample}\t{wildcards.ref}\t{roh_total}\t{genome_len}\t{f_roh:.6f}\n")
        
        roh_path.unlink()

rule roh_summary:
    input:
        expand(ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES)
    output:
        summary=ROH_OUTDIR / "f_roh_summary.tsv"
    conda:
        "envs/roh.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    run:
        from pathlib import Path
        
        out_path = Path(output.summary)
        rows = ["sample\tref\troh_bp\tgenome_bp\tf_roh"]
        
        for path in input:
            text = Path(path).read_text().splitlines()
            if len(text) < 2:
                continue
            rows.append(text[1])
        
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(rows) + "\n")

rule ref_lengths:
    input:
        fasta=lambda wildcards: _gatk_ref_fasta(wildcards),
        fai=lambda wildcards: _gatk_ref_fai(wildcards)
    output:
        lengths=ROH_OUTDIR / "lengths" / "{ref}.lengths.tsv"
    conda:
        "envs/roh.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=2000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_OUTDIR}/lengths
        awk -F '\t' 'BEGIN {{OFS="\t"}} {{print $1,$2}}' {input.fai} > {output.lengths}
        """
