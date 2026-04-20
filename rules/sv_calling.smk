# ============================================================================
# SV Calling Rules
# ============================================================================
# Variables (SV_OUTDIR, LONG_SAMPLES, ALIGN_CONSPEC, etc.) are defined in
# the main Snakefile and are available via Snakemake's include mechanism.


# ---------------------------------------------------------------------------
# Per-sample: Align ONT reads to reference with minimap2
# ---------------------------------------------------------------------------
rule sv_align_reads:
    input:
        reads=lambda wc: LONG_READS[wc.sample],
        reference=ALIGN_CONSPEC,
    output:
        bam=SV_OUTDIR / "read_alignments/{sample}.sorted.bam",
        bai=SV_OUTDIR / "read_alignments/{sample}.sorted.bam.bai",
    conda:
        "../envs/sv_calling.yaml"
    threads: SV_THREADS
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=16000,
        cpus=SV_THREADS,
    params:
        preset=lambda wc: "map-ont" if PLATFORM_MAP.get(wc.sample, "") == "ONT" else "map-hifi",
    shell:
        r"""
        set -euo pipefail
        minimap2 -ax {params.preset} -t {threads} --MD {input.reference} {input.reads} \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """


# ---------------------------------------------------------------------------
# Per-sample: Call SVs with Sniffles2
# ---------------------------------------------------------------------------
rule sv_call_sniffles2:
    input:
        bam=SV_OUTDIR / "read_alignments/{sample}.sorted.bam",
        bai=SV_OUTDIR / "read_alignments/{sample}.sorted.bam.bai",
        reference=ALIGN_CONSPEC,
    output:
        vcf=SV_OUTDIR / "read_sv_calls/{sample}.sniffles2.vcf",
    conda:
        "../envs/sv_calling.yaml"
    threads: SV_THREADS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=SV_THREADS,
    params:
        min_sv_size=SV_MIN_SIZE,
        min_support=SV_MIN_SUPPORT,
    shell:
        r"""
        set -euo pipefail
        sniffles \
            --input {input.bam} \
            --vcf {output.vcf} \
            --reference {input.reference} \
            --threads {threads} \
            --minsvlen {params.min_sv_size} \
            --minsupport {params.min_support} \
            --allow-overwrite
        """


# ---------------------------------------------------------------------------
# Per-sample: Call SVs with cuteSV
# ---------------------------------------------------------------------------
rule sv_call_cutesv:
    input:
        bam=SV_OUTDIR / "read_alignments/{sample}.sorted.bam",
        bai=SV_OUTDIR / "read_alignments/{sample}.sorted.bam.bai",
        reference=ALIGN_CONSPEC,
    output:
        vcf=SV_OUTDIR / "read_sv_calls/{sample}.cutesv.vcf",
    conda:
        "../envs/sv_calling.yaml"
    threads: SV_THREADS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=SV_THREADS,
    params:
        min_sv_size=SV_MIN_SIZE,
        min_support=SV_MIN_SUPPORT,
        tmpdir=lambda wc: str(SV_OUTDIR / f"read_sv_calls/{wc.sample}_cutesv_temp"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.tmpdir}
        cuteSV \
            {input.bam} \
            {input.reference} \
            {output.vcf} \
            {params.tmpdir} \
            --threads {threads} \
            --min_size {params.min_sv_size} \
            --min_support {params.min_support} \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --genotype
        rm -rf {params.tmpdir}
        """


# ---------------------------------------------------------------------------
# Per-sample: Align assembly to reference with minimap2 (asm5 preset)
# ---------------------------------------------------------------------------
rule sv_align_assembly:
    input:
        assembly=str(SV_ASSEMBLY_DIR) + "/{sample}/assembly.fasta",
        reference=ALIGN_CONSPEC,
    output:
        paf=SV_OUTDIR / "assembly_alignments/{sample}.sorted.paf",
    conda:
        "../envs/sv_calling.yaml"
    threads: SV_THREADS
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=16000,
        cpus=SV_THREADS,
    shell:
        r"""
        set -euo pipefail
        minimap2 -cx asm5 -t {threads} --cs {input.reference} {input.assembly} \
            | sort -k6,6 -k8,8n > {output.paf}
        """


# ---------------------------------------------------------------------------
# Per-sample: Call SVs with paftools.js
# ---------------------------------------------------------------------------
rule sv_call_paftools:
    input:
        paf=SV_OUTDIR / "assembly_alignments/{sample}.sorted.paf",
        reference=ALIGN_CONSPEC,
    output:
        vcf=SV_OUTDIR / "assembly_sv_calls/{sample}.paftools.vcf",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        min_sv_size=SV_MIN_SIZE,
    shell:
        r"""
        set -euo pipefail
        paftools.js call \
            -f {input.reference} \
            -L {params.min_sv_size} \
            {input.paf} > {output.vcf}
        """


# ---------------------------------------------------------------------------
# Per-sample: Merge read-based callers (Sniffles2 + cuteSV) with SURVIVOR
# ---------------------------------------------------------------------------
rule sv_merge_read_calls:
    input:
        sniffles=SV_OUTDIR / "read_sv_calls/{sample}.sniffles2.vcf",
        cutesv=SV_OUTDIR / "read_sv_calls/{sample}.cutesv.vcf",
    output:
        merged=SV_OUTDIR / "merged_per_sample/{sample}.read_based.merged.vcf",
        vcf_list=temp(SV_OUTDIR / "merged_per_sample/{sample}_read_vcfs.txt"),
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1,
    params:
        breakpoint_slop=SV_BREAKPOINT_SLOP,
        min_sv_size=SV_MIN_SIZE,
        survivor=SV_SURVIVOR,
    shell:
        r"""
        set -euo pipefail
        printf '%s\n' {input.sniffles} {input.cutesv} > {output.vcf_list}
        {params.survivor} merge {output.vcf_list} \
            {params.breakpoint_slop} 1 1 1 0 {params.min_sv_size} \
            {output.merged}
        """


# ---------------------------------------------------------------------------
# Per-sample: Merge assembly-based callers (paftools) with SURVIVOR
# ---------------------------------------------------------------------------
rule sv_merge_asm_calls:
    input:
        paftools=SV_OUTDIR / "assembly_sv_calls/{sample}.paftools.vcf",
    output:
        merged=SV_OUTDIR / "merged_per_sample/{sample}.assembly_based.merged.vcf",
        vcf_list=temp(SV_OUTDIR / "merged_per_sample/{sample}_asm_vcfs.txt"),
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        breakpoint_slop=SV_BREAKPOINT_SLOP,
        min_sv_size=SV_MIN_SIZE,
        survivor=SV_SURVIVOR,
    shell:
        r"""
        set -euo pipefail
        printf '%s\n' {input.paftools} > {output.vcf_list}
        {params.survivor} merge {output.vcf_list} \
            {params.breakpoint_slop} 1 1 1 0 {params.min_sv_size} \
            {output.merged}
        """


# ---------------------------------------------------------------------------
# Per-sample: Merge read-based + assembly-based into high-confidence VCF
# ---------------------------------------------------------------------------
rule sv_merge_per_sample:
    input:
        read_based=SV_OUTDIR / "merged_per_sample/{sample}.read_based.merged.vcf",
        asm_based=SV_OUTDIR / "merged_per_sample/{sample}.assembly_based.merged.vcf",
    output:
        merged=SV_OUTDIR / "merged_per_sample/{sample}.high_confidence.vcf",
        vcf_list=temp(SV_OUTDIR / "merged_per_sample/{sample}_all_vcfs.txt"),
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        breakpoint_slop=SV_BREAKPOINT_SLOP,
        min_sv_size=SV_MIN_SIZE,
        survivor=SV_SURVIVOR,
    shell:
        r"""
        set -euo pipefail
        printf '%s\n' {input.read_based} {input.asm_based} > {output.vcf_list}
        {params.survivor} merge {output.vcf_list} \
            {params.breakpoint_slop} 1 1 1 0 {params.min_sv_size} \
            {output.merged}
        """


# ---------------------------------------------------------------------------
# Pan-sample: Merge all sample high-confidence VCFs with SURVIVOR
# ---------------------------------------------------------------------------
rule sv_merge_survivor:
    input:
        vcfs=expand(SV_OUTDIR / "merged_per_sample/{sample}.high_confidence.vcf", sample=LONG_SAMPLES),
    output:
        vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        vcf_list=SV_OUTDIR / "pan_sample_catalog/all_samples.txt",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    params:
        breakpoint_slop=SV_BREAKPOINT_SLOP,
        min_sv_size=SV_MIN_SIZE,
        survivor=SV_SURVIVOR,
    shell:
        r"""
        set -euo pipefail
        printf '%s\n' {input.vcfs} > {output.vcf_list}
        {params.survivor} merge {output.vcf_list} \
            {params.breakpoint_slop} 1 1 1 1 {params.min_sv_size} \
            {output.vcf}
        """


# ---------------------------------------------------------------------------
# Pan-sample: Merge all sample VCFs with Jasmine (alternative merge strategy)
# ---------------------------------------------------------------------------
rule sv_merge_jasmine:
    input:
        vcfs=expand(SV_OUTDIR / "merged_per_sample/{sample}.high_confidence.vcf", sample=LONG_SAMPLES),
        vcf_list=SV_OUTDIR / "pan_sample_catalog/all_samples.txt",
        genome=ALIGN_CONSPEC,
    output:
        vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.jasmine.vcf",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=16000,
        cpus=1,
    params:
        jasmine_slop=SV_JASMINE_SLOP,
    shell:
        r"""
        set -euo pipefail
        jasmine \
            file_list={input.vcf_list} \
            out_file={output.vcf} \
            genome_file={input.genome} \
            max_dist={params.jasmine_slop} \
            --output_genotypes \
            --normalize_type
        """


# ---------------------------------------------------------------------------
# Pan-sample: Generate catalog statistics and SV support matrix
# ---------------------------------------------------------------------------
rule sv_catalog_stats:
    input:
        survivor_vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        jasmine_vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.jasmine.vcf",
        per_sample_vcfs=expand(SV_OUTDIR / "merged_per_sample/{sample}.high_confidence.vcf", sample=LONG_SAMPLES),
    output:
        stats=SV_OUTDIR / "pan_sample_catalog/catalog_stats.txt",
        support=SV_OUTDIR / "pan_sample_catalog/sv_support_matrix.txt",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        survivor=SV_SURVIVOR,
    shell:
        r"""
        set -euo pipefail

        echo "SURVIVOR catalog statistics:" > {output.stats}
        {params.survivor} stats {input.survivor_vcf} >> {output.stats}

        echo "" >> {output.stats}
        echo "Jasmine catalog statistics:" >> {output.stats}
        bcftools stats {input.jasmine_vcf} >> {output.stats}

        echo "" >> {output.stats}
        echo "SV counts per sample (high-confidence callset):" >> {output.stats}
        for vcf in {input.per_sample_vcfs}; do
            sample=$(basename "$vcf" .high_confidence.vcf)
            count=$(grep -v "^#" "$vcf" | wc -l)
            echo "$sample: $count" >> {output.stats}
        done

        grep -oP 'SUPP_VEC=\K[^,;]+' {input.survivor_vcf} \
            | sed 's/./& /g' > {output.support}
        """


# ---------------------------------------------------------------------------
# Per-sample: Novel sequence summary (bp of novel insertions per sample)
# ---------------------------------------------------------------------------
rule sv_novel_sequence_summary:
    input:
        per_sample_vcfs=expand(SV_OUTDIR / "merged_per_sample/{sample}.high_confidence.vcf", sample=LONG_SAMPLES),
    output:
        summary=SV_OUTDIR / "pan_sample_catalog/novel_sequence_summary.tsv",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        samples=LONG_SAMPLES,
        min_sv_size=SV_MIN_SIZE,
    script:
        "../scripts/sv_novel_sequence_summary.py"


# ---------------------------------------------------------------------------
# Pan-sample: SV counts per sample and shared/unique breakdown
# ---------------------------------------------------------------------------
rule sv_sharing_summary:
    input:
        survivor_vcf=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
    output:
        summary=SV_OUTDIR / "pan_sample_catalog/sv_sharing_summary.tsv",
    conda:
        "../envs/sv_calling.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    params:
        samples=LONG_SAMPLES,
    script:
        "../scripts/sv_sharing_summary.py"
