# ============================================================================
# Plotting Rules
# ============================================================================

rule plot_pca:
    input:
        cov=PCANGSD_OUTDIR / "{ref}/pcangsd.cov",
        samples=PCANGSD_OUTDIR / "{ref}/samples.txt",
    output:
        pdf=PLOT_OUTDIR / "pca/{ref}_pca.pdf",
        png=PLOT_OUTDIR / "pca/{ref}_pca.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        sample_pop="|".join(s + "," + SAMPLE_TO_POP.get(s, "unknown") for s in SHORT_SAMPLES),
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_pca.py"


rule plot_selection_scan:
    input:
        outliers=PCANGSD_OUTDIR / "{ref}/fst_outliers.tsv",
    output:
        pdf=PLOT_OUTDIR / "selection/{ref}_selection_scan.pdf",
        png=PLOT_OUTDIR / "selection/{ref}_selection_scan.png",
        top_outliers=PLOT_OUTDIR / "selection/{ref}_top_outliers.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        fdr=0.05,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_selection_scan.py"


rule plot_fst:
    input:
        fst=FST_OUTDIR / "fst/{ref}/{pop1}_vs_{pop2}.weir.fst",
    output:
        pdf=PLOT_OUTDIR / "fst/{ref}_{pop1}_vs_{pop2}_fst.pdf",
        png=PLOT_OUTDIR / "fst/{ref}_{pop1}_vs_{pop2}_fst.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        pop1=lambda wildcards: wildcards.pop1,
        pop2=lambda wildcards: wildcards.pop2,
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_fst.py"


rule plot_fst_heatmap:
    """Pairwise mean-FST heatmap across all population pairs for one reference —
    the population-differentiation summary. Replaces reading N separate per-pair
    manhattans with one matrix showing which populations are most divergent."""
    input:
        fsts=[FST_OUTDIR / f"fst/{{ref}}/{p1}_vs_{p2}.weir.fst"
              for (p1, p2) in POP_PAIR_TUPLES],
    output:
        pdf=PLOT_OUTDIR / "fst/{ref}_fst_heatmap.pdf",
        png=PLOT_OUTDIR / "fst/{ref}_fst_heatmap.png",
        table=PLOT_OUTDIR / "fst/{ref}_fst_pairwise.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        pairs=POP_PAIR_TUPLES,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_fst_heatmap.py"


rule plot_afs:
    input:
        afs=FST_OUTDIR / "afs/{ref}.afs.tsv",
    output:
        pdf=PLOT_OUTDIR / "afs/{ref}_afs.pdf",
        png=PLOT_OUTDIR / "afs/{ref}_afs.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    script:
        "../scripts/plot_afs.py"


rule plot_pi:
    input:
        pi=PI_OUTDIR / "{ref}.windowed.pi",
    output:
        pdf=PLOT_OUTDIR / "pi/{ref}_pi.pdf",
        png=PLOT_OUTDIR / "pi/{ref}_pi.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=300,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_pi.py"


rule plot_tajimas_d:
    input:
        tajima=PI_OUTDIR / "{ref}.Tajima.D",
    output:
        pdf=PLOT_OUTDIR / "tajimas_d/{ref}_tajimas_d.pdf",
        png=PLOT_OUTDIR / "tajimas_d/{ref}_tajimas_d.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_tajimas_d.py"


rule plot_diversity_demography:
    """π vs Tajima's D, one point per population — the demographic-quadrant plot.
    Low π + negative D = bottleneck then expansion; high π + positive D =
    structure/balancing; etc. Combines the per-population π and Tajima's D
    genome-wide means for one reference."""
    input:
        pi=PI_OUTDIR / "{ref}.windowed.pi",
        tajima=PI_OUTDIR / "{ref}.Tajima.D",
    output:
        pdf=PLOT_OUTDIR / "diversity/{ref}_pi_vs_tajd.pdf",
        png=PLOT_OUTDIR / "diversity/{ref}_pi_vs_tajd.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_diversity_demography.py"


rule plot_allelic_balance:
    input:
        expand(
            AB_OUTDIR / "{{ref}}/{sample}.allelic_balance.tsv",
            sample=SHORT_SAMPLES,
        ),
    output:
        pdf=PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.pdf",
        png=PLOT_OUTDIR / "allelic_balance/{ref}_allelic_balance.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        samples=SHORT_SAMPLES,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_allelic_balance.py"


# ROH removed (see note near ROH_CFG in the main Snakefile).
_PLOT_ROH_DISABLED = r'''
rule plot_roh:
    input:
        summary=ROH_OUTDIR / "f_roh_summary.tsv",
    output:
        pdf=PLOT_OUTDIR / "roh/{ref}_f_roh.pdf",
        png=PLOT_OUTDIR / "roh/{ref}_f_roh.png",
    conda:
        "../envs/plotting.yaml"
    params:
        ref=lambda wildcards: wildcards.ref,
        sample_pop="|".join(s + "," + SAMPLE_TO_POP.get(s, "unknown") for s in SHORT_SAMPLES),
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    script:
        "../scripts/plot_roh.py"
'''


rule plot_coverage_qc:
    """Per-sample mean depth across all linear refs, with a min-depth threshold
    line. Samples below threshold are highlighted so degraded libraries
    surface in QC at a glance."""
    input:
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv",
    output:
        pdf=PLOT_OUTDIR / "qc/coverage_qc.pdf",
        png=PLOT_OUTDIR / "qc/coverage_qc.png",
    conda:
        "../envs/plotting.yaml"
    params:
        min_depth=QC_MIN_DEPTH,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_coverage_qc.py"


rule plot_alignment_rates:
    input:
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv",
    output:
        pdf=PLOT_OUTDIR / "alignment/alignment_rates.pdf",
        png=PLOT_OUTDIR / "alignment/alignment_rates.png",
    conda:
        "../envs/plotting.yaml"
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_alignment_rates.py"


rule plot_assembly_alignment:
    input:
        metrics=METRICS_OUTDIR / "short_to_assembly_metrics.tsv",
    output:
        pdf=PLOT_OUTDIR / "alignment/assembly_alignment.pdf",
        png=PLOT_OUTDIR / "alignment/assembly_alignment.png",
    conda:
        "../envs/plotting.yaml"
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/plot_assembly_alignment.py"
