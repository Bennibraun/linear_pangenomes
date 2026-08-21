# ============================================================================
# Population Genetics Rules
# ============================================================================
# Variables (VC_OUTDIR, FST_OUTDIR, PI_OUTDIR, AB_OUTDIR, ROH_OUTDIR, etc.)
# and helper functions (_ref_fasta, _ref_fai) are defined in the
# main Snakefile and are available via Snakemake's include mechanism.


# ROH removed (see note near ROH_CFG in the main Snakefile).
# def _ref_lengths_path(wildcards):
#     if wildcards.ref in ROH_GENOME_LENGTHS:
#         return ROH_GENOME_LENGTHS[wildcards.ref]
#     return str(ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv")


# ---------------------------------------------------------------------------
# Allele Frequency Spectrum (VCF-based for all refs)
# ---------------------------------------------------------------------------
rule afs_per_ref:
    """VCF-based folded SFS from bcftools merged VCF. Uses vcftools --freq2
    to compute per-site allele frequencies, then bins into the canonical
    AFS histogram. Uniform method across all refs (including mc_graph)
    ensures the spectra are directly comparable."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        afs=FST_OUTDIR / "afs/{ref}.afs.tsv",
    conda:
        "../envs/fst_afs.yaml"
    params:
        bins=AFS_BINS,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=24000,
        cpus=1,
    script:
        "../scripts/afs_per_ref.py"


rule afs_by_region_augref:
    """Core-vs-accessory (SV_*) folded SFS for augref only. Distinguishes the
    two Tajima's-D hypotheses directly: a real excess of shared ancestral
    variation in the accessory sequence should show up as an SFS shifted
    toward common variants relative to core; a calling-depth artifact would
    instead show a rare-variant deficit tracking with the depth/MAPQ gap in
    region_qc_summary.tsv rather than a clean shift in the whole spectrum."""
    input:
        vcf=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz.tbi",
    output:
        afs_by_region=FST_OUTDIR / "afs/augref.afs_by_region.tsv",
    conda:
        "../envs/fst_afs.yaml"
    params:
        bins=AFS_BINS,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=24000,
        cpus=1,
    script:
        "../scripts/afs_per_ref.py"


# ---------------------------------------------------------------------------
# FST between population pairs
# ---------------------------------------------------------------------------
rule fst_per_ref_pair:
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        FST_OUTDIR / "fst/{ref}/{pop1}_vs_{pop2}.weir.fst"
    conda:
        "../envs/fst_afs.yaml"
    params:
        window_args=FST_WINDOW_ARGS,
        pop1_samples=lambda wildcards: " ".join(POP_SAMPLES[wildcards.pop1]),
        pop2_samples=lambda wildcards: " ".join(POP_SAMPLES[wildcards.pop2])
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FST_OUTDIR}/fst/{wildcards.ref}
        prefix={FST_OUTDIR}/fst/{wildcards.ref}/{wildcards.pop1}_vs_{wildcards.pop2}

        pop1_tmp=$(mktemp -p {FST_OUTDIR}/fst/{wildcards.ref})
        pop2_tmp=$(mktemp -p {FST_OUTDIR}/fst/{wildcards.ref})
        trap "rm -f $pop1_tmp $pop2_tmp" EXIT

        for s in {params.pop1_samples}; do echo "$s"; done > "$pop1_tmp"
        for s in {params.pop2_samples}; do echo "$s"; done > "$pop2_tmp"

        vcftools --gzvcf {input.vcf} \
          --weir-fst-pop "$pop1_tmp" \
          --weir-fst-pop "$pop2_tmp" \
          {params.window_args} \
          --out $prefix > /dev/null
        """


# ---------------------------------------------------------------------------
# F1 concordance between reference SNP call sets (Jeon et al. 2026 Table S2)
# ---------------------------------------------------------------------------
rule f1_concordance_pair:
    """bcftools isec between a pair of refs that share coordinates, then
    F1 from (shared, a_only, b_only). Restricted to contigs present in both
    VCFs so refs with extra contigs (e.g. augref's SV_*) still work."""
    input:
        vcf_a=VC_OUTDIR / "bcftools/{ref_a}/combined/merged.vcf.gz",
        vcf_b=VC_OUTDIR / "bcftools/{ref_b}/combined/merged.vcf.gz",
        tbi_a=VC_OUTDIR / "bcftools/{ref_a}/combined/merged.vcf.gz.tbi",
        tbi_b=VC_OUTDIR / "bcftools/{ref_b}/combined/merged.vcf.gz.tbi",
    output:
        f1=F1_OUTDIR / "{ref_a}__vs__{ref_b}.f1.tsv",
    conda:
        "../envs/plotting.yaml"
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=8000,
        cpus=1,
    params:
        ref_a=lambda wc: wc.ref_a,
        ref_b=lambda wc: wc.ref_b,
    script:
        "../scripts/f1_concordance.py"


rule f1_concordance_summary:
    """Aggregate all pairwise F1 results into a single table."""
    input:
        [F1_OUTDIR / f"{a}__vs__{b}.f1.tsv" for a, b in F1_PAIRS],
    output:
        summary=F1_OUTDIR / "f1_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    run:
        from pathlib import Path

        out_path = Path(output.summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rows, header = [], None
        for path in input:
            lines = Path(path).read_text().strip().splitlines()
            if not lines:
                continue
            if header is None:
                header = lines[0]
            if len(lines) > 1:
                rows.append(lines[1])
        if header is None:
            out_path.write_text("ref_a\tref_b\tshared\ta_only\tb_only\tprecision\trecall\tf1\n")
        else:
            out_path.write_text("\n".join([header] + rows) + "\n")


# ---------------------------------------------------------------------------
# Per-reference summary table (Jeon et al. 2026 Table 1 equivalent)
# ---------------------------------------------------------------------------
rule graph_total_length:
    """Total sequence length of the pangenome graph (sum of node lengths) via
    `vg stats -l`. Runs in VG_IMAGE so the summary rule's env stays light. This
    is the graph's true size, distinct from mc_graph's SNP-calling coordinate
    space (which is conspec after surjection)."""
    input:
        gbz=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.gbz",
    output:
        length=CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.total_length.txt",
    container:
        VG_IMAGE
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        # `vg stats -l` prints "length <N>"; keep the number.
        vg stats -l {input.gbz} | awk '{{print $NF}}' > {output.length}
        """


# ---------------------------------------------------------------------------
rule per_ref_summary:
    """One row per reference with the headline comparison metrics: SNP counts
    (raw + filtered), SV count, mean mapping rate/quality, mean π, mean Tajima's
    D, mean FST, genome length, size vs conspec, and (graph only) total graph
    sequence. SV counts come from the ref-appropriate source: augref from the
    long-read SV catalog (sniffles/cuteSV -> SURVIVOR), mc_graph from its own
    vg-call VCFs; these are independent discovery paths."""
    input:
        vcfs=expand(
            VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
            ref=REFERENCE_NAMES,
        ),
        # mc_graph SVs: vg-call VCFs (graph only; empty when graph off).
        sv_vcfs=expand(VC_OUTDIR / "sv/vg/{sample}.vcf.gz", sample=SHORT_SAMPLES) if INCLUDE_GRAPH else [],
        # augref SVs: the long-read pan-sample SV catalog that builds augref.
        sv_catalog=SV_OUTDIR / "pan_sample_catalog/pan_sample_catalog.survivor.vcf",
        # Total graph sequence (precomputed via `vg stats -l`; graph only).
        graph_length=[CACTUS_OUTDIR / f"{CACTUS_OUTNAME}.total_length.txt"] if INCLUDE_GRAPH else [],
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv",
        pi=expand(PI_OUTDIR / "{ref}.windowed.pi", ref=REFERENCE_NAMES),
        tajima=expand(PI_OUTDIR / "{ref}.Tajima.D", ref=REFERENCE_NAMES),
        fst=[FST_OUTDIR / f"fst/{ref}/{p1}_vs_{p2}.weir.fst"
             for ref in REFERENCE_NAMES for (p1, p2) in POP_PAIR_TUPLES],
    output:
        summary=QC_OUTDIR / "per_ref_summary.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        refs=REFERENCE_NAMES,
        fastas=[REFS_NESTED[r]["fasta"] for r in REFERENCE_NAMES],
        fais=[f"{REFS_NESTED[r]['fasta']}.fai" for r in REFERENCE_NAMES],
        conspec_ref="conspec",
        pop_pairs=POP_PAIR_TUPLES,
        include_graph=INCLUDE_GRAPH,
        min_sv_size=SV_MIN_SIZE,
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=8000,
        cpus=1,
    script:
        "../scripts/per_ref_summary.py"


# ---------------------------------------------------------------------------
# Per-sample QC summary (report-only; failed samples are flagged but not
# auto-excluded from downstream analyses)
# ---------------------------------------------------------------------------
rule sample_qc_summary:
    """Aggregate per-sample metrics across all refs into a single QC TSV.
    Flags samples that violate any configured threshold but doesn't exclude
    them — the user inspects this table and decides what to drop."""
    input:
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv",
        heterozygosity=expand(
            HET_OUTDIR / "{ref}/per_sample_heterozygosity.tsv",
            ref=REFERENCE_NAMES,
        ),
    output:
        qc=QC_OUTDIR / "sample_qc.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        min_depth=QC_MIN_DEPTH,
        min_mapping_rate=QC_MIN_MAPPING_RATE,
        min_het_rate=QC_MIN_HET_RATE,
        max_het_rate=QC_MAX_HET_RATE,
        max_missing_rate=QC_MAX_MISSING_RATE,
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/sample_qc.py"


# ---------------------------------------------------------------------------
# Per-sample heterozygosity (bcftools stats -s)
# ---------------------------------------------------------------------------
rule heterozygosity_per_ref:
    """Per-sample non-reference het and hom counts from bcftools stats -s.
    Het rate is computed against the count of called sites for that sample
    (not the genome length) to make it depth-robust. Goes into the
    per-sample QC summary."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        tsv=HET_OUTDIR / "{ref}/per_sample_heterozygosity.tsv",
    conda:
        "../envs/bcftools.yaml"
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=24000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        # bcftools stats emits a PSC block (Per-Sample Counts) with columns:
        #   PSC  id  sample  nRefHom  nNonRefHom  nHets  nTransitions
        #   nTransversions  nIndels  averageDepth  nSingletons  nHapRef
        #   nHapAlt  nMissing
        # Pull sample, nRefHom, nNonRefHom, nHets, nMissing, avgDepth; compute
        # het_rate = nHets / (nRefHom + nNonRefHom + nHets) so missingness
        # doesn't drag down the denominator.
        {{
            echo -e "sample\tref\tn_ref_hom\tn_nonref_hom\tn_het\tn_missing\tmean_depth\thet_rate"
            bcftools stats -s - {input.vcf} \
              | awk -v ref={wildcards.ref} '
                  /^# PSC/ {{ next }}
                  /^PSC\t/ {{
                      sample = $3; nRefHom = $4; nNonRefHom = $5; nHet = $6
                      nMissing = $14; avgDepth = $10
                      called = nRefHom + nNonRefHom + nHet
                      het_rate = (called > 0) ? nHet / called : 0
                      printf "%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.6f\n",
                             sample, ref, nRefHom, nNonRefHom, nHet,
                             nMissing, avgDepth, het_rate
                  }}
              '
        }} > {output.tsv}
        """


# ---------------------------------------------------------------------------
# Per-sample relatedness (KING-robust kinship via vcftools --relatedness2)
# ---------------------------------------------------------------------------
rule relatedness_per_ref:
    """Pairwise kinship coefficients across the cohort using the KING-robust
    estimator. Output is a long-format TSV (INDV1\\tINDV2\\tRELATEDNESS_PHI).
    Phi ~0.25 = parent-offspring/full-sibs; ~0.125 = half-sibs/avuncular;
    < ~0.04 = unrelated."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        rel=RELATEDNESS_OUTDIR / "{ref}/relatedness.tsv",
    conda:
        "../envs/fst_afs.yaml"
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=64000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.rel})
        mkdir -p "$outdir"
        prefix="$outdir/relatedness"
        vcftools --gzvcf {input.vcf} --relatedness2 --out "$prefix" > /dev/null
        # vcftools writes <prefix>.relatedness2 with columns:
        #   INDV1 INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI
        mv "$prefix.relatedness2" {output.rel}
        """


rule relatedness_summary:
    """Cohort-wide relatedness table across all refs. One row per (ref,
    sample_pair), tagged by the maximum kinship that pair shows on any ref so
    related pairs surface even if one ref is noisy."""
    input:
        expand(RELATEDNESS_OUTDIR / "{ref}/relatedness.tsv", ref=REFERENCE_NAMES),
    output:
        summary=RELATEDNESS_OUTDIR / "relatedness_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    run:
        from pathlib import Path

        out_path = Path(output.summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rows = []
        for path in input:
            ref = Path(path).parts[-2]
            with open(path) as fh:
                header = fh.readline().strip().split("\t")
                try:
                    phi_idx = header.index("RELATEDNESS_PHI")
                    i1_idx = header.index("INDV1")
                    i2_idx = header.index("INDV2")
                except ValueError:
                    continue
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) <= phi_idx:
                        continue
                    if parts[i1_idx] == parts[i2_idx]:
                        continue
                    rows.append((ref, parts[i1_idx], parts[i2_idx], parts[phi_idx]))
        with out_path.open("w") as fh:
            fh.write("ref\tsample1\tsample2\trelatedness_phi\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")


# ---------------------------------------------------------------------------
# Nucleotide diversity (π) in sliding windows
# ---------------------------------------------------------------------------
rule pi_per_ref:
    """Per-population windowed nucleotide diversity. Runs vcftools --window-pi
    once per population (via --keep) and concatenates into one long-format table
    with a POP column. A single genome-wide file per ref is preserved (same
    path) for backward compatibility; per-population comparison happens in the
    plot. Pooling all populations would distort π/Tajima's D via structure, so
    diversity is measured within each population."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        PI_OUTDIR / "{ref}.windowed.pi"
    conda:
        "../envs/pi.yaml"
    params:
        window_size=PI_WINDOW_SIZE,
        window_step=PI_WINDOW_STEP,
        pop_samples=POP_SAMPLES,
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=8000,
        cpus=1
    script:
        "../scripts/compute_pi_per_pop.py"


# ---------------------------------------------------------------------------
# Tajima's D in non-overlapping windows (vcftools)
# ---------------------------------------------------------------------------
rule tajimas_d_per_ref:
    """Per-population Tajima's D in non-overlapping windows. Runs vcftools
    --TajimaD once per population (via --keep) and concatenates into one
    long-format table with a POP column. Per-population is essential here:
    pooling populations inflates Tajima's D through structure and masks the
    demographic signal we want to read per group."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        PI_OUTDIR / "{ref}.Tajima.D"
    conda:
        "../envs/pi.yaml"
    params:
        # vcftools --TajimaD only takes a single non-overlapping window size.
        window_size=PI_WINDOW_SIZE,
        pop_samples=POP_SAMPLES,
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=8000,
        cpus=1
    script:
        "../scripts/compute_tajimas_d_per_pop.py"


# ---------------------------------------------------------------------------
# Allelic balance (ref bias) per sample
# ---------------------------------------------------------------------------
rule allelic_balance_per_sample:
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/per_sample/{sample}.vcf.gz.tbi"
    output:
        summary=AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv",
        raw=AB_OUTDIR / "{ref}/{sample}.allelic_balance.raw.tsv"
    conda:
        "../envs/allelic_balance.yaml"
    params:
        bins=AB_BINS,
        ref=lambda wildcards: wildcards.ref,
        min_depth=AB_MIN_DEPTH,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    script:
        "../scripts/allelic_balance_per_sample.py"


rule allelic_balance_summary:
    input:
        expand(AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES)
    output:
        summary=AB_OUTDIR / "allelic_balance_summary.tsv"
    conda:
        "../envs/allelic_balance.yaml"
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


# ---------------------------------------------------------------------------
# Core vs accessory depth / missingness / MAPQ (augref only)
# ---------------------------------------------------------------------------
rule region_depth_missingness:
    """Per-sample depth, mean MAPQ, and GT missingness split core (non-SV_*)
    vs accessory (SV_*) on augref. Distinguishes whether a core-vs-accessory
    Tajima's D / FST gap reflects real shared variation (per user hypothesis)
    from a coverage/mappability artifact: if accessory depth/MAPQ is
    systematically lower, rare variants there are undercalled, which alone
    inflates Tajima's D and can suppress FST."""
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai",
        vcf=VC_OUTDIR / "bcftools/augref/per_sample/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/augref/per_sample/{sample}.vcf.gz.tbi",
        fai=f"{ALIGN_AUGREF}.fai",
    output:
        tsv=QC_OUTDIR / "region_qc/{sample}.region_qc.tsv",
    conda:
        "../envs/bcftools.yaml"
    params:
        sample=lambda wildcards: wildcards.sample,
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/region_depth_missingness.py"


rule region_qc_summary:
    """Cohort-wide core-vs-accessory depth/MAPQ/missingness table (augref)."""
    input:
        expand(QC_OUTDIR / "region_qc/{sample}.region_qc.tsv", sample=SHORT_SAMPLES),
    output:
        summary=QC_OUTDIR / "region_qc_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=2000,
        cpus=1,
    run:
        from pathlib import Path

        out_path = Path(output.summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rows, header = [], None
        for path in input:
            lines = Path(path).read_text().strip().splitlines()
            if not lines:
                continue
            if header is None:
                header = lines[0]
            rows.extend(lines[1:])
        if header is None:
            header = "sample\tregion\tn_contigs\tmean_depth\tmean_mapq\tn_sites\tmissing_rate"
        out_path.write_text("\n".join([header] + rows) + "\n")


# ---------------------------------------------------------------------------
# Runs of homozygosity (F_ROH) — REMOVED
# ---------------------------------------------------------------------------
# ROH was cut entirely: with SV_* accessory contigs excluded from ROH calling,
# augref ROH is identical to conspec by construction (the added sequence is
# never assessed), so augref added no information and the whole stage was
# redundant. The rules below are disabled (wrapped in a module-level string so
# Snakemake never registers them). To revive, unwrap this block and restore the
# config vars / target lines / summary-table columns noted in the main Snakefile.
_ROH_DISABLED = r'''
rule ref_lengths:
    input:
        fasta=lambda wildcards: _ref_fasta(wildcards),
        fai=lambda wildcards: _ref_fai(wildcards)
    output:
        lengths=ROH_OUTDIR / "lengths" / "{ref}.lengths.tsv"
    conda:
        "../envs/roh.yaml"
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


rule roh_inputs:
    """ROH input pair (BCF + AF table) for all refs. Uses the bcftools
    merged VCF with fill-tags to compute cohort AF from genotypes."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        bcf=ROH_OUTDIR / "{ref}/roh_input.bcf",
        csi=ROH_OUTDIR / "{ref}/roh_input.bcf.csi",
        freqs=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz",
        tbi=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz.tbi",
    conda:
        "../envs/roh.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=8000,
        cpus=4,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_OUTDIR}/{wildcards.ref}

        bcftools +fill-tags {input.vcf} -- -t AF \
            | bcftools view -Ob --threads {threads} -o {output.bcf}

        bcftools index {output.bcf}

        bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' {output.bcf} \
            | bgzip -c > {output.freqs}
        tabix -s1 -b2 -e2 {output.freqs}
        """


rule roh_per_sample:
    """Per-sample ROH from population AF + bcftools roh PL mode (paper-style
    workflow). Source of the BCF + AF table differs by ref: ANGSD on BAMs
    for linear refs, vg-call cohort VCF for mc_graph."""
    input:
        bcf=ROH_OUTDIR / "{ref}/roh_input.bcf",
        csi=ROH_OUTDIR / "{ref}/roh_input.bcf.csi",
        freqs=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz",
        tbi=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz.tbi",
        lengths=lambda wildcards: ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv",
    output:
        roh=ROH_OUTDIR / "{ref}/{sample}.roh.tsv",
        froh=ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv",
    conda:
        "../envs/roh.yaml"
    params:
        min_roh_length=ROH_MIN_LENGTH,
        recomb_rate=ROH_RECOMB_RATE,
        bcftools_args=ROH_BCFTOOLS_ARGS,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/roh_per_sample.py"


rule roh_summary:
    input:
        expand(ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv", ref=REFERENCE_NAMES, sample=SHORT_SAMPLES)
    output:
        summary=ROH_OUTDIR / "f_roh_summary.tsv"
    conda:
        "../envs/roh.yaml"
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
'''


# ---------------------------------------------------------------------------
# LD pruning + PCAngsd PCA + selection scan
# ---------------------------------------------------------------------------
rule split_ld_prune_contigs:
    """Split the merged VCF into three contig groups ahead of ld_prune_vcf /
    ld_prune_accessory_small: core (normal chromosomes), accessory_big (SV_*
    contigs long enough for real LD-pruning to be worth it), and
    accessory_small (SV_* contigs too short for that -- see ld_prune_vcf).

    SV_* contigs are flanked SV sequence (sv_calling.flank bp on each side,
    per config), so they are NOT guaranteed to carry only ~1 variant:
    multi-SNP/indel clusters within one contig are tightly linked (same short
    fragment, no realistic within-contig recombination), but real LD-pruning
    still matters once a contig is long enough for that linkage assumption to
    be worth checking rather than just assumed. accessory_big contigs
    (>= pcangsd.accessory_ld_min_len) go through the same plink2 LD-prune
    pass as core; accessory_small contigs go through
    ld_prune_accessory_small's collapse-to-1-variant-per-contig instead.

    For refs with no SV_* contigs (conspec, mc_graph) both accessory subsets
    are empty and downstream is equivalent to pruning core alone.
    """
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        core_vcf=temp(PCANGSD_OUTDIR / "{ref}/core_tmp.vcf.gz"),
        core_tbi=temp(PCANGSD_OUTDIR / "{ref}/core_tmp.vcf.gz.tbi"),
        accessory_small_vcf=temp(PCANGSD_OUTDIR / "{ref}/accessory_small_tmp.vcf.gz"),
        accessory_small_tbi=temp(PCANGSD_OUTDIR / "{ref}/accessory_small_tmp.vcf.gz.tbi"),
        orig_samples=temp(PCANGSD_OUTDIR / "{ref}/orig_samples.txt"),
    conda:
        "../envs/ld_prune.yaml"
    params:
        accessory_ld_min_len=PCANGSD_ACCESSORY_LD_MIN_LEN,
    resources:
        slurm_partition="short",
        runtime=240,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.core_vcf})
        mkdir -p "$outdir"

        # Assign unique IDs and remove duplicate positions.
        dedup_vcf="$outdir/dedup_tmp.vcf.gz"
        bcftools norm --rm-dup any -Oz -o "$dedup_vcf" {input.vcf}
        bcftools index -t "$dedup_vcf"

        # plink's --double-id writes VCF sample names as FID_IID =
        # "SAMPLE_SAMPLE". Save the original (clean) names now so downstream
        # rules can reheader back to them before recombining.
        bcftools query -l {input.vcf} > {output.orig_samples}

        # Contig name + length come from the VCF's own ##contig header lines
        # (no separate .fai needed -- this rule is generic over all refs).
        # A 3-column BED (CHROM/START/END) is the format bcftools view -R
        # reliably parses via a file; a bare one-column contig-name list
        # doesn't (see split_ldpruned_vcf_by_region for the same gotcha).
        core_bed="$outdir/core_contigs.bed"
        accessory_big_bed="$outdir/accessory_big_contigs.bed"
        accessory_small_bed="$outdir/accessory_small_contigs.bed"
        bcftools view -h "$dedup_vcf" \
            | grep '^##contig=' \
            | sed -E 's/.*ID=([^,>]+).*length=([0-9]+).*/\1\t\2/' \
            > "$outdir/contigs.tsv"
        awk -F'\t' '$1 !~ /^SV_/ {{print $1"\t0\t"$2}}' "$outdir/contigs.tsv" > "$core_bed"
        awk -F'\t' -v minlen={params.accessory_ld_min_len} \
            '$1 ~ /^SV_/ && $2 >= minlen {{print $1"\t0\t"$2}}' \
            "$outdir/contigs.tsv" > "$accessory_big_bed"
        awk -F'\t' -v minlen={params.accessory_ld_min_len} \
            '$1 ~ /^SV_/ && $2 < minlen {{print $1"\t0\t"$2}}' \
            "$outdir/contigs.tsv" > "$accessory_small_bed"

        # accessory_big sites are pruned together with core (same plink2
        # pass in ld_prune_vcf), so fold their BED into core's here.
        cat "$core_bed" "$accessory_big_bed" > "$outdir/core_plus_accessory_big.bed"
        bcftools view -R "$outdir/core_plus_accessory_big.bed" -Oz -o {output.core_vcf} "$dedup_vcf"
        bcftools index -t {output.core_vcf}
        if [ -s "$accessory_small_bed" ]; then
            bcftools view -R "$accessory_small_bed" -Oz -o {output.accessory_small_vcf} "$dedup_vcf"
        else
            bcftools view -h "$dedup_vcf" -Oz -o {output.accessory_small_vcf}
        fi
        bcftools index -t {output.accessory_small_vcf}

        rm -f "$dedup_vcf" "$dedup_vcf.tbi" "$outdir/contigs.tsv" \
              "$core_bed" "$accessory_big_bed" "$accessory_small_bed" \
              "$outdir/core_plus_accessory_big.bed"
        """


rule ld_prune_vcf:
    """LD-prune the core+accessory_big VCF with plink2. accessory_small sites
    (too short for LD-pruning to be worth it) are handled separately by
    ld_prune_accessory_small and recombined here."""
    input:
        core_vcf=PCANGSD_OUTDIR / "{ref}/core_tmp.vcf.gz",
        core_tbi=PCANGSD_OUTDIR / "{ref}/core_tmp.vcf.gz.tbi",
        accessory_small_pruned=PCANGSD_OUTDIR / "{ref}/pruned_accessory_small.vcf.gz",
        accessory_small_pruned_tbi=PCANGSD_OUTDIR / "{ref}/pruned_accessory_small.vcf.gz.tbi",
        orig_samples=PCANGSD_OUTDIR / "{ref}/orig_samples.txt",
    output:
        vcf=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz.tbi",
        pruned_flag=PCANGSD_OUTDIR / "{ref}/pruned.flag",
    conda:
        "../envs/ld_prune.yaml"
    params:
        maf=PCANGSD_MAF,
        window=PCANGSD_LD_WINDOW,
        step=PCANGSD_LD_STEP,
        r2=PCANGSD_LD_R2,
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.vcf})
        mkdir -p "$outdir"

        # plink2 refuses --indep-pairwise below 50 samples/founders: the r²
        # estimates it prunes on are too noisy to trust at that cohort size.
        # Below that threshold, skip LD pruning rather than override the
        # guardrail (--bad-ld) -- keep MAF/SNP-only filtering only, and leave
        # correlated markers in. PCA/selection-scan results from an unpruned
        # small cohort should be caveated accordingly downstream.
        n_samples=$(bcftools query -l {input.core_vcf} | wc -l)
        if [ "$n_samples" -ge 50 ]; then
            echo "true" > {output.pruned_flag}
            # Step 1: compute LD prune list
            plink2 \
                --vcf {input.core_vcf} \
                --double-id --allow-extra-chr \
                --set-missing-var-ids @:#_\$r_\$a \
                --snps-only just-acgt \
                --maf {params.maf} \
                --memory {resources.mem_mb} \
                --indep-pairwise {params.window} {params.step} {params.r2} \
                --out "$outdir/prune_tmp"

            # Step 2: extract pruned sites → VCF
            plink2 \
                --vcf {input.core_vcf} \
                --double-id --allow-extra-chr \
                --set-missing-var-ids @:#_\$r_\$a \
                --snps-only just-acgt \
                --maf {params.maf} \
                --memory {resources.mem_mb} \
                --extract "$outdir/prune_tmp.prune.in" \
                --export vcf \
                --out "$outdir/pruned_core_tmp"
        else
            echo "WARNING: $n_samples samples (<50) for {wildcards.ref} -- skipping LD pruning, MAF/SNP-only filter only" >&2
            echo "false" > {output.pruned_flag}
            plink2 \
                --vcf {input.core_vcf} \
                --double-id --allow-extra-chr \
                --set-missing-var-ids @:#_\$r_\$a \
                --snps-only just-acgt \
                --maf {params.maf} \
                --memory {resources.mem_mb} \
                --export vcf \
                --out "$outdir/pruned_core_tmp"
        fi
        bgzip -f "$outdir/pruned_core_tmp.vcf"
        bcftools reheader -s {input.orig_samples} \
            -o "$outdir/pruned_core_tmp.reheader.vcf.gz" "$outdir/pruned_core_tmp.vcf.gz"
        mv "$outdir/pruned_core_tmp.reheader.vcf.gz" "$outdir/pruned_core_tmp.vcf.gz"
        bcftools index -t "$outdir/pruned_core_tmp.vcf.gz"

        # Recombine with accessory_small (already filtered, collapsed to 1
        # variant/contig, and reheadered by ld_prune_accessory_small).
        bcftools concat -a \
            "$outdir/pruned_core_tmp.vcf.gz" {input.accessory_small_pruned} \
            -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}

        rm -f "$outdir/pruned_core_tmp.vcf.gz" "$outdir/pruned_core_tmp.vcf.gz.tbi"
        rm -f "$outdir/prune_tmp".* "$outdir/pruned_core_tmp.log"
        """


rule ld_prune_accessory_small:
    """Filter and collapse accessory_small SV_* sites (see
    split_ld_prune_contigs and scripts/ld_prune_accessory_small.py for
    rationale) into 1 representative variant per contig. bcftools-only
    (no plink2), so sample names are untouched -- no reheader needed before
    ld_prune_vcf recombines this with the pruned core+accessory_big VCF."""
    input:
        vcf=PCANGSD_OUTDIR / "{ref}/accessory_small_tmp.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/accessory_small_tmp.vcf.gz.tbi",
    output:
        vcf=PCANGSD_OUTDIR / "{ref}/pruned_accessory_small.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/pruned_accessory_small.vcf.gz.tbi",
    conda:
        "../envs/ld_prune.yaml"
    params:
        min_dp=PCANGSD_ACCESSORY_MIN_DP,
        min_gq=PCANGSD_ACCESSORY_MIN_GQ,
        maf=PCANGSD_MAF,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/ld_prune_accessory_small.py"


rule pcangsd:
    """Run PCAngsd on LD-pruned VCF to get covariance matrix and selection scores."""
    input:
        vcf=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz.tbi",
        pruned_flag=PCANGSD_OUTDIR / "{ref}/pruned.flag",
    output:
        cov=PCANGSD_OUTDIR / "{ref}/pcangsd.cov",
        selection=PCANGSD_OUTDIR / "{ref}/pcangsd.selection",
        samples=PCANGSD_OUTDIR / "{ref}/samples.txt",
        snp_coords=PCANGSD_OUTDIR / "{ref}/snp_coords.tsv",
    conda:
        "../envs/pcangsd.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=480,
        mem_mb=16000,
        cpus=4,
    params:
        prefix=lambda wildcards: str(PCANGSD_OUTDIR / wildcards.ref / "pcangsd"),
        n_pcs=PCANGSD_N_PCS,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        tmpdir=$(mktemp -d -p $(dirname {params.prefix}))
        trap "rm -rf $tmpdir" EXIT

        # pcangsd >=1.x dropped --vcf / --n_eig and now consumes PLINK BED via
        # -p PREFIX. Convert with plink2 first. augref's LD-pruned VCF carries
        # one representative variant per accessory SV_* contig, and there can be
        # far more than plink2's ~65k distinct-nonstandard-contig cap, so the
        # helper rewrites CHROM/POS into a bounded set of placeholder
        # chromosomes before --make-bed (see scripts/vcf_to_plink_bucketed.py).
        # It returns the plink prefix and a real-coords table (true contig
        # names, in .bim row order) on two lines.
        conv_out=$(python {workflow.basedir}/scripts/vcf_to_plink_bucketed.py \
            --vcf {input.vcf} --tmpdir "$tmpdir" --memory-mb {resources.mem_mb})
        plink_prefix=$(echo "$conv_out" | sed -n '1p')
        real_coords=$(echo "$conv_out" | sed -n '2p')

        # Sample names for the PCA plot. Take them from the VCF (bcftools gives
        # them in the same order plink uses). Older LD-pruned VCFs carry plink's
        # doubled "SAMPLE_SAMPLE" names (from --double-id upstream); collapse any
        # "NAME_NAME" back to "NAME" so labels match the manifest sample_ids.
        # (New ld_prune_vcf runs already reheader to clean names, so this is a
        # no-op there.)
        bcftools query -l {input.vcf} \
            | sed -E 's/^(.*)_\1$/\1/' > {output.samples}

        pcangsd \
            -p "$plink_prefix" \
            -t {threads} \
            -e {params.n_pcs} \
            --selection \
            --sites-save \
            -o {params.prefix}

        # pcangsd --sites-save writes a boolean mask (0/1 per BIM site).
        # Subset the real-coords table (true SV_* names, in .bim row order --
        # the chunk placeholders fed to plink2 never appear here) to only the
        # sites pcangsd kept. real_coords is already CHROM\tPOS\tID.
        paste "$real_coords" {params.prefix}.sites \
            | awk -F'\t' '$4 == 1 {{OFS="\t"; print $1, $2, $3}}' \
            > {output.snp_coords}
        """


rule split_ldpruned_vcf_by_region:
    """Subset the augref LD-pruned VCF into core (non-SV_*) and accessory
    (SV_*) site sets for a separate PCA per region. If population structure
    (the 3 short-read populations) shows up clearly in the core-only PCA but
    collapses in the accessory-only PCA, that supports the "accessory carries
    variation shared across populations" hypothesis directly, rather than
    inferring it from Tajima's D / FST alone."""
    input:
        vcf=PCANGSD_OUTDIR / "augref/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "augref/merged.ldpruned.vcf.gz.tbi",
        fai=f"{ALIGN_AUGREF}.fai",
    output:
        vcf=PCANGSD_OUTDIR / "augref_{region}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "augref_{region}/merged.ldpruned.vcf.gz.tbi",
    conda:
        "../envs/ld_prune.yaml"
    wildcard_constraints:
        region="core|accessory",
    resources:
        slurm_partition="short",
        runtime=30,
        mem_mb=4000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})
        bed_file=$(dirname {output.vcf})/keep_contigs.bed
        if [ "{wildcards.region}" = "core" ]; then
            awk -F'\t' '$1 !~ /^SV_/ {{print $1"\t0\t"$2}}' {input.fai} > "$bed_file"
        else
            awk -F'\t' '$1 ~ /^SV_/ {{print $1"\t0\t"$2}}' {input.fai} > "$bed_file"
        fi
        # A comma-joined -t/-r argument string blows the shell's ARG_MAX with
        # augref's thousands of SV_* contigs ("Argument list too long"), and
        # -R/-T with a bare one-column contig-name file fails to parse
        # ("Could not parse ... using the columns 1,2[,-1]") — bcftools wants
        # CHROM+POS columns there. A 3-column BED (CHROM/START/END) is the
        # documented format -R accepts via a file, so use that instead.
        bcftools view -R "$bed_file" {input.vcf} \
            -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        rm -f "$bed_file"
        """


rule pcangsd_by_region:
    """PCAngsd restricted to augref core-only or accessory-only sites. Same
    body as the main `pcangsd` rule (plink BED conversion, --selection,
    --sites-save), applied to the region-subset VCF from
    split_ldpruned_vcf_by_region instead of the whole-genome LD-pruned VCF."""
    input:
        vcf=PCANGSD_OUTDIR / "augref_{region}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "augref_{region}/merged.ldpruned.vcf.gz.tbi",
    output:
        cov=PCANGSD_OUTDIR / "augref_{region}/pcangsd.cov",
        samples=PCANGSD_OUTDIR / "augref_{region}/samples.txt",
    conda:
        "../envs/pcangsd.yaml"
    wildcard_constraints:
        region="core|accessory",
    threads: 4
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=16000,
        cpus=4,
    params:
        prefix=lambda wildcards: str(PCANGSD_OUTDIR / f"augref_{wildcards.region}" / "pcangsd"),
        n_pcs=PCANGSD_N_PCS,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        tmpdir=$(mktemp -d -p $(dirname {params.prefix}))
        trap "rm -rf $tmpdir" EXIT

        # Convert to PLINK BED, remapping CHROM/POS into bounded placeholder
        # chunks so plink2 stays under its ~65k distinct-contig cap (the
        # augref_accessory region is *all* SV_* contigs, so this is essential
        # here). See
        # scripts/vcf_to_plink_bucketed.py; it returns the plink prefix and a
        # real-coords table (unused here -- no snp_coords output).
        conv_out=$(python {workflow.basedir}/scripts/vcf_to_plink_bucketed.py \
            --vcf {input.vcf} --tmpdir "$tmpdir" --memory-mb {resources.mem_mb})
        plink_prefix=$(echo "$conv_out" | sed -n '1p')

        bcftools query -l {input.vcf} \
            | sed -E 's/^(.*)_\1$/\1/' > {output.samples}

        pcangsd \
            -p "$plink_prefix" \
            -t {threads} \
            -e {params.n_pcs} \
            -o {params.prefix}
        """


rule fst_outliers:
    """
    Flag per-SNP selection outliers from PCAngsd chi-squared scores with
    Benjamini-Hochberg FDR correction. Outputs a TSV of all loci with
    coordinates, chi-squared score, p-value, q-value, and outlier flag.
    """
    input:
        selection=PCANGSD_OUTDIR / "{ref}/pcangsd.selection",
        snp_coords=PCANGSD_OUTDIR / "{ref}/snp_coords.tsv",
    output:
        outliers=PCANGSD_OUTDIR / "{ref}/fst_outliers.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        fdr=0.05,
    resources:
        slurm_partition="short",
        runtime=60,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/fst_outliers.py"
