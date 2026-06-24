# ============================================================================
# Population Genetics Rules
# ============================================================================
# Variables (VC_OUTDIR, FST_OUTDIR, PI_OUTDIR, AB_OUTDIR, ROH_OUTDIR, etc.)
# and helper functions (_ref_fasta, _ref_fai) are defined in the
# main Snakefile and are available via Snakemake's include mechanism.


def _ref_lengths_path(wildcards):
    if wildcards.ref in ROH_GENOME_LENGTHS:
        return ROH_GENOME_LENGTHS[wildcards.ref]
    return str(ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv")


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
        mem_mb=4000,
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
        mem_mb=4000,
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
        runtime=120,
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
rule per_ref_summary:
    """One row per reference with the headline comparison metrics:
    SNP counts (raw + filtered), SV count (mc_graph only), mean mapping rate,
    mean mapping quality, mean π, mean F_ROH, genome length, size vs conspec.
    Mirrors Jeon et al. 2026 Table 1."""
    input:
        vcfs=expand(
            VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
            ref=REFERENCE_NAMES,
        ),
        sv_vcfs=expand(VC_OUTDIR / "sv/vg/{sample}.vcf.gz", sample=SHORT_SAMPLES),
        metrics=METRICS_OUTDIR / "alignment_metrics.tsv",
        froh=ROH_OUTDIR / "f_roh_summary.tsv",
        pi=expand(PI_OUTDIR / "{ref}.windowed.pi", ref=REFERENCE_NAMES),
    output:
        summary=QC_OUTDIR / "per_ref_summary.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        refs=REFERENCE_NAMES,
        fastas=[REFS_NESTED[r]["fasta"] for r in REFERENCE_NAMES],
        fais=[f"{REFS_NESTED[r]['fasta']}.fai" for r in REFERENCE_NAMES],
        conspec_ref="conspec",
    resources:
        slurm_partition="short",
        runtime=120,
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
        froh_summary=ROH_OUTDIR / "f_roh_summary.tsv",
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
        mem_mb=4000,
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
        mem_mb=4000,
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
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        PI_OUTDIR / "{ref}.windowed.pi"
    conda:
        "../envs/pi.yaml"
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


# ---------------------------------------------------------------------------
# Tajima's D in non-overlapping windows (vcftools)
# ---------------------------------------------------------------------------
rule tajimas_d_per_ref:
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        PI_OUTDIR / "{ref}.Tajima.D"
    conda:
        "../envs/pi.yaml"
    params:
        # vcftools --TajimaD only takes a single non-overlapping window size.
        window_size=PI_WINDOW_SIZE
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
          --TajimaD {params.window_size} \
          --out {PI_OUTDIR}/{wildcards.ref} > /dev/null
        """


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
# Runs of homozygosity (F_ROH)
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# LD pruning + PCAngsd PCA + selection scan
# ---------------------------------------------------------------------------
rule ld_prune_vcf:
    """LD-prune the merged VCF with plink for use in PCAngsd."""
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi",
    output:
        vcf=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz.tbi",
    conda:
        "../envs/ld_prune.yaml"
    params:
        maf=PCANGSD_MAF,
        window=PCANGSD_LD_WINDOW,
        step=PCANGSD_LD_STEP,
        r2=PCANGSD_LD_R2,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.vcf})
        mkdir -p "$outdir"

        # Pre-filter: assign unique IDs and remove duplicate positions with bcftools
        dedup_vcf="$outdir/dedup_tmp.vcf.gz"
        bcftools norm --rm-dup any -Oz -o "$dedup_vcf" {input.vcf}
        bcftools index -t "$dedup_vcf"

        # Step 1: compute LD prune list
        plink \
            --vcf "$dedup_vcf" \
            --double-id --allow-extra-chr \
            --set-missing-var-ids @:#_$$1_$$2 \
            --snps-only just-acgt \
            --maf {params.maf} \
            --memory {resources.mem_mb} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out "$outdir/prune_tmp"

        # Step 2: extract pruned sites → VCF
        plink \
            --vcf "$dedup_vcf" \
            --double-id --allow-extra-chr \
            --set-missing-var-ids @:#_$$1_$$2 \
            --snps-only just-acgt \
            --maf {params.maf} \
            --memory {resources.mem_mb} \
            --extract "$outdir/prune_tmp.prune.in" \
            --recode vcf \
            --out "$outdir/pruned_tmp"

        bgzip -f "$outdir/pruned_tmp.vcf"
        mv "$outdir/pruned_tmp.vcf.gz" {output.vcf}
        tabix -p vcf {output.vcf}

        rm -f "$dedup_vcf" "$dedup_vcf.tbi"
        rm -f "$outdir/prune_tmp".* "$outdir/pruned_tmp.log" "$outdir/pruned_tmp.nosex"
        """


rule pcangsd:
    """Run PCAngsd on LD-pruned VCF to get covariance matrix and selection scores."""
    input:
        vcf=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz.tbi",
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
        runtime=120,
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
        # -p PREFIX. Convert with plink first.
        # --const-fid 0: VCF has no family info; assign FID=0 to every sample.
        # --allow-extra-chr: tolerate non-numeric contig names (NC_*, NW_*).
        plink --vcf {input.vcf} --make-bed --const-fid 0 --allow-extra-chr \
              --memory {resources.mem_mb} \
              --out "$tmpdir/plink"

        awk '{{print $2}}' "$tmpdir/plink.fam" > {output.samples}
        # Save all BIM coords; filter to pcangsd-kept sites after --sites-save.
        awk -v OFS='\t' '{{print $1, $4, $2}}' "$tmpdir/plink.bim" > "$tmpdir/all_coords.tsv"

        pcangsd \
            -p "$tmpdir/plink" \
            -t {threads} \
            -e {params.n_pcs} \
            --selection \
            --sites-save \
            -o {params.prefix}

        # pcangsd --sites-save writes a boolean mask (0/1 per BIM site).
        # Subset coords to only the sites pcangsd kept.
        paste "$tmpdir/all_coords.tsv" {params.prefix}.sites \
            | awk -F'\t' '$4 == 1 {{OFS="\t"; print $1, $2, $3}}' \
            > {output.snp_coords}
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
