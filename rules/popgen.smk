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
# Allele Frequency Spectrum
# ---------------------------------------------------------------------------
rule afs_per_ref:
    input:
        vcf=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/{ref}/combined/merged.vcf.gz.tbi"
    output:
        afs=FST_OUTDIR / "afs/{ref}.afs.tsv"
    conda:
        "../envs/fst_afs.yaml"
    params:
        bins=AFS_BINS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
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

        pop1_tmp=$(mktemp)
        pop2_tmp=$(mktemp)
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


rule angsd_roh_genotypes:
    """Genotype likelihoods from ANGSD across the cohort, producing a BCF
    with population allele frequencies in INFO/AF (per Jeon et al. 2026).
    bcftools roh on this BCF (downstream) uses --AF-file with these AFs and
    runs in PL mode, which is what the paper does and what works robustly
    at low/medium coverage. Runs only on linear refs that have a real BAM —
    mc_graph routes through mc_graph_roh_inputs instead, since vg call
    works directly on the GAM and there's no BAM to feed ANGSD."""
    input:
        bams=_bams_for_ref,
        bais=_bais_for_ref,
        fasta=lambda wildcards: _ref_fasta(wildcards),
        fai=lambda wildcards: _ref_fai(wildcards),
    output:
        bcf=ROH_OUTDIR / "{ref}/roh_input.bcf",
        csi=ROH_OUTDIR / "{ref}/roh_input.bcf.csi",
        freqs=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz",
        tbi=ROH_OUTDIR / "{ref}/roh_input.freqs.tab.gz.tbi",
    conda:
        "../envs/roh.yaml"
    threads: 8
    resources:
        slurm_partition="long",
        runtime=720,
        mem_mb=16000,
        cpus=8,
    wildcard_constraints:
        # mc_graph has no BAM; its ROH inputs come from mc_graph_roh_inputs.
        ref="|".join(_BCFTOOLS_REFS) if _BCFTOOLS_REFS else "x^"
    params:
        # Paper used setMinDepth=100, setMaxDepth=375 for ~25 samples — i.e.
        # ~4× and ~15× the cohort sample count. Scale with current cohort.
        min_depth=lambda wc: 4 * len(SHORT_SAMPLES),
        max_depth=lambda wc: 15 * len(SHORT_SAMPLES),
        # 80% of samples must have data at a site for it to be used.
        min_ind=lambda wc: max(1, int(0.8 * len(SHORT_SAMPLES))),
        prefix=lambda wc: str(ROH_OUTDIR / wc.ref / "angsd"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_OUTDIR}/{wildcards.ref}
        bam_list=$(mktemp)
        trap "rm -f $bam_list" EXIT
        printf '%s\n' {input.bams} > "$bam_list"

        # only_proper_pairs=0: graph-derived alignments (in case any survive
        # here) don't always set proper-pair flags. Setting this to 1 can
        # silently drop everything.
        angsd -bam "$bam_list" \
              -ref {input.fasta} -anc {input.fasta} \
              -GL 1 -snp_pval 1e-6 \
              -dobcf 1 -dopost 1 -domajorminor 5 -domaf 1 -docounts 1 \
              -minQ 30 -minMapQ 20 \
              -minInd {params.min_ind} \
              -setMinDepth {params.min_depth} -setMaxDepth {params.max_depth} \
              -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 \
              -baq 2 -C 50 \
              -P {threads} \
              -out {params.prefix}

        # ANGSD writes <prefix>.bcf — move to the canonical roh_input name.
        if [ -f "{params.prefix}.bcf" ]; then
            mv "{params.prefix}.bcf" {output.bcf}
        fi

        bcftools index {output.bcf}

        bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' {output.bcf} \
            | bgzip -c > {output.freqs}
        tabix -s1 -b2 -e2 {output.freqs}
        """


rule mc_graph_roh_inputs:
    """ROH input pair (BCF + AF table) for mc_graph, derived from the
    cohort vg call merged VCF rather than from BAMs. vg call's per-sample
    GTs are already in conspec coordinates, so AF computed across samples
    here is directly comparable to the ANGSD AF used for linear refs."""
    input:
        vcf=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz.tbi",
    output:
        bcf=ROH_OUTDIR / "mc_graph/roh_input.bcf",
        csi=ROH_OUTDIR / "mc_graph/roh_input.bcf.csi",
        freqs=ROH_OUTDIR / "mc_graph/roh_input.freqs.tab.gz",
        tbi=ROH_OUTDIR / "mc_graph/roh_input.freqs.tab.gz.tbi",
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
        mkdir -p {ROH_OUTDIR}/mc_graph

        # bcftools +fill-tags adds INFO/AF computed from cohort GTs so the
        # downstream --AF-file path matches what ANGSD provides for linear
        # refs.
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
        "../envs/pcangsd.yaml"
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
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        # pcangsd >=1.x dropped --vcf / --n_eig and now consumes PLINK BED via
        # -p PREFIX. Convert with plink first.
        # --const-fid 0: VCF has no family info; assign FID=0 to every sample.
        # --allow-extra-chr: tolerate non-numeric contig names (NC_*, NW_*).
        plink --vcf {input.vcf} --make-bed --const-fid 0 --allow-extra-chr \
              --out "$tmpdir/plink"

        # Derive samples.txt and snp_coords.tsv from the PLINK files, NOT the
        # VCF — pcangsd's row/column order matches the .fam/.bim, which can
        # differ from the input VCF if plink reordered chromosomes.
        awk '{{print $2}}' "$tmpdir/plink.fam" > {output.samples}
        awk -v OFS='\t' '{{print $1, $4, $2}}' "$tmpdir/plink.bim" > {output.snp_coords}

        pcangsd \
            -p "$tmpdir/plink" \
            -t {threads} \
            -e {params.n_pcs} \
            --selection \
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
