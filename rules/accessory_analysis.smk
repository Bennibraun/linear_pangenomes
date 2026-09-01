# ============================================================================
# Accessory-sequence deep-dive rules (E1/E2/E3/E5)
# ============================================================================
# All rules here are LEAVES that consume already-materialised outputs
# (the combined merged VCF, the region-split VCFs, and the region PCAngsd
# covariances). Nothing here re-runs alignment, calling, LD-pruning, or the
# graph -- so adding them to `rule all` triggers only cheap downstream work.
#
# Motivation (see manuscript/RESPONSE_TO_NOTES.md): whole-genome FST/PCA
# concordance between augref and the graph is largely INHERITED from the shared
# core (only ~3.7% of augref sites are accessory). The claim that augmentation
# adds real, graph-concordant signal has to be made on ACCESSORY-ONLY sites,
# which conspec cannot inherit (conspec has no accessory sites at all). These
# rules produce the accessory-restricted statistics for that comparison.
#
# Variables (FST_OUTDIR, PCANGSD_OUTDIR, VC_OUTDIR, ALIGN_OUTDIR, POP_SAMPLES,
# POP_NAMES, POP_PAIR_TUPLES, SHORT_SAMPLES, ALIGN_AUGREF, FST_WINDOW_ARGS,
# AFS_BINS, config) come from the main Snakefile.

# Which contig regions to compute accessory-restricted stats over. "core" =
# non-SV_* contigs, "accessory" = SV_* contigs. augref only (conspec/hetspec
# have no SV_* contigs; mc_graph's accessory partition is a separate TODO that
# needs the graph, not a linear BED filter).
ACCESSORY_REGIONS = ["core", "accessory"]


# ---------------------------------------------------------------------------
# E1 (part 1): region-restricted pairwise FST on the FULL merged VCF
# ---------------------------------------------------------------------------
rule fst_by_region_pair:
    """Pairwise Weir-Cockerham FST restricted to core or accessory (SV_*)
    contigs, from augref's already-built combined merged VCF. Mirrors
    fst_per_ref_pair but adds a BED restriction (same core/accessory split the
    existing split_ldpruned_vcf_by_region rule uses). Leaf rule: the merged VCF
    is the only (already-materialised) heavy input."""
    input:
        vcf=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz.tbi",
        fai=f"{ALIGN_AUGREF}.fai",
    output:
        fst=FST_OUTDIR / "fst_by_region/{region}/{pop1}_vs_{pop2}.weir.fst",
    conda:
        "../envs/fst_afs.yaml"
    wildcard_constraints:
        region="core|accessory",
    params:
        window_args=FST_WINDOW_ARGS,
        pop1_samples=lambda w: " ".join(POP_SAMPLES[w.pop1]),
        pop2_samples=lambda w: " ".join(POP_SAMPLES[w.pop2]),
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.fst})
        mkdir -p "$outdir"
        prefix="${{outdir}}/{wildcards.pop1}_vs_{wildcards.pop2}"

        # 3-column BED of the contigs for this region (see
        # split_ldpruned_vcf_by_region for why a BED, not a name list).
        bed="${{outdir}}/{wildcards.pop1}_vs_{wildcards.pop2}.region.bed"
        if [ "{wildcards.region}" = "core" ]; then
            awk -F'\t' '$1 !~ /^(SV_|UNMAP_)/ {{print $1"\t0\t"$2}}' {input.fai} > "$bed"
        else
            awk -F'\t' '$1 ~ /^(SV_|UNMAP_)/ {{print $1"\t0\t"$2}}' {input.fai} > "$bed"
        fi

        region_vcf="${{outdir}}/{wildcards.pop1}_vs_{wildcards.pop2}.{wildcards.region}.vcf"
        pop1_tmp=$(mktemp -p "$outdir")
        pop2_tmp=$(mktemp -p "$outdir")
        trap "rm -f $pop1_tmp $pop2_tmp $bed $region_vcf" EXIT
        for s in {params.pop1_samples}; do echo "$s"; done > "$pop1_tmp"
        for s in {params.pop2_samples}; do echo "$s"; done > "$pop2_tmp"

        # Restrict to the region contigs first (bcftools view -R is strict:
        # a real failure here aborts the rule). Materialise a small region VCF
        # rather than piping, so vcftools failure and bcftools failure can be
        # told apart.
        bcftools view -R "$bed" {input.vcf} -Ov -o "$region_vcf"

        # vcftools may legitimately exit non-zero / write nothing when a pair has
        # no usable sites in this region -- tolerate ONLY vcftools here, then
        # guarantee a headed (possibly empty) output for the summary rule.
        vcftools --vcf "$region_vcf" \
              --weir-fst-pop "$pop1_tmp" \
              --weir-fst-pop "$pop2_tmp" \
              {params.window_args} \
              --out "$prefix" > /dev/null 2>&1 || true

        if [ ! -f "${{prefix}}.weir.fst" ]; then
            printf 'CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\n' > "${{prefix}}.weir.fst"
        fi
        """


rule fst_by_region_summary:
    """Collapse per-pair region FST into one tidy table:
    region, pop1, pop2, mean_fst, n_sites. Feeds the accessory-vs-graph
    comparison figure. Pure aggregation, trivially cheap."""
    input:
        fst=lambda w: [
            FST_OUTDIR / f"fst_by_region/{region}/{p1}_vs_{p2}.weir.fst"
            for region in ACCESSORY_REGIONS
            for (p1, p2) in POP_PAIR_TUPLES
        ],
    output:
        summary=FST_OUTDIR / "fst_by_region/fst_by_region_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        import re
        from pathlib import Path

        import numpy as np
        import pandas as pd

        na = ["-nan", "nan", "-NaN", "NaN", "-inf", "inf"]
        rows = []
        for path in input.fst:
            p = Path(path)
            region = p.parent.name
            m = re.match(r"(.+)_vs_(.+)\.weir\.fst$", p.name)
            pop1, pop2 = m.group(1), m.group(2)
            df = pd.read_csv(p, sep="\t", na_values=na)
            col = (
                "WEIGHTED_FST" if "WEIGHTED_FST" in df.columns
                else "WEIR_AND_COCKERHAM_FST"
            )
            if col in df.columns:
                v = pd.to_numeric(df[col], errors="coerce").dropna().clip(lower=0)
            else:
                v = pd.Series([], dtype=float)
            rows.append({
                "region": region, "pop1": pop1, "pop2": pop2,
                "mean_fst": float(v.mean()) if len(v) else np.nan,
                "n_sites": int(len(v)),
            })
        out = pd.DataFrame(rows).sort_values(["region", "pop1", "pop2"])
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output.summary, sep="\t", index=False, float_format="%.6f")


# ---------------------------------------------------------------------------
# E1 (part 2): the GRAPH arm -- mc_graph structural (SV) vs SNP FST
# ---------------------------------------------------------------------------
# augref's accessory = its SV_* contigs. The graph's equivalent isn't a BED
# filter (mc_graph is surjected to conspec coords, no SV_* contigs); its
# "accessory" content is its STRUCTURAL variation -- the SV-sized records that
# are exactly what the graph adds over a linear SNP-only view. So we split
# mc_graph's merged VCF by record size: variant="sv" (max allele-length change
# >= min_sv_size) vs "snp" (the rest), and compute pairwise FST for each. If
# accessory-augref FST tracks sv-mc_graph FST the way core tracks snp, that's
# the accessory signal being recovered by the linear representation.
GRAPH_VARIANT_CLASSES = ["sv", "snp"]


rule fst_mc_graph_by_variant_pair:
    """Pairwise FST on mc_graph's merged VCF, split into SV vs SNP records by
    allele-length change. Leaf rule off the already-built mc_graph merged VCF."""
    input:
        vcf=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/mc_graph/combined/merged.vcf.gz.tbi",
    output:
        fst=FST_OUTDIR / "fst_mc_graph_variant/{vclass}/{pop1}_vs_{pop2}.weir.fst",
    conda:
        "../envs/fst_afs.yaml"
    wildcard_constraints:
        vclass="sv|snp",
    params:
        window_args=FST_WINDOW_ARGS,
        min_sv=SV_MIN_SIZE,
        pop1_samples=lambda w: " ".join(POP_SAMPLES[w.pop1]),
        pop2_samples=lambda w: " ".join(POP_SAMPLES[w.pop2]),
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    shell:
        r"""
        set -euo pipefail
        outdir=$(dirname {output.fst})
        mkdir -p "$outdir"
        prefix="${{outdir}}/{wildcards.pop1}_vs_{wildcards.pop2}"
        sub="${{outdir}}/{wildcards.pop1}_vs_{wildcards.pop2}.{wildcards.vclass}.vcf"
        pop1_tmp=$(mktemp -p "$outdir"); pop2_tmp=$(mktemp -p "$outdir")
        trap "rm -f $pop1_tmp $pop2_tmp $sub" EXIT
        for s in {params.pop1_samples}; do echo "$s"; done > "$pop1_tmp"
        for s in {params.pop2_samples}; do echo "$s"; done > "$pop2_tmp"

        # Split by allele-length change vs REF. Normalise to biallelic first so
        # strlen(ALT)/strlen(REF) is one unambiguous value per record; -i keeps
        # SVs (indel length >= min_sv), -e keeps everything else (SNPs/small).
        len_expr='abs(strlen(ALT)-strlen(REF))>={params.min_sv}'
        if [ "{wildcards.vclass}" = "sv" ]; then
            bcftools norm -m- {input.vcf} -Ou | bcftools view -i "$len_expr" -Ov -o "$sub"
        else
            bcftools norm -m- {input.vcf} -Ou | bcftools view -e "$len_expr" -Ov -o "$sub"
        fi

        vcftools --vcf "$sub" \
            --weir-fst-pop "$pop1_tmp" --weir-fst-pop "$pop2_tmp" \
            {params.window_args} --out "$prefix" > /dev/null 2>&1 || true
        if [ ! -f "${{prefix}}.weir.fst" ]; then
            printf 'CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\n' > "${{prefix}}.weir.fst"
        fi
        """


rule fst_mc_graph_by_variant_summary:
    """Aggregate mc_graph SV/SNP pairwise FST into one table (same shape as
    fst_by_region_summary, so the accessory-vs-graph figure can join them)."""
    input:
        fst=lambda w: [
            FST_OUTDIR / f"fst_mc_graph_variant/{vc}/{p1}_vs_{p2}.weir.fst"
            for vc in GRAPH_VARIANT_CLASSES
            for (p1, p2) in POP_PAIR_TUPLES
        ],
    output:
        summary=FST_OUTDIR / "fst_mc_graph_variant/fst_mc_graph_variant_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        import re
        from pathlib import Path

        import numpy as np
        import pandas as pd

        na = ["-nan", "nan", "-NaN", "NaN", "-inf", "inf"]
        rows = []
        for path in input.fst:
            p = Path(path)
            vclass = p.parent.name
            m = re.match(r"(.+)_vs_(.+)\.weir\.fst$", p.name)
            pop1, pop2 = m.group(1), m.group(2)
            df = pd.read_csv(p, sep="\t", na_values=na)
            col = ("WEIGHTED_FST" if "WEIGHTED_FST" in df.columns
                   else "WEIR_AND_COCKERHAM_FST")
            v = (pd.to_numeric(df[col], errors="coerce").dropna().clip(lower=0)
                 if col in df.columns else pd.Series([], dtype=float))
            rows.append({"vclass": vclass, "pop1": pop1, "pop2": pop2,
                         "mean_fst": float(v.mean()) if len(v) else np.nan,
                         "n_sites": int(len(v))})
        out = pd.DataFrame(rows).sort_values(["vclass", "pop1", "pop2"])
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output.summary, sep="\t", index=False, float_format="%.6f")


# ---------------------------------------------------------------------------
# E5: core-vs-accessory structure concordance (Procrustes + Mantel)
# ---------------------------------------------------------------------------
rule region_pca_concordance:
    """Quantify how concordant the accessory-only and core-only sample
    configurations are (Procrustes disparity + Mantel r/p), from the already-
    built region PCAngsd covariances. Also emits per-morph accessory depth so
    the coverage confound is on record. Local, seconds."""
    input:
        cov_core=PCANGSD_OUTDIR / "augref_core/pcangsd.cov",
        samples_core=PCANGSD_OUTDIR / "augref_core/samples.txt",
        cov_accessory=PCANGSD_OUTDIR / "augref_accessory/pcangsd.cov",
        samples_accessory=PCANGSD_OUTDIR / "augref_accessory/samples.txt",
        region_qc=QC_OUTDIR / "region_qc_summary.tsv",
    output:
        stats=PCANGSD_OUTDIR / "augref_region_concordance.tsv",
    conda:
        "../envs/plotting.yaml"
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=4000,
        cpus=1,
    script:
        "../scripts/region_pca_concordance.py"


# ---------------------------------------------------------------------------
# E3: accessory missingness sweep (is the SFS/D skew a missingness artifact?)
# ---------------------------------------------------------------------------
# Cutoffs are the MAX allowed per-site missing-genotype fraction. Tightening
# should pull the accessory SFS back toward core if the intermediate-frequency
# skew is a missingness artifact rather than real shared variation.
MISSINGNESS_CUTOFFS = [0.5, 0.3, 0.1]


rule accessory_missingness_sweep:
    """Recompute the accessory folded SFS at a series of max-missingness cutoffs.
    Reads the RAW combined merged VCF (all SV_* sites, before LD-pruning and the
    MAF>=0.05 filter) -- that's where the intermediate-frequency skew lives; the
    pruned VCF has the rare bins already removed, so it can't test the artifact.
    The script subsets to SV_* contigs itself. Emits: cutoff, bin_low, bin_high,
    count, prop, n_sites_kept. Leaf rule, cheap."""
    input:
        vcf=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "bcftools/augref/combined/merged.vcf.gz.tbi",
    output:
        sweep=FST_OUTDIR / "afs/augref_accessory_missingness_sweep.tsv",
    conda:
        "../envs/fst_afs.yaml"
    params:
        cutoffs=MISSINGNESS_CUTOFFS,
        bins=AFS_BINS,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    script:
        "../scripts/accessory_missingness_sweep.py"


# ---------------------------------------------------------------------------
# E2 (OPT-IN, NOT in `rule all`): core-downsampling Tajima's D control
# ---------------------------------------------------------------------------
# This one genuinely re-genotypes thinned reads, so it is INTENTIONALLY left out
# of the default target set. Build it explicitly when you want the control, e.g.
#   snakemake results/pi/augref_core_downsampled/tajimas_d.tsv
# It subsamples core-region alignments to accessory-like depth (~1.7x) and
# recomputes Tajima's D; if core D also goes positive, the augref accessory D
# shift is proven to be depth, not biology.
#
# Left as a documented stub rather than a full rule because it needs the
# per-sample BAMs (heavy inputs) and a target-depth calc that should be tuned to
# your accessory mean. Wire it up only if E3 doesn't already settle the artifact
# question. See RESPONSE_TO_NOTES.md (E2) for the intended shape.


# ---------------------------------------------------------------------------
# Depth-aware PAV (presence/absence of accessory contigs) -- Bug 4b fix
# ---------------------------------------------------------------------------
# Replaces the old missingness-based presence (depth-confounded: cave accessory
# depth is ~2x surface) with a normalised-coverage call per contig per sample
# (contig depth / that sample's own core-genome baseline). See
# scripts/accessory_pav_coverage.py. PAV done right is the most promising
# augref-native analysis (Jeon's SV-PCA separated subspecies; a clean PAV matrix
# could do the same where SNPs can't).

PAV_CFG = config.get("accessory_pav", {})
PAV_MIN_BREADTH = PAV_CFG.get("min_breadth", 0.5)
PAV_MIN_NORM_COV = PAV_CFG.get("min_norm_cov", 0.3)


rule accessory_pav_per_sample:
    input:
        bam=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam",
        bai=ALIGN_OUTDIR / "{sample}/{sample}.augref.bam.bai",
        fai=f"{ALIGN_AUGREF}.fai",
    output:
        long_tsv=FST_OUTDIR / "pav/per_sample/{sample}.pav.tsv",
        matrix_tsv=temp(FST_OUTDIR / "pav/per_sample/{sample}.pav.col.tsv"),
    conda:
        "../envs/bcftools.yaml"
    params:
        sample=lambda wc: wc.sample,
        min_breadth=PAV_MIN_BREADTH,
        min_norm_cov=PAV_MIN_NORM_COV,
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=8000,
        cpus=1,
    script:
        "../scripts/accessory_pav_coverage.py"


rule accessory_pav_matrix:
    """Merge per-sample presence columns into one contig x sample 0/1 matrix, plus
    per-population presence frequencies and pairwise Hudson PAV FST WITH a
    label-permutation null (the negative control). Substrate for PAV PCA / morph
    classification. See scripts/accessory_pav_matrix.py."""
    input:
        cols=expand(FST_OUTDIR / "pav/per_sample/{sample}.pav.col.tsv", sample=SHORT_SAMPLES),
    output:
        matrix=FST_OUTDIR / "pav/accessory_pav_matrix.tsv",
        pop_freq=FST_OUTDIR / "pav/accessory_pav_pop_freq.tsv",
        fst=FST_OUTDIR / "pav/accessory_pav_fst.tsv",
    conda:
        "../envs/plotting.yaml"
    params:
        pop_samples=POP_SAMPLES,
        pop_pairs=POP_PAIR_TUPLES,
        n_perm=config.get("accessory_pav_perm", {}).get("n_perm", 1000),
        seed=config.get("accessory_pav_perm", {}).get("seed", 17),
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=8000,
        cpus=1,
    script:
        "../scripts/accessory_pav_matrix.py"
