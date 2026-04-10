# ============================================================================
# Population Genetics Rules
# ============================================================================
# Variables (VC_OUTDIR, FST_OUTDIR, PI_OUTDIR, AB_OUTDIR, ROH_OUTDIR, etc.)
# and helper functions (_gatk_ref_fasta, _gatk_ref_fai) are defined in the
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
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
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
    run:
        from pathlib import Path
        import subprocess

        out_path = Path(output.afs)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        freq_prefix = Path(output.afs).parent / "freq_tmp"
        subprocess.run(
            f"vcftools --gzvcf {input.vcf} --freq2 --max-alleles 2 --out {freq_prefix}",
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
            for item in parts[4:]:  # columns 4+ are raw frequencies with --freq2
                try:
                    freqs.append(float(item))
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


# ---------------------------------------------------------------------------
# FST between population pairs
# ---------------------------------------------------------------------------
rule fst_per_ref_pair:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
    output:
        FST_OUTDIR / "fst/{ref}/{pop1}_vs_{pop2}.weir.fst"
    conda:
        "../envs/fst_afs.yaml"
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

        pop1_tmp=$(mktemp)
        pop2_tmp=$(mktemp)
        trap "rm -f $pop1_tmp $pop2_tmp" EXIT

        printf '%s\n' {params.pop1_samples} > "$pop1_tmp"
        printf '%s\n' {params.pop2_samples} > "$pop2_tmp"

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
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi"
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
# Allelic balance (ref bias) per sample
# ---------------------------------------------------------------------------
rule allelic_balance_per_sample:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi"
    output:
        summary=AB_OUTDIR / "{ref}/{sample}.allelic_balance.tsv",
        raw=AB_OUTDIR / "{ref}/{sample}.allelic_balance.raw.tsv"
    conda:
        "../envs/allelic_balance.yaml"
    params:
        bins=lambda wildcards: " ".join(str(b) for b in AB_BINS),
        ref=lambda wildcards: wildcards.ref,
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.summary})
        python - << 'PYEOF'
import statistics
import pysam

vcf_path   = "{input.vcf}"
out_summary = "{output.summary}"
out_raw     = "{output.raw}"
ref_name    = "{params.ref}"
bins        = [float(x) for x in "{params.bins}".split()]

vcf     = pysam.VariantFile(vcf_path)
samples = list(vcf.header.samples)
if len(samples) != 1:
    raise ValueError("Expected 1 sample in " + vcf_path + ", found " + str(len(samples)))
sample = samples[0]

counts   = [0] * (len(bins) - 1)
ratios   = []
raw_data = []

for rec in vcf.fetch():
    if "AD" not in rec.format or "GT" not in rec.format:
        continue
    gt = rec.samples[sample].get("GT")
    if gt is None or len(gt) < 2 or gt[0] == gt[1]:
        continue
    ad = rec.samples[sample].get("AD")
    if ad is None or len(ad) < 2:
        continue
    ref_depth, alt_depth = ad[0], ad[1]
    if ref_depth is None or alt_depth is None:
        continue
    total = ref_depth + alt_depth
    if total == 0:
        continue
    ratio = ref_depth / total
    ratios.append(ratio)
    raw_data.append((str(rec.contig) + ":" + str(rec.pos), ref_depth, alt_depth, ratio))
    for i in range(len(bins) - 1):
        if bins[i] <= ratio < bins[i + 1]:
            counts[i] += 1
            break

mean_ratio   = statistics.mean(ratios)   if ratios else 0.0
median_ratio = statistics.median(ratios) if ratios else 0.0
total_hets   = len(ratios)

with open(out_summary, "w") as h:
    h.write("sample\tref\ttotal_hets\tmean_ref_ratio\tmedian_ref_ratio\n")
    h.write(sample + "\t" + ref_name + "\t" + str(total_hets) + "\t"
            + format(mean_ratio, ".6f") + "\t" + format(median_ratio, ".6f") + "\n")
    h.write("bin_low\tbin_high\tcount\n")
    for i in range(len(counts)):
        h.write(str(bins[i]) + "\t" + str(bins[i + 1]) + "\t" + str(counts[i]) + "\n")

with open(out_raw, "w") as h:
    h.write("site\tref_depth\talt_depth\tref_ratio\n")
    for site, rd, ad, ratio in raw_data:
        h.write(site + "\t" + str(rd) + "\t" + str(ad) + "\t" + format(ratio, ".6f") + "\n")
PYEOF
        """


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
        fasta=lambda wildcards: _gatk_ref_fasta(wildcards),
        fai=lambda wildcards: _gatk_ref_fai(wildcards)
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


rule roh_per_sample:
    input:
        vcf=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/merged/{sample}.vcf.gz.tbi",
        lengths=lambda wildcards: ROH_OUTDIR / "lengths" / f"{wildcards.ref}.lengths.tsv"
    output:
        roh=ROH_OUTDIR / "{ref}/{sample}.roh.tsv",
        froh=ROH_OUTDIR / "{ref}/{sample}.f_roh.tsv"
    conda:
        "../envs/roh.yaml"
    params:
        canonical_chroms=ROH_CANONICAL_CHROMS,
        min_roh_length=ROH_MIN_LENGTH,
        bcftools_args=ROH_BCFTOOLS_ARGS
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=4000,
        cpus=1
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
        vcf=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz",
        tbi=VC_OUTDIR / "gatk/{ref}/combined/merged.vcf.gz.tbi",
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

        # Step 1: compute LD prune list
        plink \
            --vcf {input.vcf} \
            --double-id --allow-extra-chr \
            --set-missing-var-ids @:#$$1$$2 \
            --snps-only just-acgt \
            --maf {params.maf} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out "$outdir/prune_tmp"

        # Step 2: extract pruned sites → VCF
        plink \
            --vcf {input.vcf} \
            --double-id --allow-extra-chr \
            --set-missing-var-ids @:#$$1$$2 \
            --snps-only just-acgt \
            --maf {params.maf} \
            --extract "$outdir/prune_tmp.prune.in" \
            --recode vcf \
            --out "$outdir/pruned_tmp"

        bgzip -f "$outdir/pruned_tmp.vcf"
        mv "$outdir/pruned_tmp.vcf.gz" {output.vcf}
        tabix -p vcf {output.vcf}

        rm -f "$outdir/prune_tmp".* "$outdir/pruned_tmp.log" "$outdir/pruned_tmp.nosex"
        """


rule pcangsd:
    """Run PCAngsd on LD-pruned VCF to get covariance matrix and selection scores."""
    input:
        vcf=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz",
        tbi=PCANGSD_OUTDIR / "{ref}/merged.ldpruned.vcf.gz.tbi",
    output:
        cov=PCANGSD_OUTDIR / "{ref}/pcangsd.cov",
        selection=PCANGSD_OUTDIR / "{ref}/pcangsd.selection.npy",
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

        # Sample order in VCF = row order in covariance matrix
        bcftools query -l {input.vcf} > {output.samples}

        # SNP coordinates in same order PCAngsd sees them
        bcftools query -f '%CHROM\t%POS\t%ID\n' {input.vcf} > {output.snp_coords}

        pcangsd \
            --vcf {input.vcf} \
            --threads {threads} \
            --n_eig {params.n_pcs} \
            --selection \
            --out {params.prefix}
        """


rule fst_outliers:
    """
    Flag per-SNP selection outliers from PCAngsd chi-squared scores with
    Benjamini-Hochberg FDR correction. Outputs a TSV of all loci with
    coordinates, chi-squared score, p-value, q-value, and outlier flag.
    """
    input:
        selection=PCANGSD_OUTDIR / "{ref}/pcangsd.selection.npy",
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
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.outliers})
        python - << 'PYEOF'
import numpy as np
import pandas as pd
from scipy import stats

scores = np.load("{input.selection}")
coords = pd.read_csv("{input.snp_coords}", sep="\t", names=["chrom", "pos", "id"])
fdr    = {params.fdr}

# PCAngsd selection scores ~ chi-squared(df=1) under the null
pvals = stats.chi2.sf(scores, df=1)

n        = len(pvals)
ranks    = np.argsort(pvals)
sorted_p = pvals[ranks]
below    = sorted_p <= (np.arange(1, n + 1) / n) * fdr
cutoff   = sorted_p[below].max() if below.any() else 0.0

qvals_sorted = np.minimum(1.0, sorted_p * n / np.arange(1, n + 1))
qvals        = np.empty(n)
qvals[ranks] = qvals_sorted

coords["chi2"]    = scores
coords["pval"]    = pvals
coords["qval"]    = qvals
coords["outlier"] = pvals <= cutoff

coords.sort_values(["chrom", "pos"]).to_csv("{output.outliers}", sep="\t", index=False)
PYEOF
        """
