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

        echo "{params.pop1_samples}" > /tmp/pop1_samples.txt
        echo "{params.pop2_samples}" > /tmp/pop2_samples.txt

        vcftools --gzvcf {input.vcf} \
          --weir-fst-pop /tmp/pop1_samples.txt \
          --weir-fst-pop /tmp/pop2_samples.txt \
          {params.window_args} \
          --out $prefix > /dev/null

        rm /tmp/pop1_samples.txt /tmp/pop2_samples.txt
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

        autosome_map = params.autosomes
        autosomes = autosome_map.get(wildcards.ref, [])
        if autosomes:
            genome_len = sum(lengths.get(c, 0) for c in autosomes)
        else:
            genome_len = sum(lengths.values())

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
