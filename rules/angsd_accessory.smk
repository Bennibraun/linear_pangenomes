# ===========================================================================
# ANGSD genotype-likelihood arm for the augref ACCESSORY (SV_* + UNMAP_*)
# ===========================================================================
# The accessory contigs sit at ~0.8-1.9x depth (vs ~7x core). At that depth
# bcftools hard calls (`call -mv`) are single-read noise, which is what made the
# hard-call accessory Tajima's D / SFS / FST artifactual. ANGSD works from
# genotype LIKELIHOODS and integrates over genotype uncertainty per site, so it
# is the right estimator for low-coverage accessory sequence -- and it is what
# Jeon et al. use throughout.
#
# Ancestral = the augref itself and everything is FOLDED (-fold 1): accessory
# contigs are absent from conspec, so the inserted/novel allele is "reference"
# on its own contig and no polarised ancestral state is available.
#
# Scope: accessory contigs only. Core SNP stats stay on the existing bcftools
# path (depth there is fine and matches Jeon's own hard-call core usage).
# ===========================================================================

ANGSD_OUTDIR = FST_OUTDIR / "angsd_accessory"
AUGREF_FAI = f"{ALIGN_AUGREF}.fai"

# ANGSD site/quality filters, shared across the arm.
ANGSD_FILTERS = "-GL 1 -minMapQ 20 -minQ 20 -baq 1 -remove_bads 1 -uniqueOnly 1"

# Constrain popn/pop1/pop2 to actual population names so the accessory
# rules ({popn}/...) don't greedily match the E2 core paths (core_{cond}/{popn}/...)
# or the fst/ pair paths. Regex-escape names (they're plain, but be safe).
import re as _re
_POP_RE = "|".join(sorted((_re.escape(p) for p in POP_NAMES), key=len, reverse=True))


wildcard_constraints:
    popn=_POP_RE,
    pop1=_POP_RE,
    pop2=_POP_RE,
    cond="full|ds",


def _augref_bams_for_pop(pop):
    # POP_SHORT_SAMPLES (not POP_SAMPLES): only short-read samples have augref
    # BAMs. Using POP_SAMPLES would pull in the long-read assembly samples and
    # demand a short-read pipeline output for them (KeyError in downsample).
    return [str(ALIGN_OUTDIR / s / f"{s}.augref.bam") for s in POP_SHORT_SAMPLES[pop]]


# ---------------------------------------------------------------------------
# Region file: the accessory contigs (SV_* + UNMAP_*) as ANGSD -rf lines.
# ---------------------------------------------------------------------------
rule angsd_accessory_regions:
    input:
        fai=AUGREF_FAI,
    output:
        rf=ANGSD_OUTDIR / "accessory_regions.txt",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        Path(output.rf).parent.mkdir(parents=True, exist_ok=True)
        with open(input.fai) as fh, open(output.rf, "w") as out:
            for line in fh:
                name = line.split("\t", 1)[0]
                if name.startswith("SV_") or name.startswith("UNMAP_"):
                    out.write(name + ":\n")


# ---------------------------------------------------------------------------
# Per-population SAF (site allele frequency likelihoods) over the accessory.
# ---------------------------------------------------------------------------
rule angsd_saf_per_pop:
    input:
        bams=lambda wc: _augref_bams_for_pop(wc.popn),
        bais=lambda wc: [b + ".bai" for b in _augref_bams_for_pop(wc.popn)],
        ref=ALIGN_AUGREF,
        fai=AUGREF_FAI,
        rf=ANGSD_OUTDIR / "accessory_regions.txt",
    output:
        # angsd -doSaf writes a 3-file SAF set; declare all so Snakemake tracks
        # them together (realSFS/saf2theta need all three side by side).
        saf=ANGSD_OUTDIR / "{popn}/{popn}.saf.idx",
        saf_pos=ANGSD_OUTDIR / "{popn}/{popn}.saf.pos.gz",
        saf_gz=ANGSD_OUTDIR / "{popn}/{popn}.saf.gz",
        bamlist=temp(ANGSD_OUTDIR / "{popn}/{popn}.bamlist"),
    conda:
        "../envs/angsd.yaml"
    threads: 8
    resources:
        slurm_partition="long",
        runtime=720,
        mem_mb=24000,
        cpus=8,
    params:
        prefix=lambda wc: str(ANGSD_OUTDIR / wc.popn / wc.popn),
        filters=ANGSD_FILTERS,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        printf '%s\n' {input.bams} > {output.bamlist}
        angsd -bam {output.bamlist} \
            -ref {input.ref} -anc {input.ref} -fold 1 \
            -rf {input.rf} \
            {params.filters} \
            -doSaf 1 \
            -nThreads {threads} \
            -out {params.prefix}
        """


# ---------------------------------------------------------------------------
# Per-population folded SFS + thetas (pi, Tajima's D) from the SAF.
# ---------------------------------------------------------------------------
rule angsd_sfs_per_pop:
    input:
        saf=ANGSD_OUTDIR / "{popn}/{popn}.saf.idx",
        saf_pos=ANGSD_OUTDIR / "{popn}/{popn}.saf.pos.gz",
        saf_gz=ANGSD_OUTDIR / "{popn}/{popn}.saf.gz",
    output:
        sfs=ANGSD_OUTDIR / "{popn}/{popn}.sfs",
    conda:
        "../envs/angsd.yaml"
    threads: 8
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=16000,
        cpus=8,
    shell:
        r"""
        set -euo pipefail
        realSFS {input.saf} -fold 1 -P {threads} > {output.sfs}
        """


rule angsd_thetas_per_pop:
    input:
        saf=ANGSD_OUTDIR / "{popn}/{popn}.saf.idx",
        saf_pos=ANGSD_OUTDIR / "{popn}/{popn}.saf.pos.gz",
        saf_gz=ANGSD_OUTDIR / "{popn}/{popn}.saf.gz",
        sfs=ANGSD_OUTDIR / "{popn}/{popn}.sfs",
    output:
        pestPG=ANGSD_OUTDIR / "{popn}/{popn}.thetas.pestPG",
    conda:
        "../envs/angsd.yaml"
    threads: 4
    resources:
        slurm_partition="short",
        runtime=120,
        mem_mb=16000,
        cpus=4,
    params:
        prefix=lambda wc: str(ANGSD_OUTDIR / wc.popn / wc.popn),
    shell:
        r"""
        set -euo pipefail
        realSFS saf2theta {input.saf} -sfs {input.sfs} -fold 1 \
            -P {threads} -outnames {params.prefix}
        # Genome-wide (single-window) theta summary: pi, Watterson, Tajima's D.
        # thetaStat writes to <input>.pestPG (i.e. {prefix}.thetas.idx.pestPG);
        # rename to the declared output so the rule is self-consistent.
        thetaStat do_stat {params.prefix}.thetas.idx
        mv {params.prefix}.thetas.idx.pestPG {output.pestPG}
        """


# ---------------------------------------------------------------------------
# Per-population-pair accessory FST from 2D-SFS (Bhatia/Reynolds via realSFS).
# ---------------------------------------------------------------------------
rule angsd_fst_pair:
    input:
        saf1=ANGSD_OUTDIR / "{pop1}/{pop1}.saf.idx",
        saf2=ANGSD_OUTDIR / "{pop2}/{pop2}.saf.idx",
        # sibling SAF files realSFS reads by convention (kept explicit so the DAG
        # rebuilds them if missing, not just the .idx)
        saf1_pos=ANGSD_OUTDIR / "{pop1}/{pop1}.saf.pos.gz",
        saf1_gz=ANGSD_OUTDIR / "{pop1}/{pop1}.saf.gz",
        saf2_pos=ANGSD_OUTDIR / "{pop2}/{pop2}.saf.pos.gz",
        saf2_gz=ANGSD_OUTDIR / "{pop2}/{pop2}.saf.gz",
    output:
        fst=ANGSD_OUTDIR / "fst/{pop1}__{pop2}.fst.global",
    conda:
        "../envs/angsd.yaml"
    threads: 8
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=24000,
        cpus=8,
    params:
        prefix=lambda wc: str(ANGSD_OUTDIR / "fst" / f"{wc.pop1}__{wc.pop2}"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        realSFS {input.saf1} {input.saf2} -fold 1 -P {threads} \
            > {params.prefix}.2dsfs
        realSFS fst index {input.saf1} {input.saf2} \
            -sfs {params.prefix}.2dsfs -fold 1 \
            -fstout {params.prefix}
        # Global weighted + unweighted FST for the accessory.
        realSFS fst stats {params.prefix}.fst.idx > {output.fst}
        """


# ---------------------------------------------------------------------------
# Collect per-pop thetas and per-pair FST into two tidy tables.
# ---------------------------------------------------------------------------
rule angsd_accessory_summary:
    input:
        thetas=expand(ANGSD_OUTDIR / "{popn}/{popn}.thetas.pestPG", popn=POP_NAMES),
        fsts=expand(
            ANGSD_OUTDIR / "fst/{pair[0]}__{pair[1]}.fst.global",
            pair=POP_PAIR_TUPLES,
        ),
    output:
        thetas_tsv=ANGSD_OUTDIR / "accessory_thetas_summary.tsv",
        fst_tsv=ANGSD_OUTDIR / "accessory_fst_summary.tsv",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        import csv

        # thetaStat pestPG columns (tab-sep, one data row): the informative
        # fields are tW (Watterson theta), tP (pairwise theta = pi*nSites),
        # Tajima's D, and nSites. Header line starts with '#'.
        trows = []
        for f in input.thetas:
            pop = Path(f).stem.split(".")[0]
            with open(f) as fh:
                lines = [ln for ln in fh if ln.strip()]
            header = lines[0].lstrip("#").strip().split("\t")
            data = lines[1].strip().split("\t")
            rec = dict(zip(header, data))
            nsites = float(rec.get("nSites", 0) or 0)
            tP = float(rec.get("tP", 0) or 0)
            trows.append({
                "pop": pop,
                "tajimas_d": rec.get("Tajima", rec.get("Tajima's D", "")),
                "watterson_theta": rec.get("tW", ""),
                "pairwise_theta": rec.get("tP", ""),
                "pi_per_site": (tP / nsites) if nsites else 0.0,
                "n_sites": rec.get("nSites", ""),
            })
        with open(output.thetas_tsv, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=list(trows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(trows)

        # realSFS fst stats prints (one line, tab/space separated):
        #   FST.Unweight[nObs:N]:<u>  Fst.Weight:<w>
        # The label itself can contain a bracketed "[...:...]" so the numeric
        # value is always after the LAST colon; strip any bracketed part first.
        import re

        frows = []
        for f in input.fsts:
            name = Path(f).name.replace(".fst.global", "")
            pop_a, pop_b = name.split("__")
            txt = Path(f).read_text().strip()
            unweighted = weighted = ""
            for tok in txt.replace("\n", " ").split():
                if ":" not in tok:
                    continue
                label = re.sub(r"\[.*?\]", "", tok).lower()   # drop [nObs:N]
                value = tok.rsplit(":", 1)[1]                 # value after last colon
                if "unweight" in label:
                    unweighted = value
                elif "weight" in label:
                    weighted = value
            frows.append({
                "pop_a": pop_a, "pop_b": pop_b,
                "fst_unweighted": unweighted, "fst_weighted": weighted,
            })
        with open(output.fst_tsv, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=list(frows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(frows)


# ===========================================================================
# E2 control: core Tajima's D at full vs accessory-like (downsampled) depth
# ===========================================================================
# Proves the accessory Tajima's D shift is a depth artifact rather than biology.
# We run the SAME ANGSD SAF->SFS->thetaStat path on the CORE contigs at full
# depth and again with `-downSample frac` thinning core to accessory-like depth.
# If downsampled-core D tracks the accessory D (and full-depth-core D does not),
# the accessory signal was depth all along.

ANGSD_DS_FRAC = config.get("angsd_core_downsample", {}).get("frac", 0.15)


rule angsd_core_regions:
    input:
        fai=AUGREF_FAI,
    output:
        rf=ANGSD_OUTDIR / "core_regions.txt",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        Path(output.rf).parent.mkdir(parents=True, exist_ok=True)
        with open(input.fai) as fh, open(output.rf, "w") as out:
            for line in fh:
                name = line.split("\t", 1)[0]
                if not name.startswith(("SV_", "UNMAP_")):
                    out.write(name + ":\n")


rule angsd_core_saf_per_pop:
    """Core SAF per pop, at full depth (cond='full') or thinned to accessory-like
    depth (cond='ds'). The {cond} wildcard selects the -downSample argument."""
    input:
        bams=lambda wc: _augref_bams_for_pop(wc.popn),
        bais=lambda wc: [b + ".bai" for b in _augref_bams_for_pop(wc.popn)],
        ref=ALIGN_AUGREF,
        fai=AUGREF_FAI,
        rf=ANGSD_OUTDIR / "core_regions.txt",
    output:
        saf=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.idx",
        saf_pos=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.pos.gz",
        saf_gz=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.gz",
        bamlist=temp(ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.bamlist"),
    conda:
        "../envs/angsd.yaml"
    wildcard_constraints:
        cond="full|ds",
    threads: 8
    resources:
        slurm_partition="long",
        runtime=1440,
        mem_mb=24000,
        cpus=8,
    params:
        prefix=lambda wc: str(ANGSD_OUTDIR / f"core_{wc.cond}" / wc.popn / wc.popn),
        filters=ANGSD_FILTERS,
        downsample=lambda wc: f"-downSample {ANGSD_DS_FRAC}" if wc.cond == "ds" else "",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        printf '%s\n' {input.bams} > {output.bamlist}
        angsd -bam {output.bamlist} \
            -ref {input.ref} -anc {input.ref} -fold 1 \
            -rf {input.rf} \
            {params.filters} {params.downsample} \
            -doSaf 1 \
            -nThreads {threads} \
            -out {params.prefix}
        """


rule angsd_core_thetas_per_pop:
    input:
        saf=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.idx",
        saf_pos=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.pos.gz",
        saf_gz=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.saf.gz",
    output:
        pestPG=ANGSD_OUTDIR / "core_{cond}/{popn}/{popn}.thetas.pestPG",
    conda:
        "../envs/angsd.yaml"
    wildcard_constraints:
        cond="full|ds",
    threads: 8
    resources:
        slurm_partition="short",
        runtime=180,
        mem_mb=16000,
        cpus=8,
    params:
        prefix=lambda wc: str(ANGSD_OUTDIR / f"core_{wc.cond}" / wc.popn / wc.popn),
    shell:
        r"""
        set -euo pipefail
        realSFS {input.saf} -fold 1 -P {threads} > {params.prefix}.sfs
        realSFS saf2theta {input.saf} -sfs {params.prefix}.sfs -fold 1 \
            -P {threads} -outnames {params.prefix}
        thetaStat do_stat {params.prefix}.thetas.idx
        mv {params.prefix}.thetas.idx.pestPG {output.pestPG}
        """


rule angsd_core_downsample_summary:
    """One table: per-pop Tajima's D for accessory, core-full, core-downsampled.
    The decisive comparison for 'accessory D is depth, not biology'."""
    input:
        accessory=expand(ANGSD_OUTDIR / "{popn}/{popn}.thetas.pestPG", popn=POP_NAMES),
        core_full=expand(ANGSD_OUTDIR / "core_full/{popn}/{popn}.thetas.pestPG", popn=POP_NAMES),
        core_ds=expand(ANGSD_OUTDIR / "core_ds/{popn}/{popn}.thetas.pestPG", popn=POP_NAMES),
    output:
        tsv=ANGSD_OUTDIR / "core_downsample_tajima_control.tsv",
    resources:
        slurm_partition="short",
        runtime=20,
        mem_mb=2000,
        cpus=1,
    run:
        import csv

        def tajima_by_pop(files):
            out = {}
            for f in files:
                pop = Path(f).stem.split(".")[0]
                lines = [ln for ln in Path(f).read_text().splitlines() if ln.strip()]
                header = lines[0].lstrip("#").strip().split("\t")
                data = lines[1].strip().split("\t")
                rec = dict(zip(header, data))
                out[pop] = rec.get("Tajima", rec.get("Tajima's D", ""))
            return out

        acc = tajima_by_pop(input.accessory)
        cf = tajima_by_pop(input.core_full)
        cd = tajima_by_pop(input.core_ds)
        rows = []
        for pop in POP_NAMES:
            rows.append({
                "pop": pop,
                "tajimas_d_accessory": acc.get(pop, ""),
                "tajimas_d_core_full": cf.get(pop, ""),
                "tajimas_d_core_downsampled": cd.get(pop, ""),
            })
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
