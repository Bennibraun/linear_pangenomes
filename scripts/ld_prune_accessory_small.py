"""Filter and collapse "accessory_small" SV_* sites for ld_prune_vcf.

accessory_small contigs (SV_* contigs shorter than pcangsd.accessory_ld_min_len,
see rules/popgen.smk::ld_prune_vcf) are too short for real LD-pruning to be
worth it: all variants on one short flanked-SV fragment are essentially
guaranteed to be in near-total LD with each other (same fragment, no realistic
within-contig recombination). Longer accessory contigs (accessory_big) go
through actual plink2 LD-pruning instead, folded in with core in the calling
rule's shell block.

This script:
  1. Restricts to biallelic SNPs (bcftools view -m2 -M2 -v snps), matching the
     core branch's plink2 --snps-only just-acgt filter.
  2. Masks (sets to missing) individual genotypes below accessory_min_dp /
     accessory_min_gq. SV_* contigs are short and often low/no coverage per
     sample, so bcftools call emits single-read genotypes at many sites
     (DP=1, GQ<20) alongside ./. for other samples at the same position --
     noise, not real presence/absence signal. Masking happens BEFORE the MAF
     filter so a site's frequency is computed from confidently-genotyped
     samples only.
  3. Applies the MAF filter (matching the core branch's plink2 --maf).
  4. Sets missing variant IDs (matching plink2 --set-missing-var-ids).
  5. Collapses to 1 representative variant per contig: keeps the site with
     the most non-missing genotypes (post-masking) as a quality tiebreak.
     This is NOT LD-pruning -- it's a stand-in that stops a multi-SNP contig
     from being overweighted relative to a true single-variant contig in
     PCA/selection-scan input.

This pipeline is bcftools-only (no plink2), so sample names are never
touched -- unlike the core/accessory_big branch, no reheader step is needed
here before recombining in ld_prune_vcf.
"""

import subprocess
from pathlib import Path

import pysam

in_vcf = snakemake.input.vcf
out_vcf = Path(snakemake.output.vcf)
min_dp = snakemake.params.min_dp
min_gq = snakemake.params.min_gq
maf = snakemake.params.maf
tmp_dir = out_vcf.parent

filtered_vcf = tmp_dir / "accessory_small_filtered.vcf.gz"


def run(cmd, **kwargs):
    subprocess.run(cmd, check=True, **kwargs)


# --- Steps 1-4: biallelic SNPs -> depth/GQ genotype masking -> MAF -> set-ID
biallelic = subprocess.Popen(
    ["bcftools", "view", "-m2", "-M2", "-v", "snps", in_vcf],
    stdout=subprocess.PIPE,
)
masked = subprocess.Popen(
    ["bcftools", "filter", "-e", f"FMT/DP<{min_dp} | FMT/GQ<{min_gq}", "--set-GTs", "."],
    stdin=biallelic.stdout,
    stdout=subprocess.PIPE,
)
biallelic.stdout.close()
maf_filtered = subprocess.Popen(
    ["bcftools", "view", "-i", f"MAF>={maf}"],
    stdin=masked.stdout,
    stdout=subprocess.PIPE,
)
masked.stdout.close()
with filtered_vcf.open("wb") as out_fh:
    id_set = subprocess.Popen(
        # Underscores immediately after a %TAG are parsed as part of the tag
        # name (bcftools tried to look up INFO/REF_ here) unless escaped.
        ["bcftools", "annotate", "--set-id", r"+%CHROM:%POS\_%REF\_%FIRST_ALT", "-Oz"],
        stdin=maf_filtered.stdout,
        stdout=out_fh,
    )
    maf_filtered.stdout.close()
    id_set.communicate()
    if id_set.returncode != 0:
        raise RuntimeError("bcftools annotate --set-id failed")

run(["bcftools", "index", "-t", str(filtered_vcf)])

# --- Step 5: collapse to 1 representative variant per contig -------------
# For each contig, keep the record with the most non-missing genotypes.
# Done with pysam (not `bcftools view -i 'ID=@file'`) -- that filter
# expression failed outright on this bcftools build ("Could not read:
# @file"), so record selection is done directly here instead of relying on
# bcftools's -i list-membership syntax.
def _is_missing_gt(sample):
    """True if a sample's genotype is fully missing (./. or equivalent).
    pysam's allele_indices is either None (whole-genotype-missing) or a
    tuple that may itself contain None per-allele -- handle both rather
    than assume one shape."""
    idx = sample.allele_indices
    if idx is None:
        return True
    return all(a is None for a in idx)


with pysam.VariantFile(str(filtered_vcf)) as vf:
    best_per_contig = {}  # chrom -> (n_non_missing, record.id)
    for rec in vf:
        n_non_missing = sum(
            1 for sample in rec.samples.values() if not _is_missing_gt(sample)
        )
        current_best = best_per_contig.get(rec.chrom)
        if current_best is None or n_non_missing > current_best[0]:
            best_per_contig[rec.chrom] = (n_non_missing, rec.id)

keep_ids = {var_id for _, var_id in best_per_contig.values()}

tmp_out_vcf = tmp_dir / "accessory_small_collapsed.vcf.gz"
with pysam.VariantFile(str(filtered_vcf)) as vf_in:
    with pysam.VariantFile(str(tmp_out_vcf), "wz", header=vf_in.header) as vf_out:
        for rec in vf_in:
            if rec.id in keep_ids:
                vf_out.write(rec)

tmp_out_vcf.rename(out_vcf)

# ld_prune_vcf's bcftools concat needs this indexed.
run(["bcftools", "index", "-t", str(out_vcf)])

filtered_vcf.unlink()
Path(str(filtered_vcf) + ".tbi").unlink(missing_ok=True)
