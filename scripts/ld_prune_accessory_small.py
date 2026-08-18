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
        ["bcftools", "annotate", "--set-id", "+%CHROM:%POS_%REF_%FIRST_ALT", "-Oz"],
        stdin=maf_filtered.stdout,
        stdout=out_fh,
    )
    maf_filtered.stdout.close()
    id_set.communicate()
    if id_set.returncode != 0:
        raise RuntimeError("bcftools annotate --set-id failed")

run(["bcftools", "index", "-t", str(filtered_vcf)])

# --- Step 5: collapse to 1 representative variant per contig -------------
# For each contig, keep the site ID with the most non-missing genotypes.
query = subprocess.run(
    ["bcftools", "query", "-f", "%CHROM\t%POS\t%ID\t[%GT,]\n", str(filtered_vcf)],
    check=True, capture_output=True, text=True,
)

best_per_contig = {}  # chrom -> (n_non_missing, id)
for line in query.stdout.splitlines():
    chrom, pos, var_id, gt_field = line.rstrip("\n").split("\t")
    gts = gt_field.rstrip(",").split(",")
    n_non_missing = sum(1 for gt in gts if gt != "./.")
    current_best = best_per_contig.get(chrom)
    if current_best is None or n_non_missing > current_best[0]:
        best_per_contig[chrom] = (n_non_missing, var_id)

keep_ids = [var_id for _, var_id in best_per_contig.values()]

if keep_ids:
    keep_ids_path = tmp_dir / "accessory_small_keep_ids.txt"
    keep_ids_path.write_text("\n".join(keep_ids) + "\n")
    run([
        "bcftools", "view", "-i", f"ID=@{keep_ids_path}",
        str(filtered_vcf), "-Oz", "-o", str(out_vcf),
    ])
    keep_ids_path.unlink()
else:
    # No accessory_small sites survived filtering -- emit a header-only VCF.
    run(["bcftools", "view", "-h", str(filtered_vcf), "-Oz", "-o", str(out_vcf)])

filtered_vcf.unlink()
Path(str(filtered_vcf) + ".tbi").unlink(missing_ok=True)
