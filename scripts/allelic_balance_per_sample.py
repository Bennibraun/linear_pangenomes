import statistics
import pysam
from pathlib import Path

vcf_path    = snakemake.input.vcf
out_summary = snakemake.output.summary
out_raw     = snakemake.output.raw
ref_name    = snakemake.params.ref
bins        = snakemake.params.bins

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

Path(out_summary).parent.mkdir(parents=True, exist_ok=True)
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
