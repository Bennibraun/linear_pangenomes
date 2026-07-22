import subprocess
from pathlib import Path

vcf        = snakemake.input.vcf
out_keys   = snakemake.output.keys()
out_afs    = snakemake.output.afs if "afs" in out_keys else None
out_region = snakemake.output.afs_by_region if "afs_by_region" in out_keys else None
bins       = snakemake.params.bins

freq_prefix = Path(out_afs if out_afs is not None else out_region).parent / "freq_tmp"
freq_prefix.parent.mkdir(parents=True, exist_ok=True)

subprocess.run(
    ["vcftools", "--gzvcf", vcf, "--freq2", "--max-alleles", "2", "--out", str(freq_prefix)],
    check=True,
)

freq_path = freq_prefix.with_suffix(".frq")


def bin_index(maf):
    last = len(bins) - 2
    for i in range(len(bins) - 1):
        # Include the right edge on the final bin so MAF == bins[-1] (e.g.
        # 0.5 with default bins) doesn't fall through unaccounted for.
        in_bin = bins[i] <= maf <= bins[i + 1] if i == last else bins[i] <= maf < bins[i + 1]
        if in_bin:
            return i
    return None


counts = [0] * (len(bins) - 1)
# Region-split counts are only meaningful for augref (SV_* accessory contigs);
# populated regardless, but only written out when the caller wants them.
region_counts = {"core": [0] * (len(bins) - 1), "accessory": [0] * (len(bins) - 1)}

for line in freq_path.read_text().strip().splitlines()[1:]:
    parts = line.split()
    if len(parts) < 6:
        continue
    chrom = parts[0]
    freqs = []
    for item in parts[4:]:
        try:
            freqs.append(float(item))
        except ValueError:
            continue
    if not freqs:
        continue
    maf = min(freqs)
    idx = bin_index(maf)
    if idx is None:
        continue
    counts[idx] += 1
    region = "accessory" if chrom.startswith("SV_") else "core"
    region_counts[region][idx] += 1

if out_afs is not None:
    out_path = Path(out_afs)
    out_path.write_text("bin_low\tbin_high\tcount\n")
    with out_path.open("a") as h:
        for i in range(len(counts)):
            h.write(str(bins[i]) + "\t" + str(bins[i + 1]) + "\t" + str(counts[i]) + "\n")

if out_region is not None:
    region_path = Path(out_region)
    region_path.parent.mkdir(parents=True, exist_ok=True)
    region_path.write_text("bin_low\tbin_high\tcount\tregion\n")
    with region_path.open("a") as h:
        for region in ("core", "accessory"):
            for i in range(len(bins) - 1):
                h.write(f"{bins[i]}\t{bins[i + 1]}\t{region_counts[region][i]}\t{region}\n")

freq_path.unlink(missing_ok=True)
freq_prefix.with_suffix(".log").unlink(missing_ok=True)
