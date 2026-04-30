import subprocess
from pathlib import Path

vcf      = snakemake.input.vcf
out_afs  = snakemake.output.afs
bins     = snakemake.params.bins

out_path    = Path(out_afs)
out_path.parent.mkdir(parents=True, exist_ok=True)
freq_prefix = out_path.parent / "freq_tmp"

subprocess.run(
    ["vcftools", "--gzvcf", vcf, "--freq2", "--max-alleles", "2", "--out", str(freq_prefix)],
    check=True,
)

freq_path = freq_prefix.with_suffix(".frq")
counts    = [0] * (len(bins) - 1)

for line in freq_path.read_text().strip().splitlines()[1:]:
    parts = line.split()
    if len(parts) < 6:
        continue
    freqs = []
    for item in parts[4:]:
        try:
            freqs.append(float(item))
        except ValueError:
            continue
    if not freqs:
        continue
    maf = min(freqs)
    last = len(bins) - 2
    for i in range(len(bins) - 1):
        # Include the right edge on the final bin so MAF == bins[-1] (e.g.
        # 0.5 with default bins) doesn't fall through unaccounted for.
        in_bin = bins[i] <= maf <= bins[i + 1] if i == last else bins[i] <= maf < bins[i + 1]
        if in_bin:
            counts[i] += 1
            break

out_path.write_text("bin_low\tbin_high\tcount\n")
with out_path.open("a") as h:
    for i in range(len(counts)):
        h.write(str(bins[i]) + "\t" + str(bins[i + 1]) + "\t" + str(counts[i]) + "\n")

freq_path.unlink(missing_ok=True)
freq_prefix.with_suffix(".log").unlink(missing_ok=True)
