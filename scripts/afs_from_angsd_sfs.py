"""Reshape ANGSD realSFS folded-SFS output into the canonical
bin_low / bin_high / count TSV used by plot_afs.py and the
variant-count summary.

realSFS -fold 1 emits a whitespace-separated row of (N+1) expected site
counts (where N is the number of diploid samples × 2 = chromosomes / 2),
indexed by minor allele count from 0 to N. We collapse those expected
counts into the same minor-allele-frequency bins used by the linear /
graph variant pipelines so the AFS plots are directly comparable.

Values are expected (fractional) site counts from the ML estimator, not
integer counts — we round when emitting since the downstream consumers
treat the column as a count.
"""

from pathlib import Path

sfs_path = snakemake.input.sfs
out_path = Path(snakemake.output.afs)
bins = snakemake.params.bins
n_samples = snakemake.params.n_samples

# realSFS folded output is a single whitespace-separated row (sometimes
# multiple if the SFS was per-region). We sum any multi-row outputs into
# one cohort SFS — same semantics as the unfolded → folded reduction.
rows = []
for line in Path(sfs_path).read_text().splitlines():
    parts = line.split()
    if not parts:
        continue
    try:
        rows.append([float(x) for x in parts])
    except ValueError:
        continue

if not rows:
    raise RuntimeError(f"Empty or malformed SFS: {sfs_path}")

# Sum across rows (handles per-region SFS outputs); use the first row's
# length as canonical so a stray short row doesn't truncate.
n_bins = max(len(r) for r in rows)
sfs = [0.0] * n_bins
for r in rows:
    for i, v in enumerate(r):
        sfs[i] += v

# Folded SFS is indexed by minor allele count (MAC) from 0 to n_samples
# (diploid sample count). Convert each MAC index to MAF and bin.
# MAC=0 entries are monomorphic sites — exclude from the AFS plot.
counts = [0.0] * (len(bins) - 1)
last_bin = len(bins) - 2

for mac, val in enumerate(sfs):
    if mac == 0:
        continue  # monomorphic — not a real variant
    maf = mac / (2 * n_samples) if n_samples > 0 else 0.0
    for i in range(len(bins) - 1):
        in_bin = (
            bins[i] <= maf <= bins[i + 1]
            if i == last_bin
            else bins[i] <= maf < bins[i + 1]
        )
        if in_bin:
            counts[i] += val
            break

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as h:
    h.write("bin_low\tbin_high\tcount\n")
    for i, c in enumerate(counts):
        # Round expected counts to integers for the canonical schema.
        h.write(f"{bins[i]}\t{bins[i + 1]}\t{int(round(c))}\n")
