import subprocess
import tempfile
from pathlib import Path

bcf            = snakemake.input.bcf
freqs          = snakemake.input.freqs
lengths_file   = snakemake.input.lengths
out_roh        = Path(snakemake.output.roh)
out_froh       = Path(snakemake.output.froh)

min_roh_length = snakemake.params.min_roh_length
bcftools_args  = snakemake.params.bcftools_args
sample         = snakemake.wildcards.sample
ref            = snakemake.wildcards.ref

out_roh.parent.mkdir(parents=True, exist_ok=True)

lengths = {}
for line in Path(lengths_file).read_text().splitlines():
    parts = line.split()
    if len(parts) < 2:
        continue
    try:
        lengths[parts[0]] = int(parts[1])
    except ValueError:
        continue

# Canonical chromosomes for this ref = all contigs except augref-style SV
# alt contigs. Per-ref because each reference has its own contig set.
canonical = [c for c in lengths.keys() if not str(c).startswith("SV_")]

regions_arg = ""
if canonical:
    regions_arg = "--regions " + ",".join(canonical)

with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".roh") as tmp:
    roh_tmp = tmp.name

# Paper-style ROH: PL-based HMM (-G30) using population AF from ANGSD's
# tabixed freq table. The BCF holds genotype likelihoods, and --AF-file
# overrides any in-BCF AF with the cohort frequencies we extracted.
subprocess.run(
    f"bcftools roh {bcftools_args} --AF-file {freqs} -G30 "
    f"-s {sample} {regions_arg} -o {roh_tmp} {bcf}",
    shell=True, check=True,
)

if canonical:
    genome_len = sum(lengths.get(c, 0) for c in canonical)
else:
    genome_len = sum(lengths.values())

# Parse bcftools roh output — only RG (region) lines.
# RG format: RG sample chrom start end length_bp n_markers quality
roh_rows = []
roh_total = 0
canonical_set = set(canonical)
for line in Path(roh_tmp).read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    parts = line.split()
    if parts[0] != "RG" or len(parts) < 8:
        continue
    chrom = parts[2]
    try:
        start = int(parts[3])
        end = int(parts[4])
    except ValueError:
        continue
    length = end - start + 1
    if canonical_set and chrom not in canonical_set:
        continue
    if length < min_roh_length:
        continue
    roh_total += length
    roh_rows.append((chrom, start, end, length))

out_roh.write_text("chrom\tstart\tend\tlength\n")
with out_roh.open("a") as fh:
    for chrom, start, end, length in roh_rows:
        fh.write(f"{chrom}\t{start}\t{end}\t{length}\n")

f_roh = (roh_total / genome_len) if genome_len > 0 else 0.0
out_froh.write_text("sample\tref\troh_bp\tgenome_bp\tf_roh\n")
with out_froh.open("a") as fh:
    fh.write(f"{sample}\t{ref}\t{roh_total}\t{genome_len}\t{f_roh:.6f}\n")

Path(roh_tmp).unlink()
