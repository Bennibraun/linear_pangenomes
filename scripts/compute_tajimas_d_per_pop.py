"""Per-population Tajima's D in non-overlapping windows (vcftools --TajimaD --keep).

Runs vcftools once per population on the same merged VCF, restricting to that
population's samples, then concatenates into one long-format table with a POP
column:

    CHROM  BIN_START  N_SNPS  TajimaD  POP

Populations with fewer than 2 samples present in the VCF are skipped (Tajima's D
needs multiple sequences). vcftools writes "nan" for windows with too few SNPs;
those rows are kept as-is and filtered at plot time.
"""
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

vcf          = snakemake.input.vcf
out_path     = Path(snakemake.output[0])
window_size  = int(snakemake.params.window_size)
pop_samples  = dict(snakemake.params.pop_samples)

out_path.parent.mkdir(parents=True, exist_ok=True)

present = set(
    subprocess.run(
        ["bcftools", "query", "-l", vcf],
        capture_output=True, text=True, check=True,
    ).stdout.split()
)

frames = []
with tempfile.TemporaryDirectory(dir=str(out_path.parent)) as tmp:
    tmp = Path(tmp)
    for pop, samples in pop_samples.items():
        keep = [s for s in samples if s in present]
        if len(keep) < 2:
            continue

        keep_file = tmp / f"{pop}.keep"
        keep_file.write_text("\n".join(keep) + "\n")
        prefix = tmp / pop

        subprocess.run(
            ["vcftools", "--gzvcf", vcf,
             "--keep", str(keep_file),
             "--TajimaD", str(window_size),
             "--out", str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

        # vcftools appends ".Tajima.D" to --out; concatenate (NOT with_suffix,
        # which mangles pop names containing a ".").
        td_file = Path(str(prefix) + ".Tajima.D")
        if not td_file.exists():
            continue
        df = pd.read_csv(td_file, sep="\t")
        if df.empty:
            continue
        df["POP"] = pop
        frames.append(df)

if frames:
    out = pd.concat(frames, ignore_index=True)
else:
    out = pd.DataFrame(columns=["CHROM", "BIN_START", "N_SNPS", "TajimaD", "POP"])

out.to_csv(out_path, sep="\t", index=False)
