"""Per-population windowed nucleotide diversity (vcftools --window-pi --keep).

Runs vcftools once per population on the same merged VCF, restricting to that
population's samples, then concatenates the per-population outputs into one
long-format table with a POP column:

    CHROM  BIN_START  BIN_END  N_VARIANTS  PI  POP

Populations with fewer than 2 samples present in the VCF are skipped (π is
undefined for a single haplotype set). A population is also skipped if vcftools
produces no windows (e.g. too few variants).
"""
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

vcf          = snakemake.input.vcf
out_path     = Path(snakemake.output[0])
window_size  = int(snakemake.params.window_size)
window_step  = int(snakemake.params.window_step)
pop_samples  = dict(snakemake.params.pop_samples)

out_path.parent.mkdir(parents=True, exist_ok=True)

# Which samples are actually in this VCF (a population may not be fully present).
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
             "--window-pi", str(window_size),
             "--window-pi-step", str(window_step),
             "--out", str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

        # vcftools appends ".windowed.pi" to the --out value; build the path by
        # concatenation (NOT Path.with_suffix, which would mangle pop names that
        # contain a ".").
        pi_file = Path(str(prefix) + ".windowed.pi")
        if not pi_file.exists():
            continue
        df = pd.read_csv(pi_file, sep="\t")
        if df.empty:
            continue
        df["POP"] = pop
        frames.append(df)

if frames:
    out = pd.concat(frames, ignore_index=True)
else:
    # Preserve the schema so downstream readers don't crash on an empty run.
    out = pd.DataFrame(
        columns=["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI", "POP"]
    )

out.to_csv(out_path, sep="\t", index=False)
