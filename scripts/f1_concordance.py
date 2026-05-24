"""Compute F1 concordance between two reference SNP call sets.

Both VCFs must be in the same coordinate system (we restrict to the
intersection of contigs present in both headers). bcftools isec splits the
inputs into:
  0000.vcf — sites only in A
  0001.vcf — sites only in B
  0002.vcf — sites in both (A's records)
  0003.vcf — sites in both (B's records)

F1 is computed as the geometric balance of A→B precision/recall, treating
neither call set as ground truth:
  precision = shared / (shared + b_only)   (B's perspective: how much of B is in A)
  recall    = shared / (shared + a_only)   (A's perspective: how much of A is in B)
  F1        = 2·P·R / (P + R)

The script also reports the counts so the user can re-weight however they
prefer.
"""

import subprocess
import tempfile
from pathlib import Path

import pysam

vcf_a = snakemake.input.vcf_a
vcf_b = snakemake.input.vcf_b
out_path = Path(snakemake.output.f1)
ref_a = snakemake.params.ref_a
ref_b = snakemake.params.ref_b


def contig_set(vcf_path):
    with pysam.VariantFile(vcf_path) as v:
        return set(v.header.contigs)


# Restrict to the intersection of contigs present in both VCFs. Some refs
# have extra contigs (e.g. augref has SV_* contigs absent from conspec); we
# only count concordance on the shared contigs.
shared_contigs = sorted(contig_set(vcf_a) & contig_set(vcf_b))
if not shared_contigs:
    raise RuntimeError(
        f"No shared contigs between {ref_a} and {ref_b} — coordinate systems "
        "likely differ. Skip this pair or surject onto a common reference first."
    )


def count_records(vcf_path):
    n = 0
    with pysam.VariantFile(vcf_path) as v:
        for rec in v:
            if rec.contig in set(shared_contigs):
                n += 1
    return n


with tempfile.TemporaryDirectory() as tmpdir:
    regions = ",".join(shared_contigs)
    subprocess.run(
        [
            "bcftools", "isec", "-p", tmpdir,
            "-r", regions,
            "-Oz",
            vcf_a, vcf_b,
        ],
        check=True,
    )

    a_only = count_records(Path(tmpdir) / "0000.vcf.gz")
    b_only = count_records(Path(tmpdir) / "0001.vcf.gz")
    shared = count_records(Path(tmpdir) / "0002.vcf.gz")

precision = shared / (shared + b_only) if (shared + b_only) > 0 else 0.0
recall = shared / (shared + a_only) if (shared + a_only) > 0 else 0.0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as h:
    h.write("ref_a\tref_b\tshared\ta_only\tb_only\tprecision\trecall\tf1\n")
    h.write(
        f"{ref_a}\t{ref_b}\t{shared}\t{a_only}\t{b_only}\t"
        f"{precision:.6f}\t{recall:.6f}\t{f1:.6f}\n"
    )
