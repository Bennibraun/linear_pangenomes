"""Convert an LD-pruned VCF to a PLINK BED fileset for PCAngsd, remapping
CHROM/POS into a bounded number of placeholder chromosomes so plink2 never sees
more than a bounded number of distinct chromosome names.

Why this exists
---------------
augref's LD-pruned VCF carries one representative variant per accessory SV_*
contig (see rules/popgen.smk::ld_prune_accessory_small). There can be far more
than 65k such contigs, and plink2 --make-bed refuses once the distinct
nonstandard chromosome/contig count exceeds ~65k:

    Error: Too many distinct nonstandard chromosome/contig names.

The core/accessory_big split upstream keeps the *LD-pruning* plink2 pass under
the limit, but the final BED conversion for PCAngsd still saw every contig,
because it ran on the recombined whole-genome VCF. This script fixes that last
plink2 pass.

Approach
--------
plink2 is only used here for --make-bed (no window/LD operations), so CHROM is
just a grouping label -- pcangsd consumes the genotype matrix, and --make-bed
preserves input record order. So we:

  1. bcftools norm -m- to split multiallelics to biallelic (PLINK BED can't
     hold multiallelics), matching the previous behaviour.
  2. Capture the REAL (CHROM, POS, ID) for every record, in file order, from
     the normalized VCF. This is the exact order plink2 writes .bim rows in.
  3. Rewrite CHROM and POS to placeholder coordinates in CONTIGUOUS chunks of
     file order: record i gets CHROM="chunk{i // CHUNK_SIZE}" and
     POS=(i % CHUNK_SIZE) + 1. Each placeholder chromosome is thus one
     contiguous block with strictly increasing POS -- plink2 never sees
     interleaved contigs or unsorted positions. Distinct chrom count seen by
     plink2 is ceil(n_sites / CHUNK_SIZE), far under the ~65k cap. Real
     coordinates are recovered from the step-2 table, so overwriting POS here
     is harmless. REF/ALT/genotypes are untouched, so the genotype matrix
     PCAngsd consumes is identical.
  4. plink2 --make-bed on the remapped VCF.
  5. Emit the real-coords table (step 2) so downstream snp_coords.tsv keeps the
     true SV_* contig names, not the placeholders. Callers row-align this table
     against pcangsd's --sites mask.

The remapping is transparent to PCAngsd: covariance and selection statistics
depend on the genotype matrix and its row order, both of which are preserved.
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Records per placeholder chromosome. plink2's cap is ~65k DISTINCT chrom
# names, so this must be large enough that ceil(n_sites / CHUNK_SIZE) stays
# well under that: at 50k sites/chunk, ~3.2 billion sites would be needed to
# approach the cap. POS is rewritten to a within-chunk 1..CHUNK_SIZE counter,
# so CHUNK_SIZE must also stay under INT32_MAX (trivially true).
CHUNK_SIZE = 50000


def run(cmd, **kwargs):
    # plink2 and bcftools log to STDOUT. This module's own stdout is reserved
    # for the two result paths (see main), and callers capture it with $(...),
    # so send every child's stdout to stderr unless the caller redirects it
    # (e.g. to a results file). Otherwise a plink2 banner line would be mistaken
    # for the plink prefix downstream.
    kwargs.setdefault("stdout", sys.stderr)
    subprocess.run(cmd, check=True, **kwargs)


def build_plink_bed(in_vcf, tmpdir, memory_mb):
    """Normalize + chunk-remap `in_vcf`, run plink2 --make-bed, and return
    (plink_prefix, real_coords_path).

    real_coords is a headerless CHROM\\tPOS\\tID TSV, one row per plink .bim
    site in the SAME order, carrying the true (pre-remap) contig names.
    """
    tmpdir = Path(tmpdir)
    split_vcf = tmpdir / "split.vcf.gz"
    remapped_vcf = tmpdir / "remapped.vcf.gz"
    real_coords = tmpdir / "real_coords.tsv"
    plink_prefix = tmpdir / "plink"

    # 1. Split multiallelics to biallelic (PLINK BED can't hold multiallelics).
    run(["bcftools", "norm", "-m-", "-Oz", "-o", str(split_vcf), str(in_vcf)])
    run(["bcftools", "index", "-t", str(split_vcf)])

    # 2. Real (CHROM, POS, ID) in plink .bim order, from the normalized VCF.
    with real_coords.open("w") as fh:
        run(
            ["bcftools", "query", "-f", r"%CHROM\t%POS\t%ID\n", str(split_vcf)],
            stdout=fh,
        )
    n_sites = sum(1 for _ in real_coords.open())

    # 3. Stream the VCF as text and rewrite CHROM/POS to contiguous placeholder
    #    chunks (see module docstring). The ##contig header lines are rebuilt to
    #    match, since the originals reference the real (now-unused) contigs.
    _rewrite_to_chunks(split_vcf, remapped_vcf, n_sites)

    # 4. plink2 --make-bed on the remapped VCF.
    #    --const-fid 0: VCF has no family info; assign FID=0 to every sample.
    #    --allow-extra-chr: tolerate non-numeric contig names ("chunk123").
    run(
        [
            "plink2",
            "--vcf",
            str(remapped_vcf),
            "--make-bed",
            "--const-fid",
            "0",
            "--allow-extra-chr",
            "--memory",
            str(memory_mb),
            "--out",
            str(plink_prefix),
        ]
    )

    return str(plink_prefix), str(real_coords)


def _rewrite_to_chunks(in_vcf_gz, out_vcf_gz, n_sites):
    """Read bgzipped `in_vcf_gz`, write bgzipped `out_vcf_gz` with CHROM/POS of
    record i replaced by chunk{i // CHUNK_SIZE} / (i % CHUNK_SIZE)+1. All other
    columns pass through verbatim. Original ##contig header lines are dropped
    and replaced by one per placeholder chunk (emitted just before #CHROM)."""
    n_chunks = max(1, -(-n_sites // CHUNK_SIZE))  # ceil division

    reader = subprocess.Popen(
        ["bcftools", "view", str(in_vcf_gz)], stdout=subprocess.PIPE, text=True
    )
    with out_vcf_gz.open("wb") as out_fh:
        writer = subprocess.Popen(
            ["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh, text=True
        )
        try:
            i = 0
            for line in reader.stdout:
                if line.startswith("##contig="):
                    continue  # rebuilt below, just before the #CHROM line
                if line.startswith("#CHROM"):
                    for c in range(n_chunks):
                        writer.stdin.write(
                            f"##contig=<ID=chunk{c},length={CHUNK_SIZE}>\n"
                        )
                    writer.stdin.write(line)
                    continue
                if line.startswith("#"):
                    writer.stdin.write(line)
                    continue
                # CHROM<TAB>POS<TAB>rest(ID onward). maxsplit=2 keeps genotype
                # columns (which contain tabs) untouched in fields[2].
                fields = line.rstrip("\n").split("\t", 2)
                new_chrom = f"chunk{i // CHUNK_SIZE}"
                new_pos = (i % CHUNK_SIZE) + 1
                writer.stdin.write(f"{new_chrom}\t{new_pos}\t{fields[2]}\n")
                i += 1
            writer.stdin.close()
        finally:
            reader.stdout.close()
        w_rc = writer.wait()
    r_rc = reader.wait()
    if r_rc != 0:
        raise RuntimeError(f"bcftools view failed (rc={r_rc}) reading {in_vcf_gz}")
    if w_rc != 0:
        raise RuntimeError(f"bgzip failed (rc={w_rc}) writing {out_vcf_gz}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="input (LD-pruned) VCF")
    p.add_argument("--tmpdir", required=True, help="scratch dir for intermediates")
    p.add_argument("--memory-mb", type=int, required=True, help="plink2 --memory")
    args = p.parse_args()

    plink_prefix, real_coords = build_plink_bed(
        args.vcf, args.tmpdir, args.memory_mb
    )
    # Print the two paths the calling shell needs, one per line.
    sys.stdout.write(f"{plink_prefix}\n{real_coords}\n")


if __name__ == "__main__":
    main()
