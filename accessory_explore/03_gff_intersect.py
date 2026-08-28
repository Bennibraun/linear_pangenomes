#!/usr/bin/env python3
"""Step 3: does each accessory contig's INSERTION LOCUS fall in/near a conspec
gene, and does that gene match the contig's DIAMOND protein hit?

The DIAMOND hit says the inserted SEQUENCE resembles gene X. This says WHERE on
the conspec backbone the insertion landed (parsed from the contig name
SV_<chrom>_<pos>_...). A contig is most credible when both agree: an X-like
insertion sitting in/near the X locus.

Conspec = GCF_023375975.1_AstMex3_surface (chroms NC_064408-064431).

Usage:
  03_gff_intersect.py OVERLAY_TSV GFF_GZ [OUTDIR]
GFF from:
  ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/
    GCF_023375975.1_AstMex3_surface_genomic.gff.gz
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

FLANK = 5000   # call a locus "near" a gene if within this many bp


def load_genes(gff_gz):
    """Return DataFrame chrom,start,end,gene for protein-coding genes."""
    rows = []
    op = gzip.open if str(gff_gz).endswith(".gz") else open
    with op(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"gene=([^;]+)", f[8])
            rows.append((f[0], int(f[3]), int(f[4]), m.group(1) if m else ""))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def parse_locus(contig):
    m = re.match(r"SV_(NC_[\d.]+)_(\d+)_", contig)
    return (m.group(1), int(m.group(2))) if m else (None, None)


def gene_at(genes, chrom, pos):
    """Gene the position sits in, else nearest within FLANK, else None."""
    c = genes[genes.chrom == chrom]
    inside = c[(c.start <= pos) & (pos <= c.end)]
    if len(inside):
        return inside.iloc[0]["gene"], 0
    near = c[(c.start - FLANK <= pos) & (pos <= c.end + FLANK)]
    if len(near):
        d = ((near.start - pos).abs()).combine((near.end - pos).abs(), min)
        i = d.idxmin()
        return near.loc[i, "gene"], int(d.loc[i])
    return None, None


def diamond_gene(best_hit):
    if not isinstance(best_hit, str):
        return None
    m = re.search(r"gene_symbol:(\S+)", best_hit)
    return m.group(1) if m else None


def main():
    overlay, gff = sys.argv[1], sys.argv[2]
    outdir = Path(sys.argv[3] if len(sys.argv) > 3 else "output")
    genes = load_genes(gff)
    print(f"{len(genes)} genes loaded")

    o = pd.read_csv(overlay, sep="\t", low_memory=False)
    loci = o["contig"].map(parse_locus)
    o["locus_chrom"] = [c for c, _ in loci]
    o["locus_pos"] = [p for _, p in loci]

    res = o["contig"].copy().to_frame()
    hits = [gene_at(genes, c, p) if c else (None, None)
            for c, p in zip(o["locus_chrom"], o["locus_pos"])]
    o["locus_gene"] = [g for g, _ in hits]
    o["locus_gene_dist"] = [d for _, d in hits]
    o["hit_gene"] = o["best_hit"].map(diamond_gene)
    # do the insertion-site gene and the sequence-hit gene agree?
    o["locus_matches_hit"] = (
        o["locus_gene"].notna() & o["hit_gene"].notna()
        & (o["locus_gene"].str.lower() == o["hit_gene"].str.lower())
    )

    out = outdir / "overlay_with_locus.tsv"
    o.to_csv(out, sep="\t", index=False)
    n_match = int(o["locus_matches_hit"].sum())
    n_ingene = int(o["locus_gene"].notna().sum())
    print(f"contigs with insertion in/near a gene: {n_ingene}")
    print(f"contigs where locus gene == DIAMOND hit gene: {n_match}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
