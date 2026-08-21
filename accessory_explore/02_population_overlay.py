#!/usr/bin/env python3
"""Step 2: which accessory contigs are genic and population-structured?

For each SV_* contig: presence per population (a sample "has" it if non-missing
at >=1 of its variant sites in the merged VCF -- robust to low accessory depth),
its best protein hit (TEs flagged), then rank genic + cave-shared/surface-biased
contigs to the top.

Usage:
  02_population_overlay.py VCF DIAMOND_TSV REPEAT_FRAC_TSV MANIFEST_TSV [OUTDIR]
"""
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

SURFACE = ("choy", "rio", "mante", "rascon")
PRESENT_FRAC = 0.5   # a population "has" a contig if >= this fraction of its samples do
TE = re.compile(r"transpos|reverse transcript|gag-pol|retrovir|integrase|"
                r"helitron|polyprotein|LINE|copia|gypsy|mariner|DDE", re.I)


def is_surface(pop):
    return pop.lower().startswith(SURFACE)


def contig_presence(vcf_path, samples):
    """{contig: bool array over samples}, present = non-missing at >=1 site."""
    vcf = pysam.VariantFile(vcf_path)
    idx = {s: i for i, s in enumerate(samples)}
    pres = {}
    for rec in vcf:
        if not rec.chrom.startswith("SV_"):
            continue
        row = pres.setdefault(rec.chrom, np.zeros(len(samples), bool))
        for s, smp in rec.samples.items():
            gt = smp.get("GT")
            if gt and any(a is not None for a in gt):
                row[idx[s]] = True
    return pres


def best_hits(path):
    cols = ["contig", "sseqid", "stitle", "pident", "length", "mm", "go",
            "qs", "qe", "ss", "se", "evalue", "bits"]
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=["contig", "best_hit", "best_pident", "is_TE"])
    d = pd.read_csv(path, sep="\t", names=cols)
    top = d.sort_values("bits", ascending=False).groupby("contig", as_index=False).first()
    top["is_TE"] = top["stitle"].fillna("").str.contains(TE)
    return top.rename(columns={"stitle": "best_hit", "pident": "best_pident"})[
        ["contig", "best_hit", "best_pident", "is_TE"]]


def main():
    vcf, diamond, repeat_frac, manifest = sys.argv[1:5]
    outdir = Path(sys.argv[5] if len(sys.argv) > 5 else "results/accessory_explore")
    outdir.mkdir(parents=True, exist_ok=True)

    man = pd.read_csv(manifest, sep="\t")
    sample_pop = dict(zip(man["sample_id"], man["grouping"]))
    samples = [s for s in pysam.VariantFile(vcf).header.samples if s in sample_pop]
    pop = pd.Series({s: sample_pop[s] for s in samples})

    pres = contig_presence(vcf, samples)
    print(f"{len(pres)} SV_* contigs x {len(samples)} samples")

    rows = []
    for contig, row in pres.items():
        r = {"contig": contig}
        present = {p: row[[i for i, s in enumerate(samples) if pop[s] == p]].mean() >= PRESENT_FRAC
                   for p in pop.unique()}
        cave = [p for p in present if not is_surface(p) and present[p]]
        surf = [p for p in present if is_surface(p) and present[p]]
        r.update(n_cave=len(cave), n_surface=len(surf),
                 cave_pops=",".join(cave), surface_pops=",".join(surf),
                 surface_biased=bool(surf and not cave),
                 cave_shared=len(cave) >= 2)   # present in >=2 cave pops = Roback-style reuse
        rows.append(r)
    summ = pd.DataFrame(rows)

    rep = pd.read_csv(repeat_frac, sep="\t")[["contig", "length", "masked_frac"]]
    out = summ.merge(rep, on="contig", how="left").merge(best_hits(diamond), on="contig", how="left")
    out["genic"] = out["best_hit"].notna() & ~out["is_TE"].fillna(False)
    out["score"] = 2 * out["genic"] + out["cave_shared"] + out["surface_biased"]
    out = out.sort_values(["score", "n_cave", "best_pident"], ascending=False)

    out.to_csv(outdir / "overlay_all.tsv", sep="\t", index=False)
    shortlist = out[out["genic"] & (out["cave_shared"] | out["surface_biased"])]
    shortlist.to_csv(outdir / "shortlist.tsv", sep="\t", index=False)

    print(f"genic contigs: {out.genic.sum()} | "
          f"genic+cave_shared: {(out.genic & out.cave_shared).sum()} | "
          f"genic+surface_biased: {(out.genic & out.surface_biased).sum()}")
    print(f"wrote overlay_all.tsv ({len(out)}) and shortlist.tsv ({len(shortlist)})")


if __name__ == "__main__":
    main()
