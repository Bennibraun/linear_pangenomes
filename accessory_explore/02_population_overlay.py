#!/usr/bin/env python3
"""Accessory-sequence exploration, step 2: overlay population presence/absence
onto the candidate contigs and their protein hits.

The biology we're after: accessory (SV_*) contigs that (a) carry coding
sequence and (b) are differentially present across populations -- especially
present in surface but absent/reduced in cave (the polarity of a Roback-style
cave deletion), and/or SHARED across independent cave populations.

"Presence" of a contig in a sample = the sample has a non-missing genotype at
one or more of the contig's variant sites in the merged augref VCF. This is a
robust binary signal (survives the low-coverage/missingness issues that plague
accessory allele-frequency stats -- see the paper's QC section).

Inputs (all already on the cluster after 01_mask_and_search.slurm):
  --vcf         merged augref VCF (results/variants/bcftools/augref/combined/merged.vcf.gz)
  --diamond     candidate_diamond.tsv from step 1
  --repeat-frac repeat_fraction.tsv from step 1
  --manifest    reads_manifest.tsv (sample_id, grouping) to map sample->population
Output:
  a per-contig table: presence by population, morph pattern, best protein hit,
  and a "cave_shared" / "surface_biased" flag -- ranked so the interesting
  contigs (genic + cave-shared or surface-biased) sort to the top.

Presence-call detail: SV_* contig names are CHROM in the VCF, so we group the
VCF's SV_* records by CHROM and, per sample, mark the contig present if the
sample is non-missing at ANY of its sites.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

# Population -> morph. Adjust if population naming differs in your manifest.
SURFACE_POPS = {"Riochoy", "Rio_choy", "Choy", "Mante", "Rascon"}


def morph_of(pop):
    return "surface" if pop in SURFACE_POPS or pop.lower().startswith(("choy", "rio", "mante", "rascon")) else "cave"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="merged augref VCF (bgzipped, indexed)")
    p.add_argument("--diamond", required=True, help="candidate_diamond.tsv from step 1")
    p.add_argument("--repeat-frac", required=True, help="repeat_fraction.tsv from step 1")
    p.add_argument("--manifest", required=True, help="reads_manifest.tsv (sample_id, grouping)")
    p.add_argument("--outdir", default="results/accessory_explore")
    p.add_argument("--present-frac", type=float, default=0.5,
                   help="a population 'has' a contig if this fraction of its samples are present (default 0.5)")
    return p.parse_args()


def load_sample_pop(manifest_path):
    df = pd.read_csv(manifest_path, sep="\t")
    # tolerate common column names
    id_col = "sample_id" if "sample_id" in df.columns else df.columns[0]
    grp_col = "grouping" if "grouping" in df.columns else "population"
    return dict(zip(df[id_col], df[grp_col]))


def contig_presence(vcf_path, sample_pop):
    """Return DataFrame: index=SV_* contig, columns=samples, values=presence(0/1).
    Present = sample non-missing at >=1 site on that contig."""
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    # only samples we have a population for
    keep = [s for s in samples if s in sample_pop]
    present = defaultdict(lambda: np.zeros(len(keep), dtype=bool))
    sidx = {s: i for i, s in enumerate(keep)}

    for rec in vcf:
        chrom = rec.chrom
        if not chrom.startswith("SV_"):
            continue
        row = present[chrom]
        for s in keep:
            gt = rec.samples[s].get("GT")
            # non-missing if any allele is called (not all None)
            if gt is not None and any(a is not None for a in gt):
                row[sidx[s]] = True
    df = pd.DataFrame(present, index=keep).T  # contigs x samples
    df.columns = keep
    return df


def summarize(pres, sample_pop, present_frac):
    """Per-contig: fraction present in each population, morph presence, patterns."""
    samples = pres.columns.tolist()
    pops = pd.Series({s: sample_pop[s] for s in samples})
    pop_names = sorted(pops.unique())
    rows = []
    for contig, row in pres.iterrows():
        rec = {"contig": contig}
        pop_present = {}
        for pop in pop_names:
            members = pops[pops == pop].index
            frac = row[members].mean() if len(members) else np.nan
            rec[f"present_frac__{pop}"] = round(float(frac), 3)
            pop_present[pop] = frac >= present_frac
        surface_pops = [p for p in pop_names if morph_of(p) == "surface"]
        cave_pops = [p for p in pop_names if morph_of(p) == "cave"]
        surf_hit = [p for p in surface_pops if pop_present[p]]
        cave_hit = [p for p in cave_pops if pop_present[p]]
        rec["n_surface_present"] = len(surf_hit)
        rec["n_cave_present"] = len(cave_hit)
        rec["cave_pops_present"] = ",".join(cave_hit)
        rec["surface_pops_present"] = ",".join(surf_hit)
        # patterns
        rec["surface_biased"] = len(surf_hit) > 0 and len(cave_hit) == 0
        rec["cave_biased"] = len(cave_hit) > 0 and len(surf_hit) == 0
        # "cave_shared" = present in >=2 INDEPENDENT cave populations (the
        # Roback-style reuse signal). Molino / Pachon / Tinaja are the key
        # independent lineages.
        rec["cave_shared"] = len(cave_hit) >= 2
        rows.append(rec)
    return pd.DataFrame(rows)


def best_hits(diamond_path):
    """Best protein hit per contig from DIAMOND outfmt 6, and a TE flag."""
    cols = ["qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    if not Path(diamond_path).exists() or Path(diamond_path).stat().st_size == 0:
        return pd.DataFrame(columns=["contig", "best_hit", "best_pident",
                                     "best_evalue", "looks_like_TE"])
    d = pd.read_csv(diamond_path, sep="\t", names=cols)
    d = d.sort_values(["qseqid", "bitscore"], ascending=[True, False])
    top = d.groupby("qseqid").first().reset_index()
    te_pat = re.compile(
        r"transpos|reverse transcript|gag-pol|retrovir|integrase|helitron|"
        r"pol polyprotein|LINE-|copia|gypsy|mariner|tc1|DDE", re.I)
    top["looks_like_TE"] = top["stitle"].fillna("").str.contains(te_pat)
    return top.rename(columns={
        "qseqid": "contig", "stitle": "best_hit",
        "pident": "best_pident", "evalue": "best_evalue"})[
        ["contig", "best_hit", "best_pident", "best_evalue", "looks_like_TE"]]


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sample_pop = load_sample_pop(args.manifest)

    print("Computing per-contig presence/absence from the merged VCF ...")
    pres = contig_presence(args.vcf, sample_pop)
    print(f"  {pres.shape[0]} SV_* contigs x {pres.shape[1]} samples")

    summ = summarize(pres, sample_pop, args.present_frac)
    rep = pd.read_csv(args.repeat_frac, sep="\t")[["contig", "length", "masked_frac"]]
    hits = best_hits(args.diamond)

    out = summ.merge(rep, on="contig", how="left").merge(hits, on="contig", how="left")
    out["has_protein_hit"] = out["best_hit"].notna()
    out["genic_nonTE"] = out["has_protein_hit"] & (out["looks_like_TE"] != True)  # noqa: E712

    # Rank: interesting = genic (non-TE) AND (cave_shared OR surface_biased).
    out["interest_score"] = (
        out["genic_nonTE"].astype(int) * 2
        + out["cave_shared"].astype(int)
        + out["surface_biased"].astype(int)
    )
    out = out.sort_values(
        ["interest_score", "n_cave_present", "best_pident"],
        ascending=[False, False, False])

    full = outdir / "accessory_population_overlay.tsv"
    out.to_csv(full, sep="\t", index=False)
    print(f"wrote {full}  ({len(out)} contigs)")

    # A focused shortlist for the paper: genic, non-TE, and structured.
    short = out[out["genic_nonTE"] & (out["cave_shared"] | out["surface_biased"])]
    short_path = outdir / "accessory_candidates_shortlist.tsv"
    short.to_csv(short_path, sep="\t", index=False)
    print(f"wrote {short_path}  ({len(short)} candidate genic + structured contigs)")

    # Quick console summary.
    print("\n=== summary ===")
    print(f"  contigs with protein hit:      {int(out['has_protein_hit'].sum())}")
    print(f"  genic (non-TE) contigs:        {int(out['genic_nonTE'].sum())}")
    print(f"  surface-biased presence:       {int(out['surface_biased'].sum())}")
    print(f"  cave-shared (>=2 cave pops):   {int(out['cave_shared'].sum())}")
    print(f"  genic + cave-shared:           {int((out['genic_nonTE'] & out['cave_shared']).sum())}")
    print(f"  genic + surface-biased:        {int((out['genic_nonTE'] & out['surface_biased']).sum())}")


if __name__ == "__main__":
    main()
