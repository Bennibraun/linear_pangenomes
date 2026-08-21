# Accessory-sequence exploration (standalone, for the paper)

One-off characterisation of what's *inside* the augref accessory (`SV_*`
insertion) sequence — coding potential, repeat content, and which contigs are
differentially present across populations (esp. cave-shared or surface-biased).
**Not part of the Snakemake pipeline** — run these by hand on the cluster.

## Question

Is there anything biologically interesting in the accessory sequence itself —
genes (not repeats) that are present/absent in a population-structured way, and
especially **shared across independent cave populations** (Molino / Pachón /
Tinaja)? That's the accessory-specific discovery a plain linear reference can't
make, and it's the biology payoff for the paper. Roback et al. is context, not a
coordinate join: a hit landing on one of their cave genes/pathways corroborates,
but we're not restricted to their list.

## Pipeline (2 steps)

```bash
conda env create -f env.yaml
conda activate accessory-explore

# Step 1: repeat-mask (dustmasker) + protein search (DIAMOND blastx)
#   -> edit the CONFIG block first (esp. the protein DB), then:
sbatch 01_mask_and_search.slurm

# Step 2: overlay population presence/absence + rank interesting contigs
sbatch 02_population_overlay.slurm
```

## What each step does

**01_mask_and_search.slurm**
- `dustmasker` the accessory FASTA (`extracted_flanked_sv_seqs.dedup.fasta`) →
  per-contig repeat/low-complexity fraction (`repeat_fraction.tsv`). This is a
  *filter*, not repeat modelling — fast, and we report "% of accessory sequence
  that's repetitive" as a result.
- Split into `repetitive_contigs.txt` (≥50% masked) and `candidate_contigs.txt`.
- `DIAMOND blastx` the candidate contigs vs a protein DB → `candidate_diamond.tsv`.

**02_population_overlay.py / .slurm**
- Per-contig presence/absence across the 72 samples, from the merged augref VCF
  (a sample "has" a contig if non-missing at ≥1 of its variant sites — a robust
  binary signal that survives the accessory low-coverage confound).
- Classifies each contig: `surface_biased`, `cave_biased`, `cave_shared`
  (present in ≥2 independent cave populations — the Roback-style reuse signal).
- Joins the DIAMOND best hit, flags TEs (transposase/gag-pol/etc.), and ranks so
  **genic (non-TE) + structured** contigs sort to the top.
- Outputs: `accessory_population_overlay.tsv` (all contigs, ranked) and
  `accessory_candidates_shortlist.tsv` (the genic + cave-shared/surface-biased
  shortlist to eyeball for the paper).

## You must provide: a protein DB

Set `DIAMOND_DB` (a prebuilt `.dmnd`) or `PROTEIN_FASTA` (built automatically) in
step 1's CONFIG. Recommended: a teleost proteome — **zebrafish (*Danio rerio*)
UniProt reference proteome** is a good default (well annotated, close enough to
find real fish genes and name them). A broader set (UniRef90) catches more but is
slower and noisier. The *A. mexicanus* proteome itself would mostly re-find
conspec genes; a related-species proteome is better for spotting genes in genuinely
novel accessory sequence.

## Interpreting hits (caveats to keep in mind)

- Contigs are **flanked** (200 bp conspec on each side of the novel insertion).
  A DIAMOND hit confined to the flanks just re-finds a conspec gene — the novel
  content is the middle (`ins<len>` in the contig name gives the novel length).
  Prefer hits whose `qstart..qend` fall inside `[200, length-200]`.
- Expect a large repetitive fraction — that's normal for accessory insertion
  sequence, and reporting it is itself a result.
- TEs that slip past dustmasker will hit transposase/gag-pol in the protein DB
  and get `looks_like_TE=True` — filtered from the shortlist.
- Presence/absence is robust; don't over-read exact per-population fractions at
  low accessory depth (see the paper's QC section).

## Inputs assumed (cluster paths, edit in the CONFIG blocks)

- accessory FASTA: `results/sv_calls/augref/extracted_flanked_sv_seqs.dedup.fasta`
- merged VCF: `results/variants/bcftools/augref/combined/merged.vcf.gz`
- manifest: `reads_manifest.tsv` (columns `sample_id`, `grouping`)
