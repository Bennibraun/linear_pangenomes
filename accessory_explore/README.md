# Accessory-sequence exploration

What's *inside* the augref accessory (`SV_*` insertion) sequence — genes (not
repeats) that are present/absent in a population-structured way, especially
shared across independent cave populations (Molino/Pachón/Tinaja).

**This is a SEPARATE analysis** that READS the Snakemake pipeline's results and
writes everything into its own `OUTDIR`. It never writes into the pipeline tree.
All paths live in `config.sh` — set them once.

## Setup

Edit `config.sh`:
- `PIPELINE_RESULTS` — absolute path to the finished pipeline's `results/` dir.
- `OUTDIR` — where this analysis writes (default `output`, local).
- `MANIFEST` — the sample→population TSV.

```bash
mamba env create -f env.yaml
mamba activate accessory-explore
```

## Run

```bash
bash 00_get_proteome.sh          # zebrafish proteome + DIAMOND db (needs internet)
sbatch 01_mask_and_search.slurm  # dustmasker repeat filter -> DIAMOND blastx
sbatch 02_population_overlay.slurm

# insertion-locus check vs the conspec annotation (GCF_023375975.1_AstMex3_surface)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_genomic.gff.gz
python 03_gff_intersect.py "$OUTDIR/overlay_all.tsv" GCF_023375975.1_AstMex3_surface_genomic.gff.gz "$OUTDIR"
```

(For 03, `$OUTDIR` is whatever you set in config.sh, e.g. `output`.)

## Resuming an earlier run (outputs written to the wrong place)

An earlier version of these scripts wrote outputs into the pipeline's
`results/accessory_explore/`. To move them into this analysis's `OUTDIR` and
resume WITHOUT recomputing 01/02/03:

```bash
bash migrate_outputs.sh          # moves them into $OUTDIR
# then, once confirmed:  rm -rf "$PIPELINE_RESULTS/accessory_explore"
```

After migrating, everything (`candidate_diamond.tsv`, `repeat_fraction.tsv`,
`overlay_all.tsv`, `overlay_with_locus.tsv`, …) is in `OUTDIR`; nothing needs
rerunning.

## Outputs (in `$OUTDIR`)

- `repeat_fraction.tsv` — per-contig dustmasker masked fraction.
- `candidate_diamond.tsv` — protein hits for non-repetitive contigs.
- `overlay_all.tsv` — every contig: per-population presence fractions
  (`frac_<pop>`) + flags, best hit, ranked.
- `shortlist.tsv` — genic (non-TE) + structured contigs.
- `overlay_with_locus.tsv` (after 03) — adds insertion-site gene (`locus_gene`)
  and whether it matches the DIAMOND hit gene (`locus_matches_hit`).

## Notes

- Contigs are flanked (200 bp conspec each side). A hit only in the flanks
  re-finds a conspec gene; the novel part is the middle (`ins<len>` in the name).
- Expect lots of repeats — normal for accessory sequence; masked fraction is a
  reportable number.
- Presence/absence is robust to low accessory depth; don't over-read exact
  per-population fractions. A "cave-specific" contig assembled from a surface
  donor but with high surface `frac_*` is a false call.
