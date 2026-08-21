# Accessory-sequence exploration

What's *inside* the augref accessory (`SV_*` insertion) sequence — genes (not
repeats) that are present/absent in a population-structured way, especially
shared across independent cave populations (Molino/Pachón/Tinaja). Standalone,
not part of the Snakemake pipeline.

## Run

```bash
mamba env create -f env.yaml
mamba activate accessory-explore

bash 00_get_proteome.sh          # zebrafish proteome + DIAMOND db (needs internet)
sbatch 01_mask_and_search.slurm  # dustmasker repeat filter -> DIAMOND blastx
sbatch 02_population_overlay.slurm
```

Edit paths at the top of each script if yours differ (accessory FASTA, merged
VCF, `reads_manifest.tsv`, and `DIAMOND_DB` in step 1).

## Outputs (in `results/accessory_explore/`)

- `repeat_fraction.tsv` — per-contig dustmasker masked fraction.
- `candidate_diamond.tsv` — protein hits for non-repetitive contigs.
- `overlay_all.tsv` — every contig: presence per population, cave_shared /
  surface_biased flags, best hit, ranked.
- `shortlist.tsv` — genic (non-TE) + structured contigs. The list to eyeball.

## Notes

- Contigs are flanked (200 bp conspec each side). A hit only in the flanks
  re-finds a conspec gene; the novel part is the middle (`ins<len>` in the name).
- Expect lots of repeats — that's normal for accessory sequence, and the masked
  fraction is itself a reportable number.
- Presence/absence is robust to low accessory depth; don't over-read exact
  per-population fractions.
