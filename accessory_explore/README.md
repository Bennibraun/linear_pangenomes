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

# Optional but recommended: confirm insertion loci vs the conspec annotation.
# conspec = GCF_023375975.1_AstMex3_surface.
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_genomic.gff.gz
python 03_gff_intersect.py results/accessory_explore/overlay_all.tsv \
  GCF_023375975.1_AstMex3_surface_genomic.gff.gz results/accessory_explore
```

`overlay_all.tsv` now carries per-population presence FRACTIONS (`frac_<pop>`
columns) as well as the boolean flags — use those to sanity-check calls (a
"cave-specific" contig assembled from a surface donor with high surface
`frac_*` is a false call, e.g. the OCA2 artifact).

Edit paths at the top of each script if yours differ (accessory FASTA, merged
VCF, `reads_manifest.tsv`, and `DIAMOND_DB` in step 1).

## Outputs (in `results/accessory_explore/`)

- `repeat_fraction.tsv` — per-contig dustmasker masked fraction.
- `candidate_diamond.tsv` — protein hits for non-repetitive contigs.
- `overlay_all.tsv` — every contig: per-population presence fractions
  (`frac_<pop>`) + flags, best hit, ranked.
- `shortlist.tsv` — genic (non-TE) + structured contigs. The list to eyeball.
- `overlay_with_locus.tsv` (after step 3) — adds the insertion-site gene
  (`locus_gene`) and whether it matches the DIAMOND hit gene
  (`locus_matches_hit`). A match = an X-like insertion sitting at the X locus.

## Notes

- Contigs are flanked (200 bp conspec each side). A hit only in the flanks
  re-finds a conspec gene; the novel part is the middle (`ins<len>` in the name).
- Expect lots of repeats — that's normal for accessory sequence, and the masked
  fraction is itself a reportable number.
- Presence/absence is robust to low accessory depth; don't over-read exact
  per-population fractions.
