"""Consolidate per-sample QUAST + BUSCO outputs into one TSV.

Inputs:
    snakemake.params.samples       List[str] of long-read sample IDs
    snakemake.params.quast_dirs    List[str] paralleled with samples
    snakemake.params.busco_dirs    List[str] paralleled with samples
    snakemake.params.lineage       BUSCO lineage name (for filename matching)
Output:
    snakemake.output.summary       single TSV
"""
import json
import re
from pathlib import Path

samples     = list(snakemake.params.samples)
quast_dirs  = list(snakemake.params.quast_dirs)
busco_dirs  = list(snakemake.params.busco_dirs)
lineage     = snakemake.params.lineage
out_path    = Path(snakemake.output.summary)


def parse_quast(quast_dir: Path) -> dict:
    """Pull a few headline metrics from QUAST report.tsv."""
    report = quast_dir / "report.tsv"
    if not report.exists():
        return {}
    metrics = {}
    for line in report.read_text().splitlines():
        if "\t" not in line:
            continue
        key, _, val = line.partition("\t")
        metrics[key.strip()] = val.strip()
    return {
        "n_contigs":      metrics.get("# contigs", "NA"),
        "largest_contig": metrics.get("Largest contig", "NA"),
        "total_length":   metrics.get("Total length", "NA"),
        "gc_pct":         metrics.get("GC (%)", "NA"),
        "n50":            metrics.get("N50", "NA"),
        "l50":            metrics.get("L50", "NA"),
        "n75":            metrics.get("N75", "NA"),
        "n_per_100kb":    metrics.get("# N's per 100 kbp", "NA"),
    }


def parse_busco(busco_dir: Path, lineage: str, sample: str) -> dict:
    """Pull C/S/D/F/M percentages from BUSCO short_summary.json (preferred) or
    fall back to the .txt summary."""
    candidates = list(busco_dir.glob("short_summary.specific.*.json"))
    if not candidates:
        candidates = list(busco_dir.glob("short_summary.*.json"))
    for cand in candidates:
        try:
            data    = json.loads(cand.read_text())
            results = data.get("results", {})
            return {
                "busco_complete":           results.get("Complete percentage", "NA"),
                "busco_single_copy":        results.get("Single copy percentage", "NA"),
                "busco_duplicated":         results.get("Multi copy percentage",
                                                        results.get("Duplicated percentage", "NA")),
                "busco_fragmented":         results.get("Fragmented percentage", "NA"),
                "busco_missing":            results.get("Missing percentage", "NA"),
                "busco_total":              results.get("n_markers",
                                                        results.get("Total markers", "NA")),
                "busco_lineage":            data.get("lineage_dataset", {}).get("name", lineage),
            }
        except (json.JSONDecodeError, OSError):
            continue

    txt_candidates = list(busco_dir.glob("short_summary.specific.*.txt"))
    if not txt_candidates:
        txt_candidates = list(busco_dir.glob("short_summary.*.txt"))
    for cand in txt_candidates:
        text = cand.read_text()
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", text)
        if m:
            return {
                "busco_complete":    m.group(1),
                "busco_single_copy": m.group(2),
                "busco_duplicated":  m.group(3),
                "busco_fragmented":  m.group(4),
                "busco_missing":     m.group(5),
                "busco_total":       m.group(6),
                "busco_lineage":     lineage,
            }
    return {}


cols = [
    "sample",
    "n_contigs", "largest_contig", "total_length", "gc_pct",
    "n50", "l50", "n75", "n_per_100kb",
    "busco_complete", "busco_single_copy", "busco_duplicated",
    "busco_fragmented", "busco_missing", "busco_total", "busco_lineage",
]

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as fh:
    fh.write("\t".join(cols) + "\n")
    for sample, qdir, bdir in zip(samples, quast_dirs, busco_dirs):
        row = {"sample": sample}
        row.update(parse_quast(Path(qdir)))
        row.update(parse_busco(Path(bdir), lineage, sample))
        fh.write("\t".join(str(row.get(c, "NA")) for c in cols) + "\n")
