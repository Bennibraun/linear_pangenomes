import pysam
from pathlib import Path

samples   = snakemake.params.samples
vcf_paths = snakemake.input.per_sample_vcfs
out_path  = snakemake.output.summary

rows = ["sample\tn_insertions\tnovel_bp\tn_deletions\tdel_bp\tn_other\ttotal_svs"]
for sample, vcf_path in zip(samples, vcf_paths):
    vcf = pysam.VariantFile(vcf_path)
    n_ins = n_del = n_other = 0
    ins_bp = del_bp = 0
    for rec in vcf.fetch():
        svtype = rec.info.get("SVTYPE", "")
        svlen  = abs(rec.info.get("SVLEN", 0))
        if isinstance(svlen, tuple):
            svlen = abs(svlen[0]) if svlen else 0
        if svtype == "INS":
            n_ins += 1
            ins_bp += svlen
        elif svtype == "DEL":
            n_del += 1
            del_bp += svlen
        else:
            n_other += 1
    rows.append("\t".join([sample, str(n_ins), str(ins_bp), str(n_del), str(del_bp), str(n_other), str(n_ins + n_del + n_other)]))

Path(out_path).parent.mkdir(parents=True, exist_ok=True)
Path(out_path).write_text("\n".join(rows) + "\n")
