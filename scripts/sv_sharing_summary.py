import pysam
from pathlib import Path
from collections import defaultdict

samples          = snakemake.params.samples
survivor_vcf     = snakemake.input.survivor_vcf
out_path         = snakemake.output.summary
n                = len(samples)

per_sample_total  = defaultdict(int)
per_sample_unique = defaultdict(int)
per_sample_shared = defaultdict(int)
svtype_counts     = defaultdict(lambda: defaultdict(int))

vcf = pysam.VariantFile(survivor_vcf)
for rec in vcf.fetch():
    supp_vec = rec.info.get("SUPP_VEC", "")
    if not supp_vec:
        continue
    supporters = [i for i, c in enumerate(supp_vec) if c == "1"]
    svtype = rec.info.get("SVTYPE", "OTHER")
    for i in supporters:
        if i < n:
            per_sample_total[samples[i]] += 1
            svtype_counts[samples[i]][svtype] += 1
            if len(supporters) == 1:
                per_sample_unique[samples[i]] += 1
            else:
                per_sample_shared[samples[i]] += 1

all_svtypes   = sorted({st for s in svtype_counts for st in svtype_counts[s]})
svtype_header = "\t".join("n_" + st.lower() for st in all_svtypes)

rows = ["sample\ttotal_svs\tunique_svs\tshared_svs\t" + svtype_header]
for s in samples:
    svtype_cols = "\t".join(str(svtype_counts[s].get(st, 0)) for st in all_svtypes)
    rows.append("\t".join([s, str(per_sample_total[s]), str(per_sample_unique[s]), str(per_sample_shared[s]), svtype_cols]))

Path(out_path).parent.mkdir(parents=True, exist_ok=True)
Path(out_path).write_text("\n".join(rows) + "\n")
