import pysam

vcf_path   = snakemake.input.vcf
ref_path   = snakemake.input.reference
out_fasta  = snakemake.output.fasta
flank_size = snakemake.params.flank
min_ins    = snakemake.params.min_ins
samples    = snakemake.params.samples.split()

vcf       = pysam.VariantFile(vcf_path)
reference = pysam.FastaFile(ref_path)

count = 0
seen_headers = {}
with open(out_fasta, "w") as out:
    for record in vcf:
        alt_seq    = record.alts[0]
        ref_allele = record.ref
        if len(alt_seq) <= len(ref_allele) or alt_seq.startswith("<"):
            continue
        ins_len = len(alt_seq) - len(ref_allele)
        if ins_len < min_ins:
            continue
        novel_part = alt_seq[len(ref_allele):]
        chrom = record.chrom
        pos   = record.pos
        try:
            left_flank  = reference.fetch(chrom, max(0, pos - flank_size), pos)
            right_flank = reference.fetch(chrom, pos, pos + flank_size)
        except KeyError:
            continue
        supp_vec    = record.info.get("SUPP_VEC", "")
        sample_name = "unknown"
        for i, val in enumerate(supp_vec):
            if val == "1" and i < len(samples):
                sample_name = samples[i]
                break
        full_seq  = left_flank + novel_part + right_flank
        base_name = "SV_" + chrom + "_" + str(pos) + "_" + sample_name + "_ins" + str(ins_len)
        if base_name in seen_headers:
            seen_headers[base_name] += 1
            header = ">" + base_name + "_dup" + str(seen_headers[base_name])
        else:
            seen_headers[base_name] = 0
            header = ">" + base_name
        out.write(header + "\n" + full_seq + "\n")
        count += 1

print("Extracted", count, "flanked sequences to", out_fasta)
