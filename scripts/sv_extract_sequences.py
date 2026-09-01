"""Emit one FASTA record per insertion SV: the BARE novel inserted sequence,
no reference flanks.

Earlier this padded each insert with `flank` bp of conspec sequence on each
side. That was a mistake: the flanks are duplicated core, so (a) short reads
multi-map between the accessory contig and the real conspec locus, and (b) any
FST/SFS/PAV computed over the contig partly measures core sequence, which is
exactly why the accessory signal collapsed toward the core. Bare inserts map
cleanly (a read places on the insert or on conspec, not both) and every base of
an accessory contig is genuinely novel sequence.

Trade-off accepted: inserts shorter than a read can't be uniquely placed by
short reads and will simply go ungenotyped. That is honest -- those variants
aren't short-read genotypable -- rather than papered over with core padding.
"""
import pysam

vcf_path  = snakemake.input.vcf
out_fasta = snakemake.output.fasta
min_ins   = snakemake.params.min_ins
samples   = snakemake.params.samples.split()

vcf = pysam.VariantFile(vcf_path)

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

        supp_vec    = record.info.get("SUPP_VEC", "")
        sample_name = "unknown"
        for i, val in enumerate(supp_vec):
            if val == "1" and i < len(samples):
                sample_name = samples[i]
                break

        base_name = "SV_" + record.chrom + "_" + str(record.pos) + "_" + sample_name + "_ins" + str(ins_len)
        if base_name in seen_headers:
            seen_headers[base_name] += 1
            header = ">" + base_name + "_dup" + str(seen_headers[base_name])
        else:
            seen_headers[base_name] = 0
            header = ">" + base_name
        out.write(header + "\n" + novel_part + "\n")
        count += 1

print("Extracted", count, "bare insertion sequences to", out_fasta)
