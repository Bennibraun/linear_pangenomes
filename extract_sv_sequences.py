import pysam
import argparse

# Sample order from your VCF header
SAMPLES = [
    "florida1", "florida2", "florida3", "florida4", 
    "thailand1", "thailand2", "thailand4", 
    "tokyo1", "tokyo2", "tokyo3"
]

def extract_flanked_insertions(vcf_path, ref_path, output_fasta, flank_size=500, min_ins=50):
    vcf = pysam.VariantFile(vcf_path)
    reference = pysam.FastaFile(ref_path)
    
    with open(output_fasta, "w") as out:
        count = 0
        seen_headers = {}
        for record in vcf:
            alt_seq = record.alts[0]
            ref_allele = record.ref
            
            # Filter for Insertions (ALT > REF and not symbolic)
            if len(alt_seq) > len(ref_allele) and not alt_seq.startswith('<'):
                ins_len = len(alt_seq) - len(ref_allele)
                if ins_len < min_ins:
                    continue

                # 1. Get the novel part
                novel_part = alt_seq[len(ref_allele):]
                
                # 2. Pull Reference Flanks
                chrom = record.chrom
                pos = record.pos
                
                start_flank = max(0, pos - flank_size)
                end_flank = pos + flank_size
                
                try:
                    left_flank = reference.fetch(chrom, start_flank, pos)
                    right_flank = reference.fetch(chrom, pos, end_flank)
                    
                    # 3. Identify Sample
                    supp_vec = record.info.get("SUPP_VEC", "")
                    sample_name = "unknown"
                    for i, val in enumerate(supp_vec):
                        if val == "1" and i < len(SAMPLES):
                            sample_name = SAMPLES[i]
                            break

                    # 4. Construct Chimeric Sequence
                    full_seq = left_flank + novel_part + right_flank
                    base_name = f"SV_{chrom}_{pos}_{sample_name}_ins{ins_len}"
                    # Ensure unique FASTA headers when multiple records map to the same name.
                    if base_name in seen_headers:
                        seen_headers[base_name] += 1
                        header = f">{base_name}_dup{seen_headers[base_name]}"
                    else:
                        seen_headers[base_name] = 0
                        header = f">{base_name}"
                    out.write(f"{header}\n{full_seq}\n")
                    count += 1
                except KeyError:
                    # Skip if chromosome name doesn't match reference
                    continue
                
    print(f"Extracted {count} flanked sequences to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcf", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-f", "--flank", type=int, default=500)
    args = parser.parse_args()
    
    extract_flanked_insertions(args.vcf, args.ref, args.out, args.flank)

# Usage: python 2.1.extract_svs.py -v /Users/bebr1814/projects/bee_paper/data/sv_calls/pan_sample_catalog/pan_sample_catalog.survivor.vcf -r /Users/bebr1814/projects/bee_paper/reference/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna --flank 200 -o /Users/bebr1814/projects/bee_paper/data/augref/extracted_flanked_sv_seqs.fasta