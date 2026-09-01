"""Extract assembly sequence that does NOT align to the conspec reference --
the accessory content Jeon's linear pangenome is built from, which our
insertion-only SV catalog misses entirely.

For one assembly and its minimap2 PAF vs conspec, we walk each query (assembly)
contig, mark the intervals covered by primary alignments, and emit every
uncovered run >= min_len as its own FASTA record. This captures both wholly
novel contigs (no alignment at all) and large unaligned segments of contigs
that otherwise map -- e.g. a big insertion carried on an assembled contig.

Records are named UNMAP_<sample>_<contig>_<start>_<end> (0-based, end-exclusive)
so provenance and coordinates are recoverable.
"""
import pysam

paf_path   = snakemake.input.paf
asm_path   = snakemake.input.assembly
out_fasta  = snakemake.output.fasta
sample     = snakemake.params.sample
min_len    = int(snakemake.params.min_len)
min_mapq   = int(snakemake.params.min_mapq)

# PAF columns: qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq ...
# Collect covered query intervals per contig from primary alignments.
covered = {}      # qname -> list of (qstart, qend)
qlen = {}         # qname -> length
with open(paf_path) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qname = f[0]
        qlen[qname] = int(f[1])
        if int(f[11]) < min_mapq:      # ignore ambiguous alignments
            continue
        covered.setdefault(qname, []).append((int(f[2]), int(f[3])))


def uncovered_runs(length, intervals):
    """Yield (start, end) gaps >= min_len not covered by any interval."""
    if not intervals:
        if length >= min_len:
            yield (0, length)
        return
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    cursor = 0
    for s, e in merged:
        if s - cursor >= min_len:
            yield (cursor, s)
        cursor = max(cursor, e)
    if length - cursor >= min_len:
        yield (cursor, length)


asm = pysam.FastaFile(asm_path)
count = 0
with open(out_fasta, "w") as out:
    for contig in asm.references:
        length = asm.get_reference_length(contig)
        for start, end in uncovered_runs(length, covered.get(contig, [])):
            seq = asm.fetch(contig, start, end)
            header = ">UNMAP_" + sample + "_" + contig + "_" + str(start) + "_" + str(end)
            out.write(header + "\n" + seq + "\n")
            count += 1

print("Extracted", count, "unmapped segments (>=", min_len, "bp) from", sample)
