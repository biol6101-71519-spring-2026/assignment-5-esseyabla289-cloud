#!/usr/bin/env python3
"""
Concatenate trimmed FASTA alignments into a supermatrix.
Writes the concatenated alignment and an RAxML-style partition file.

Called by Snakemake; receives input/output/params via the snakemake object.
"""

import sys
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger()


def read_fasta(filepath):
    """Parse a FASTA file into {header: sequence} preserving insertion order."""
    seqs = {}
    name = None
    buf  = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].strip()
                buf  = []
            elif line:
                buf.append(line.upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(seqs, filepath, width=60):
    with open(filepath, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


# ── Inputs from Snakemake ──────────────────────────────────────────────────
aln_files     = snakemake.input.alns
genes         = list(snakemake.params.genes)
concat_out    = snakemake.output.concat
partitions_out = snakemake.output.partitions

log.info("Reading %d alignments: %s", len(genes), genes)

# Load each alignment
gene_alns   = {}
aln_lengths = {}
for gene, aln_file in zip(genes, aln_files):
    seqs = read_fasta(aln_file)
    gene_alns[gene]   = seqs
    aln_lengths[gene] = len(next(iter(seqs.values())))
    log.info("  %s: %d taxa, %d bp", gene, len(seqs), aln_lengths[gene])

# Union of all taxa (sorted for determinism)
all_taxa = sorted({t for seqs in gene_alns.values() for t in seqs})
log.info("Total taxa: %d", len(all_taxa))

# Build concatenated sequences and partition boundaries
concat_seqs = {taxon: "" for taxon in all_taxa}
partitions  = []
pos         = 1

for gene in genes:
    seqs    = gene_alns[gene]
    aln_len = aln_lengths[gene]
    end_pos = pos + aln_len - 1
    partitions.append(f"DNA, {gene} = {pos}-{end_pos}")
    pos += aln_len

    for taxon in all_taxa:
        if taxon in seqs:
            concat_seqs[taxon] += seqs[taxon]
        else:
            # Fill with gap/missing data for taxa absent from this locus
            concat_seqs[taxon] += "-" * aln_len
            log.warning("  Taxon '%s' absent from %s — filled with gaps", taxon, gene)

total_len = pos - 1
log.info("Concatenated alignment: %d taxa × %d bp", len(all_taxa), total_len)

# Write outputs
write_fasta(concat_seqs, concat_out)
log.info("Wrote concatenated alignment: %s", concat_out)

with open(partitions_out, "w") as fh:
    fh.write("\n".join(partitions) + "\n")
log.info("Wrote partition file: %s", partitions_out)

print(f"Concatenated {len(genes)} genes → {total_len} bp total")
for p in partitions:
    print(f"  {p}")
