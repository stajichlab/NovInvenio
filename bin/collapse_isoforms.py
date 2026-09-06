#!/usr/bin/env python3
"""
Collapse a protein FASTA to one representative record per gene, using an
NCBI RefSeq "feature_table.txt.gz" for the real gene<->protein mapping
(never guessed from FASTA header text/description, which is not a
reliable isoform-grouping key across species/annotation pipelines).

For each gene, the CDS feature with the greatest product_length (amino
acids) is kept as that gene's representative; every other isoform's
protein is dropped. A protein present in the FASTA but absent from the
feature table (or vice versa) is reported on stderr and simply excluded
from the collapsed output -- not a fatal error, since the two files can
legitimately be from slightly different processing passes.

Usage:
    bin/collapse_isoforms.py \\
        --feature-table GCF_000002985.6_WBcel235_feature_table.txt.gz \\
        --protein-fasta Caenorhabditis_elegans.proteins.fa \\
        --output Caenorhabditis_elegans.proteins.collapsed.fa
"""
import argparse
import csv
import gzip
import sys


def load_gene_map(feature_table_path) -> dict[str, tuple[str, int]]:
    """product_accession -> (GeneID, product_length), from CDS rows only."""
    opener = gzip.open if str(feature_table_path).endswith('.gz') else open
    mapping = {}
    with opener(feature_table_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        # NCBI's header line starts with "# feature", not "feature" -- DictReader
        # takes it as-is, so index the first column back out by position instead
        # of relying on the leading "# " surviving into the fieldname.
        feature_col = reader.fieldnames[0]
        for row in reader:
            if row[feature_col] != 'CDS':
                continue
            protein_id = row['product_accession'].strip()
            gene_id = row['GeneID'].strip()
            length = row['product_length'].strip()
            if not protein_id or not gene_id or not length:
                continue
            mapping[protein_id] = (gene_id, int(length))
    return mapping


def read_fasta_records(path):
    """Yield (header_id, full_header_line, sequence_lines) per record."""
    header_id = None
    header_line = None
    seq_lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header_id is not None:
                    yield header_id, header_line, seq_lines
                header_line = line
                header_id = line[1:].split()[0].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header_id is not None:
            yield header_id, header_line, seq_lines


def collapse(feature_table_path, protein_fasta_path, output_path) -> dict:
    gene_map = load_gene_map(feature_table_path)

    best_per_gene: dict[str, tuple[str, int]] = {}  # gene_id -> (protein_id, length)
    records: dict[str, tuple[str, list[str]]] = {}  # protein_id -> (header_line, seq_lines)
    unmapped = []

    for protein_id, header_line, seq_lines in read_fasta_records(protein_fasta_path):
        records[protein_id] = (header_line, seq_lines)
        if protein_id not in gene_map:
            unmapped.append(protein_id)
            continue
        gene_id, length = gene_map[protein_id]
        current = best_per_gene.get(gene_id)
        if current is None or length > current[1]:
            best_per_gene[gene_id] = (protein_id, length)

    kept_protein_ids = {protein_id for protein_id, _ in best_per_gene.values()}

    with open(output_path, 'w') as out:
        for protein_id in sorted(kept_protein_ids):
            header_line, seq_lines = records[protein_id]
            out.write(header_line)
            out.writelines(seq_lines)

    return {
        'input_records': len(records),
        'genes': len(best_per_gene),
        'kept': len(kept_protein_ids),
        'dropped_isoforms': len(records) - len(kept_protein_ids) - len(unmapped),
        'unmapped_no_gene': len(unmapped),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--feature-table', required=True, help='NCBI *_feature_table.txt(.gz)')
    p.add_argument('--protein-fasta', required=True, help='Protein FASTA to collapse')
    p.add_argument('--output', required=True, help='Destination collapsed FASTA')
    args = p.parse_args()

    stats = collapse(args.feature_table, args.protein_fasta, args.output)
    print(
        f"Wrote {args.output}: {stats['genes']} genes from {stats['input_records']} "
        f"input records ({stats['dropped_isoforms']} isoforms dropped, "
        f"{stats['unmapped_no_gene']} records had no gene mapping and were excluded)",
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
