#!/usr/bin/env python3
"""Extract a protein_id -> (gene_name, description) lookup from UniProt-style
FASTA headers, so the report can show a real name/description for a protein ID
instead of just the bare accession -- used both by build_alignment_shards.py
(the TBLASTN popup's query-protein name) and by lib/report_data.py (resolving
a pairwise search hit's target protein name).

UniProt FASTA headers look like:
  >tr|A7UWL5|A7UWL5_NEUCR NRS/ER OS=Neurospora crassa (strain ...) OX=367110 GN=NCU10683 PE=4 SV=1
  >sp|P12345|NAME_ORG Some protein description OS=... OX=... GN=... PE=... SV=...

The first whitespace-delimited token is the protein_id (matches the ID diamond/
phmmer/blast report as qseqid/sseqid throughout this pipeline). The description
is everything between that token and the first " OS=" (UniProt's organism-name
field always follows the free-text description, never precedes it). GN= is
optional -- many UniProt entries (mostly TrEMBL) have no gene_name.

A header with no " OS=" segment (a non-UniProt-style FASTA) still gets a
best-effort description: the whole remainder of the header line.

Usage:
  extract_protein_descriptions.py --pep species1.pep.fa species2.pep.fa ... \
      --output descriptions.tsv
"""
import argparse
import gzip
import re
import sys

_OS_SPLIT_RE = re.compile(r'\s+OS=')
_GN_RE = re.compile(r'\bGN=(\S+)')


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def parse_fasta_headers(fh):
    """Yield (protein_id, gene_name, description) for every '>' header line."""
    for line in fh:
        if not line.startswith('>'):
            continue
        header = line[1:].rstrip('\n')
        parts = header.split(None, 1)
        pid = parts[0]
        rest = parts[1] if len(parts) > 1 else ''
        description = _OS_SPLIT_RE.split(rest, 1)[0].strip()
        m = _GN_RE.search(rest)
        gene_name = m.group(1) if m else ''
        yield pid, gene_name, description


def build_descriptions(pep_paths):
    """{protein_id: (gene_name, description)} -- first occurrence wins on a
    duplicate ID (shouldn't normally happen across distinct proteomes, since
    IDs are already globally unique UniProt accessions)."""
    descriptions = {}
    for path in pep_paths:
        with _open(path) as fh:
            for pid, gene_name, description in parse_fasta_headers(fh):
                if pid not in descriptions:
                    descriptions[pid] = (gene_name, description)
    return descriptions


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pep', nargs='+', required=True,
                    help='Peptide FASTA file(s) (plain or .gz) to scan for headers')
    ap.add_argument('--output', required=True, help='Output TSV: protein_id\\tgene_name\\tdescription')
    args = ap.parse_args()

    descriptions = build_descriptions(args.pep)

    with open(args.output, 'w') as fh:
        fh.write('protein_id\tgene_name\tdescription\n')
        for pid in sorted(descriptions):
            gene_name, description = descriptions[pid]
            fh.write(f'{pid}\t{gene_name}\t{description}\n')

    print(f'Wrote {len(descriptions)} protein description(s) to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
