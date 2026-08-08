#!/usr/bin/env python3
"""
Add gene_name, product_description, function_source, Best_Swissprot, and
Pfam_Names columns to presence_matrix.tsv.

Annotation priority per protein:
  1. Model organism gene names (via --modelorgs_config YAML)
  2. Pfam-A hmmsearch --tblout (all unique domain names per protein)
  3. SwissProt diamond blastp outfmt 6 with stitle (best hit)

Model organism lookup is configured via a YAML file (--modelorgs_config).
See lib/model_organisms.py for format details.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from model_organisms import ModelOrgAnnotator  # noqa: E402
from fasta import read_fasta  # noqa: E402


def parse_pfam_tblout(tblout_path):
    """Return {protein_id: [(name, accession, evalue), ...]} from hmmsearch --tblout.

    Column layout (hmmsearch format, target=sequence, query=HMM):
      parts[0] = protein ID, parts[2] = domain name, parts[3] = Pfam accession, parts[4] = e-value
    """
    hits: dict[str, list[tuple]] = {}
    seen: dict[str, set] = {}
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            protein_id, domain_name, domain_acc, evalue = parts[0], parts[2], parts[3], parts[4]
            entry = hits.setdefault(protein_id, [])
            seen_names = seen.setdefault(protein_id, set())
            if domain_name not in seen_names:
                seen_names.add(domain_name)
                entry.append((domain_name, domain_acc, evalue))
    return hits


def parse_swissprot_hits(tsv_path):
    """
    Return {query_id: description} from diamond blastp with
    --outfmt '6 qseqid sseqid stitle' (best hit per query).
    Strips the leading 'sp|ACCN|ID ' prefix from stitle.
    """
    hits = {}
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            qid, stitle = parts[0], parts[2]
            m = re.match(r'sp\|\S+\|\S+\s+(.*)', stitle)
            hits.setdefault(qid, m.group(1) if m else stitle)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.tsv from the pipeline')
    ap.add_argument('--modelorgs_config',
                    help='YAML file listing model organisms for gene name lookup '
                         '(see configs/modelorgs.yaml)')
    ap.add_argument('--launch_dir',
                    help='Project root directory; relative paths in the modelorgs YAML '
                         'are resolved from here (defaults to the YAML file\'s directory)')
    ap.add_argument('--pfam_hits',
                    help='hmmsearch --tblout output (candidates vs Pfam-A)')
    ap.add_argument('--swissprot_hits',
                    help='Diamond blastp outfmt "6 qseqid sseqid stitle" vs SwissProt')
    ap.add_argument('--candidates_fa',
                    help='FASTA file of candidate proteins; adds protein_sequence column when provided')
    ap.add_argument('--output', required=True,
                    help='Output TSV: presence_matrix.function.tsv')
    args = ap.parse_args()

    # --- build lookup tables ---
    annotator = None
    if args.modelorgs_config:
        annotator = ModelOrgAnnotator.from_yaml(args.modelorgs_config,
                                                launch_dir=args.launch_dir)

    pfam_hits: dict[str, list[tuple]] = {}
    if args.pfam_hits:
        pfam_hits = parse_pfam_tblout(args.pfam_hits)

    swissprot_hits: dict[str, str] = {}
    if args.swissprot_hits:
        swissprot_hits = parse_swissprot_hits(args.swissprot_hits)

    sequences: dict[str, str] = {}
    if args.candidates_fa:
        sequences = {rec_id: str(rec.seq) for rec_id, rec in read_fasta(args.candidates_fa).items()}

    # --- annotate ---
    with open(args.matrix) as fin, open(args.output, 'w', newline='') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        extra_cols = ['gene_name', 'product_description', 'function_source',
                      'Best_Swissprot', 'Pfam_Names', 'Pfam_Accessions', 'Pfam_Evalues']
        if args.candidates_fa:
            extra_cols.append('protein_sequence')
        out_fields = list(reader.fieldnames) + extra_cols
        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()

        for row in reader:
            pid = row['protein_id']
            source = row.get('source_proteome', '')
            gene_name = ''
            product = ''
            func_source = ''

            # --- Model organism lookup (via YAML config) ---
            if annotator:
                gene_name, product = annotator.annotate(pid, source)
                if product or gene_name:
                    func_source = f'ModelOrg_{source}'

            # --- Best Swissprot hit (always stored separately) ---
            best_swissprot = swissprot_hits.get(pid, '')

            # --- All Pfam domain names/accessions/evalues for this protein ---
            pfam_entries = pfam_hits.get(pid, [])
            pfam_names = ','.join(e[0] for e in pfam_entries)
            pfam_accessions = ','.join(e[1] for e in pfam_entries)
            pfam_evalues = ','.join(e[2] for e in pfam_entries)

            # --- fallback product_description: Pfam first, then SwissProt ---
            if not product and pfam_names:
                product = pfam_names
                func_source = 'Pfam'

            if not product and best_swissprot:
                product = best_swissprot
                func_source = 'SwissProt'

            row['gene_name'] = gene_name
            row['product_description'] = product
            row['function_source'] = func_source
            row['Best_Swissprot'] = best_swissprot
            row['Pfam_Names'] = pfam_names
            row['Pfam_Accessions'] = pfam_accessions
            row['Pfam_Evalues'] = pfam_evalues
            if args.candidates_fa:
                row['protein_sequence'] = sequences.get(pid, '')
            writer.writerow(row)


if __name__ == '__main__':
    main()
