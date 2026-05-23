#!/usr/bin/env python3
"""Extract candidate protein sequences from ingroup proteome FASTAs.

candidates.txt format: source_proteome_short::protein_id  (one per line)
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', required=True,
                    help='candidates.txt with lines like short::protein_id')
    ap.add_argument('--fastas', nargs='+', required=True,
                    help='Ingroup proteome FASTA files')
    ap.add_argument('--config', required=True,
                    help='Analysis CSV (used to map FASTA filenames to short IDs)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    # Build FASTA basename → short name from config
    samples = parse_config(args.config)
    basename_to_short = {s.protein: s.short for s in samples if s.group == 'IN'}

    # Parse candidates grouped by short name
    wanted: dict[str, set] = defaultdict(set)
    with open(args.candidates) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            short, prot_id = line.split('::', 1)
            wanted[short].add(prot_id)

    written = 0
    with open(args.output, 'w') as out_fh:
        for fa_path in args.fastas:
            short = basename_to_short.get(Path(fa_path).name)
            if short is None or short not in wanted:
                continue
            target_ids = wanted[short]
            for rec in SeqIO.parse(fa_path, 'fasta'):
                if rec.id in target_ids:
                    SeqIO.write(rec, out_fh, 'fasta')
                    written += 1

    print(f"Extracted {written} candidate sequences", file=sys.stderr)


if __name__ == '__main__':
    main()
