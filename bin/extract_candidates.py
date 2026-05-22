#!/usr/bin/env python3
"""Extract candidate protein sequences from ingroup proteome FASTAs.

candidates.txt format: source_proteome::protein_id  (one per line)
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', required=True,
                    help='candidates.txt with lines like proteome::protein_id')
    ap.add_argument('--fastas', nargs='+', required=True,
                    help='Ingroup proteome FASTA files')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    # Parse candidates grouped by source proteome
    wanted: dict[str, set] = defaultdict(set)
    with open(args.candidates) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            proteome, prot_id = line.split('::', 1)
            wanted[proteome].add(prot_id)

    written = 0
    with open(args.output, 'w') as out_fh:
        for fa_path in args.fastas:
            # Derive proteome short-name from filename stem
            stem = Path(fa_path).stem
            if stem not in wanted:
                continue
            target_ids = wanted[stem]
            for rec in SeqIO.parse(fa_path, 'fasta'):
                if rec.id in target_ids:
                    SeqIO.write(rec, out_fh, 'fasta')
                    written += 1

    print(f"Extracted {written} candidate sequences", file=sys.stderr)


if __name__ == '__main__':
    main()
