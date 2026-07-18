#!/usr/bin/env python3
"""Extract candidate protein sequences from the given proteome FASTAs.

candidates.txt format: source_proteome_short::protein_id  (one per line)

--fastas is whatever set of proteome FASTAs the caller staged (ingroup for
the novelty-search direction, outgroup for the loss-search direction) — the
short-name lookup is built from every config sample, not filtered by group,
so the same script serves both without a --group flag: a candidate whose
short id has no matching file among --fastas is simply never found.
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
                    help='Proteome FASTA files to search for candidate sequences')
    ap.add_argument('--config', required=True,
                    help='Analysis CSV (used to map FASTA filenames to short IDs)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    # Build FASTA basename → short name from config. Not filtered by group —
    # which species actually get searched is determined by which FASTAs the
    # caller passed via --fastas, not by this lookup.
    samples = parse_config(args.config)
    basename_to_short = {s.protein: s.short for s in samples}

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
