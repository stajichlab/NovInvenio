#!/usr/bin/env python3
"""Normalise raw search output (phmmer/diamond/blast) to a common 6-column TSV.

Output columns: query_id, target_id, evalue, bitscore, query_proteome, target_proteome
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from hits import PARSERS, filter_hits, open_input


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',  required=True, help='Raw search output file')
    ap.add_argument('--format', required=True, choices=PARSERS.keys())
    ap.add_argument('--output', required=True, help='Output TSV path')
    ap.add_argument('--evalue', type=float, default=1e-15)
    ap.add_argument('--bitscore', type=float, default=None)
    ap.add_argument('--query-proteome',  default='', dest='query_proteome')
    ap.add_argument('--target-proteome', default='', dest='target_proteome')
    args = ap.parse_args()

    parser = PARSERS[args.format]
    header = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'

    with open_input(args.input) as fh_in, open(args.output, 'w') as fh_out:
        fh_out.write(header)
        for hit in filter_hits(parser(fh_in), args.evalue, args.bitscore):
            hit.query_proteome  = args.query_proteome
            hit.target_proteome = args.target_proteome
            fh_out.write(
                f'{hit.query_id}\t{hit.target_id}\t{hit.evalue}\t'
                f'{hit.bitscore}\t{hit.query_proteome}\t{hit.target_proteome}\n'
            )


if __name__ == '__main__':
    main()
