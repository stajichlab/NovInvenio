#!/usr/bin/env python3
"""
Extract per-gene paralog e-value cutoffs from a self-vs-self search.

For each query protein the rank-2 hit (first hit where target_id != query_id)
is the best within-proteome paralog.  Its e-value becomes the cutoff used by
build_presence_matrix.py: a cross-species hit is only counted as a homolog if
its e-value is strictly less than this cutoff.  Proteins with no detectable
within-proteome paralog are omitted; build_presence_matrix.py falls back to
the global default (1e-5) for those.

Output TSV columns (one row per protein that has a detectable paralog):
  protein_ID, paralog_protein_ID, bitscore, evalue
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from hits import PARSERS, open_input


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',    required=True, help='Raw self-search output (.gz ok)')
    ap.add_argument('--format',   required=True, choices=PARSERS.keys())
    ap.add_argument('--proteome', required=True, help='Proteome short ID (unused, kept for interface)')
    ap.add_argument('--output',   required=True)
    args = ap.parse_args()

    parser = PARSERS[args.format]

    rank2_hit = {}   # query_id -> best (lowest-evalue) Hit where target_id != query_id

    with open_input(args.input) as fh:
        for hit in parser(fh):
            if hit.target_id == hit.query_id:
                continue
            prev = rank2_hit.get(hit.query_id)
            if prev is None or float(hit.evalue) < float(prev.evalue):
                rank2_hit[hit.query_id] = hit

    with open(args.output, 'w') as out:
        out.write('protein_ID\tparalog_protein_ID\tbitscore\tevalue\n')
        for qid in sorted(rank2_hit):
            r2 = rank2_hit[qid]
            out.write(f'{qid}\t{r2.target_id}\t{r2.bitscore:.1f}\t{r2.evalue:.2e}\n')


if __name__ == '__main__':
    main()
