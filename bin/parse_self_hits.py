#!/usr/bin/env python3
"""
Extract each protein's within-proteome paralog from a self-vs-self search.

For each query protein the rank-2 hit (first hit where target_id != query_id)
is the best within-proteome paralog.  build_presence_matrix.py (and the
novelty_discovery/context_search equivalents) use paralog_protein_ID for the
paralog-competition filter: a cross-species hit is disqualified if the query's
paralog scores better against the same target.  The evalue column here is
informational only (report-only, e.g. for debugging) and is NOT used as a
presence-calling significance cutoff -- see lib/singleton_presence.py's module
docstring for why a per-query self-search-derived cutoff was tried and
dropped.  Proteins with no detectable within-proteome paralog are omitted.

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
