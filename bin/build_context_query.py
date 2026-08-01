#!/usr/bin/env python3
"""Extend candidates.txt with each candidate's within-genome paralog (issue #48).

The NEAR_INGROUP/BROAD_OUTGROUP context search (bin/context_presence.py) reuses the
same paralog-competition filter as bin/build_presence_matrix.py: a candidate's hit
against a context proteome is disqualified if the candidate's own paralog scores
better on the same target. That comparison needs the *paralog's* hit data too, so
the paralog must be searched alongside the candidate -- otherwise there is nothing
to compare against. This script just extends the query list; the actual search
(workflows/context_search.nf) treats every line the same way.

Paralogs are looked up from paralog_cutoffs.tsv (protein_ID -> paralog_protein_ID),
the same self-vs-self search output build_presence_matrix.py already reads. A
candidate with no detected paralog contributes nothing extra. Output format matches
candidates.txt: source_proteome::protein_id, one per line, deduplicated.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def load_paralog_of(cutoff_files):
    """Return protein_id -> paralog_protein_id from paralog_cutoffs.tsv files."""
    paralog_of = {}
    for path in cutoff_files:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        paralog_of.update(dict(zip(df['protein_ID'], df['paralog_protein_ID'])))
    return paralog_of


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--candidates', required=True,
                    help='candidates.txt (source_proteome::protein_id, one per line)')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[], dest='paralog_cutoffs',
                    help='paralog_cutoffs.tsv files from the self-vs-self search')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    paralog_of = load_paralog_of(args.paralog_cutoffs)

    entries: set[str] = set()
    with open(args.candidates) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            entries.add(line)
            short, _, protein_id = line.partition('::')
            # A protein's within-genome paralog is always in the same proteome
            # (self-vs-self search never crosses proteomes).
            paralog_id = paralog_of.get(protein_id)
            if paralog_id and not pd.isna(paralog_id):
                entries.add(f'{short}::{paralog_id}')

    with open(args.output, 'w') as out:
        if entries:
            out.write('\n'.join(sorted(entries)) + '\n')

    print(f"Context query: {len(entries)} proteins (candidates + their paralogs)",
          file=sys.stderr)


if __name__ == '__main__':
    main()
