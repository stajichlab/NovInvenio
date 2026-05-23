#!/usr/bin/env python3
"""Build a protein × proteome presence/absence matrix and emit candidate IDs.

Presence scoring uses per-gene paralog e-value cutoffs (from self-vs-self
searches).  A cross-species hit is counted only when its e-value is strictly
less than the query protein's paralog cutoff.  Proteins absent from the
paralog cutoff files fall back to --default-evalue (1e-5).

A candidate protein (from an ingroup proteome) must be:
  - present in >= --ingroup-min-frac of all ingroup proteomes
  - absent from every outgroup proteome
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config, short_to_group

DEFAULT_EVALUE = 1e-5


def load_paralog_cutoffs(cutoff_files, default_evalue):
    """Return dict protein_id -> evalue threshold from all paralog_cutoffs.tsv files."""
    cutoffs = {}
    for path in cutoff_files:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        for _, row in df.iterrows():
            cutoffs[row['protein_ID']] = float(row['evalue'])
    return cutoffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', nargs='+', required=True,
                    help='Parsed hit TSV files (one per query-target pair)')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[],
                    dest='paralog_cutoffs',
                    help='Per-species paralog_cutoffs.tsv files from self-vs-self search')
    ap.add_argument('--config',  required=True, help='Analysis description CSV')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75,
                    dest='ingroup_min_frac')
    ap.add_argument('--default-evalue', type=float, default=DEFAULT_EVALUE,
                    dest='default_evalue',
                    help='Fallback e-value cutoff for proteins with no detectable paralog')
    ap.add_argument('--output-matrix',     required=True)
    ap.add_argument('--output-candidates', required=True)
    args = ap.parse_args()

    samples      = parse_config(args.config)
    ingroup_ids  = {s.short for s in samples if s.group == 'IN'}
    outgroup_ids = {s.short for s in samples if s.group == 'OUT'}
    all_ids      = ingroup_ids | outgroup_ids

    paralog_cutoffs = load_paralog_cutoffs(args.paralog_cutoffs, args.default_evalue)

    # protein_key (query_proteome, protein_id) -> set of target proteomes with qualifying hits
    presence: dict[tuple, set] = defaultdict(set)

    for tsv_path in args.hits:
        df = pd.read_csv(tsv_path, sep='\t')
        if df.empty:
            continue
        for _, row in df.iterrows():
            qp = row['query_proteome']
            if qp not in ingroup_ids:
                continue
            qid     = row['query_id']
            cutoff  = paralog_cutoffs.get(qid, args.default_evalue)
            if float(row['evalue']) < cutoff:
                presence[(qp, qid)].add(row['target_proteome'])

    # Build the full matrix
    rows = []
    for (qp, pid), hit_proteomes in presence.items():
        all_present = hit_proteomes | {qp}
        row = {'protein_id': pid, 'source_proteome': qp}
        for sp in sorted(all_ids):
            row[sp] = int(sp in all_present)
        rows.append(row)

    matrix = pd.DataFrame(rows).fillna(0)
    matrix.to_csv(args.output_matrix, sep='\t', index=False)

    n_ingroup  = len(ingroup_ids)
    candidates = []
    for _, row in matrix.iterrows():
        ingroup_count  = sum(int(row.get(s, 0)) for s in ingroup_ids)
        outgroup_count = sum(int(row.get(s, 0)) for s in outgroup_ids)
        if ingroup_count / n_ingroup >= args.ingroup_min_frac and outgroup_count == 0:
            candidates.append(f"{row['source_proteome']}::{row['protein_id']}")

    with open(args.output_candidates, 'w') as fh:
        fh.write('\n'.join(candidates) + '\n')

    print(f"Candidates: {len(candidates)} / {len(matrix)} ingroup proteins pass thresholds",
          file=sys.stderr)


if __name__ == '__main__':
    main()
