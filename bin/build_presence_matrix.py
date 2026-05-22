#!/usr/bin/env python3
"""Build a protein × proteome presence/absence matrix and emit candidate IDs.

A candidate protein (from an ingroup proteome) must be:
  - present (has a hit) in >= --ingroup-min-frac of all ingroup proteomes
  - absent (no hit) from every outgroup proteome
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config, short_to_group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', nargs='+', required=True,
                    help='Parsed hit TSV files (one per query-target pair)')
    ap.add_argument('--config',  required=True, help='Analysis description CSV')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75,
                    dest='ingroup_min_frac')
    ap.add_argument('--output-matrix',     required=True)
    ap.add_argument('--output-candidates', required=True)
    args = ap.parse_args()

    samples      = parse_config(args.config)
    group_map    = short_to_group(samples)
    ingroup_ids  = {s.short for s in samples if s.group == 'IN'}
    outgroup_ids = {s.short for s in samples if s.group == 'OUT'}
    all_ids      = ingroup_ids | outgroup_ids

    # protein_key → (query_proteome, protein_id)  uniquely identifies a protein
    # Presence recorded as a set of target proteomes where a hit was found.
    # The query protein is also present in its own proteome by definition.
    presence: dict[tuple, set] = defaultdict(set)

    for tsv_path in args.hits:
        df = pd.read_csv(tsv_path, sep='\t')
        if df.empty:
            continue
        for _, row in df.iterrows():
            qp = row['query_proteome']
            tp = row['target_proteome']
            if qp not in ingroup_ids:
                continue   # only track ingroup query proteins
            key = (qp, row['query_id'])
            presence[key].add(tp)

    # Build the full matrix
    rows = []
    for (qp, pid), hit_proteomes in presence.items():
        # The protein is present in its own proteome plus wherever it has hits
        all_present = hit_proteomes | {qp}
        row = {'protein_id': pid, 'source_proteome': qp}
        for pid2 in sorted(all_ids):
            row[pid2] = int(pid2 in all_present)
        rows.append(row)

    matrix = pd.DataFrame(rows).fillna(0)
    matrix.to_csv(args.output_matrix, sep='\t', index=False)

    # Apply thresholds
    n_ingroup = len(ingroup_ids)
    candidates = []
    for _, row in matrix.iterrows():
        ingroup_count  = sum(int(row.get(s, 0)) for s in ingroup_ids)
        outgroup_count = sum(int(row.get(s, 0)) for s in outgroup_ids)
        frac = ingroup_count / n_ingroup
        if frac >= args.ingroup_min_frac and outgroup_count == 0:
            candidates.append(f"{row['source_proteome']}::{row['protein_id']}")

    with open(args.output_candidates, 'w') as fh:
        fh.write('\n'.join(candidates) + '\n')

    print(f"Candidates: {len(candidates)} / {len(matrix)} ingroup proteins pass thresholds",
          file=sys.stderr)


if __name__ == '__main__':
    main()
