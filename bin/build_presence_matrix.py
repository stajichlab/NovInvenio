#!/usr/bin/env python3
"""Build a protein × proteome presence/absence matrix and emit candidate IDs.

Presence scoring uses two paralog-aware filters (both must pass):

  1. Paralog-cutoff filter: hit e-value < query protein's within-proteome
     paralog e-value (from self-vs-self search).  Proteins with no detectable
     paralog fall back to --default-evalue (1e-5).

  2. Paralog-competition filter: if the query's paralog hits the same target
     proteome with a *better* (lower) e-value than the query does, the hit is
     disqualified.  The logic is that such a hit is better explained by the
     conserved domain shared with the paralog than by the query protein itself.

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


def load_paralog_info(cutoff_files):
    """Return (cutoffs, paralog_of) from all paralog_cutoffs.tsv files.

    cutoffs:    protein_id -> float evalue threshold
    paralog_of: protein_id -> paralog_protein_id
    """
    cutoffs    = {}
    paralog_of = {}
    for path in cutoff_files:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        for _, row in df.iterrows():
            pid = row['protein_ID']
            cutoffs[pid]    = float(row['evalue'])
            paralog_of[pid] = row['paralog_protein_ID']
    return cutoffs, paralog_of


def build_hit_lookup(hit_files):
    """Build (source_proteome, protein_id, target_proteome) -> best evalue.

    Covers every ingroup query across all pairwise TSVs; used to check how
    well a paralog hits the same target as its query protein.
    """
    lookup: dict[tuple, float] = {}
    for tsv_path in hit_files:
        df = pd.read_csv(tsv_path, sep='\t')
        if df.empty:
            continue
        for _, row in df.iterrows():
            key = (row['query_proteome'], row['query_id'], row['target_proteome'])
            ev  = float(row['evalue'])
            if ev < lookup.get(key, float('inf')):
                lookup[key] = ev
    return lookup


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

    paralog_cutoffs, paralog_of = load_paralog_info(args.paralog_cutoffs)

    # Pass 1: build lookup of best evalue per (source_proteome, protein_id, target_proteome).
    # Needed by filter 2 to check how well a paralog hits the same target.
    hit_lookup = build_hit_lookup(args.hits)

    # Pass 2: apply both filters and accumulate presence.
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
            qid = row['query_id']
            tp  = row['target_proteome']
            ev  = float(row['evalue'])

            # Filter 1: hit must be better than the paralog cutoff.
            if ev >= paralog_cutoffs.get(qid, args.default_evalue):
                continue

            # Filter 2: if the paralog hits this same target better, disqualify.
            paralog_id = paralog_of.get(qid)
            if paralog_id is not None:
                paralog_ev = hit_lookup.get((qp, paralog_id, tp))
                if paralog_ev is not None and paralog_ev < ev:
                    continue

            presence[(qp, qid)].add(tp)

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
