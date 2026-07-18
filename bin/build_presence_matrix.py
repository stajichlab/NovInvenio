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

A candidate protein (from a --query-group proteome, default IN) must be:
  - present in >= --ingroup-min-frac of all --query-group proteomes
  - present in <= --other-max-frac of the other group's proteomes (default 0.0,
    i.e. absent from every other-group proteome)

--query-group OUT runs the same logic in the opposite direction — outgroup
proteomes as query, looking for genes conserved in the outgroup but absent
from the ingroup (candidate lineage-specific losses). --ingroup-min-frac is
reused as the query group's own presence threshold in both directions (pass
--ingroup-min-frac params.outgroup_min_frac from the loss-search workflow).
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
        cutoffs.update(dict(zip(df['protein_ID'], df['evalue'].astype(float))))
        paralog_of.update(dict(zip(df['protein_ID'], df['paralog_protein_ID'])))
    return cutoffs, paralog_of


HIT_COLUMNS = ['query_id', 'target_id', 'evalue', 'bitscore',
               'query_proteome', 'target_proteome']


def load_hits(hit_files):
    """Concatenate all parsed pairwise-hit TSVs into a single DataFrame."""
    frames = []
    for tsv_path in hit_files:
        df = pd.read_csv(tsv_path, sep='\t')
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=HIT_COLUMNS)
    hits = pd.concat(frames, ignore_index=True)
    hits['evalue'] = hits['evalue'].astype(float)
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', nargs='+', required=True,
                    help='Parsed hit TSV files (one per query-target pair)')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[],
                    dest='paralog_cutoffs',
                    help='Per-species paralog_cutoffs.tsv files from self-vs-self search')
    ap.add_argument('--config',  required=True, help='Analysis description CSV')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75,
                    dest='ingroup_min_frac',
                    help='Presence threshold within --query-group (fraction of that '
                         'group\'s proteomes that must contain a hit)')
    ap.add_argument('--query-group', choices=['IN', 'OUT'], default='IN',
                    dest='query_group',
                    help='Which config group supplies the query proteomes candidates are '
                         'sourced from (default: IN, the novelty-search direction). OUT '
                         'runs the loss-search direction: outgroup query, ingroup must be 0.')
    ap.add_argument('--other-max-frac', type=float, default=0.0,
                    dest='other_max_frac',
                    help='Max fraction of the *other* group a candidate may still be '
                         'present in (default: 0.0 = must be absent from every other-group '
                         'proteome). Relax it in the loss direction to allow candidates '
                         'that survive in a small fraction of the ingroup — genes nearly, '
                         'but not entirely, lost.')
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
    query_ids    = ingroup_ids if args.query_group == 'IN' else outgroup_ids
    other_ids    = outgroup_ids if args.query_group == 'IN' else ingroup_ids

    paralog_cutoffs, paralog_of = load_paralog_info(args.paralog_cutoffs)

    hits = load_hits(args.hits)

    # Best evalue per (source_proteome, protein_id, target_proteome).
    # Used by filter 2 to check how well a paralog hits the same target.
    if hits.empty:
        best_ev = {}
    else:
        best_ev = (hits.groupby(['query_proteome', 'query_id', 'target_proteome'])
                       ['evalue'].min().to_dict())

    # Restrict to --query-group queries, then apply both paralog-aware filters vectorised.
    ing = hits[hits['query_proteome'].isin(query_ids)].copy()

    if not ing.empty:
        # Filter 1: hit evalue must beat the query's paralog cutoff (fallback default).
        cutoff = ing['query_id'].map(paralog_cutoffs).fillna(args.default_evalue)
        ing = ing[ing['evalue'] < cutoff]

    if not ing.empty:
        # Filter 2: disqualify if the query's paralog hits this same target better.
        paralog_ids = ing['query_id'].map(paralog_of)
        paralog_ev = [
            best_ev.get((qp, pid, tp)) if pid is not None and not pd.isna(pid) else None
            for qp, pid, tp in zip(ing['query_proteome'], paralog_ids, ing['target_proteome'])
        ]
        paralog_ev = pd.Series(paralog_ev, index=ing.index, dtype='float64')
        disqualified = paralog_ev.notna() & (paralog_ev < ing['evalue'])
        ing = ing[~disqualified]

    # protein_key (query_proteome, protein_id) -> set of target proteomes with qualifying hits
    presence: dict[tuple, set] = defaultdict(set)
    for qp, qid, tp in zip(ing['query_proteome'], ing['query_id'], ing['target_proteome']):
        presence[(qp, qid)].add(tp)

    # Build the full matrix (always emit the id columns + one column per proteome,
    # so an empty result still writes a well-formed header).
    sorted_ids = sorted(all_ids)
    columns = ['protein_id', 'source_proteome'] + sorted_ids
    rows = []
    for (qp, pid), hit_proteomes in presence.items():
        all_present = hit_proteomes | {qp}
        row = {'protein_id': pid, 'source_proteome': qp}
        for sp in sorted_ids:
            row[sp] = int(sp in all_present)
        rows.append(row)

    matrix = pd.DataFrame(rows, columns=columns)
    matrix.to_csv(args.output_matrix, sep='\t', index=False)

    n_query     = len(query_ids)
    query_cols  = sorted(query_ids)
    other_cols  = sorted(other_ids)
    if matrix.empty:
        candidates = []
    else:
        query_count = matrix[query_cols].sum(axis=1)
        other_count = matrix[other_cols].sum(axis=1) if other_cols else 0
        other_frac = (other_count / len(other_cols)) if other_cols else 0
        keep = (query_count / n_query >= args.ingroup_min_frac) & (other_frac <= args.other_max_frac)
        kept = matrix[keep]
        candidates = (kept['source_proteome'] + '::' + kept['protein_id']).tolist()

    with open(args.output_candidates, 'w') as fh:
        if candidates:
            fh.write('\n'.join(candidates) + '\n')

    direction = 'ingroup' if args.query_group == 'IN' else 'outgroup'
    print(f"Candidates: {len(candidates)} / {len(matrix)} {direction} proteins pass thresholds",
          file=sys.stderr)


if __name__ == '__main__':
    main()
