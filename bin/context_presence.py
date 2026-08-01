#!/usr/bin/env python3
"""Score novelty candidates against NEAR_INGROUP/BROAD_OUTGROUP as report-only context
(issue #48).

These proteomes never determine novelty -- candidates.txt is already fixed by the
strict IN/OUT search (bin/build_presence_matrix.py) by the time this runs. This just
answers "does this candidate also show up in a close relative, or a distant lineage?"
for the report's benefit, using the *same* paralog-aware filtering:

  1. Paralog-cutoff filter: hit e-value < the candidate's own within-proteome paralog
     e-value (fallback --default-evalue).
  2. Paralog-competition filter (--paralog-competition-scope, mirrors
     build_presence_matrix.py): disqualify a hit the candidate's paralog out-scores.

Filter 2 needs the paralog's own hits against the same context proteomes, which is
why --hits also includes hits from bin/build_context_query.py's extended query list
(candidates + their paralogs) -- only rows whose query_id is itself a candidate are
ever scored/emitted here; paralog rows are used purely for the filter-2 lookup.

Emits two TSVs, same shape as build_presence_matrix.py's matrix/evalues sidecar but
restricted to candidate protein rows and NEAR_INGROUP/BROAD_OUTGROUP columns:
  --output-matrix:  protein_id, source_proteome, <context proteome shorts> (0/1)
  --output-evalues: same shape, qualifying hit e-value (or '') instead of 0/1
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config

DEFAULT_EVALUE = 1e-5


def load_candidates(path):
    """Return {protein_id: source_proteome} from candidates.txt (source::protein)."""
    candidates = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            short, _, protein_id = line.partition('::')
            candidates[protein_id] = short
    return candidates


def load_paralog_info(cutoff_files):
    cutoffs, paralog_of = {}, {}
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hits', nargs='+', required=True,
                    help='Parsed hit TSVs: context query (candidates + paralogs) vs '
                         'NEAR_INGROUP/BROAD_OUTGROUP proteomes')
    ap.add_argument('--candidates', required=True,
                    help='candidates.txt -- only these protein_ids are scored/emitted')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[], dest='paralog_cutoffs',
                    help='paralog_cutoffs.tsv files from the ingroup self-vs-self search')
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--paralog-competition-scope', choices=['proteome', 'target'],
                    default='proteome', dest='paralog_competition_scope',
                    help='Same semantics as build_presence_matrix.py (default: proteome)')
    ap.add_argument('--default-evalue', type=float, default=DEFAULT_EVALUE,
                    dest='default_evalue')
    ap.add_argument('--output-matrix', required=True, dest='output_matrix')
    ap.add_argument('--output-evalues', required=True, dest='output_evalues')
    args = ap.parse_args()

    samples = parse_config(args.config)
    context_ids = sorted(
        s.short for s in samples if s.group in ('NEAR_INGROUP', 'BROAD_OUTGROUP')
    )

    candidates = load_candidates(args.candidates)
    paralog_cutoffs, paralog_of = load_paralog_info(args.paralog_cutoffs)
    hits = load_hits(args.hits)

    if hits.empty or not context_ids:
        _write(args.output_matrix, candidates, context_ids, {}, as_presence=True)
        _write(args.output_evalues, candidates, context_ids, {}, as_presence=False)
        print("Context presence: 0 hits scored (no hits or no context proteomes)",
              file=sys.stderr)
        return

    # Best evalue per (query_proteome, query_id, target_key) for filter 2 -- covers
    # BOTH candidates and their paralogs, since either might be the "paralog side"
    # of a comparison.
    if args.paralog_competition_scope == 'proteome':
        best_ev = (hits.groupby(['query_proteome', 'query_id', 'target_proteome'])
                       ['evalue'].min().to_dict())
    else:
        best_ev = (hits.groupby(['query_proteome', 'query_id', 'target_id'])
                       ['evalue'].min().to_dict())

    # Only candidate rows are ever emitted; paralog rows exist in `hits` solely to
    # populate best_ev above.
    cand = hits[hits['query_id'].isin(candidates)].copy()

    if not cand.empty:
        cutoff = cand['query_id'].map(paralog_cutoffs).fillna(args.default_evalue)
        cand = cand[cand['evalue'] < cutoff]

    if not cand.empty:
        paralog_ids = cand['query_id'].map(paralog_of)
        target_key = (cand['target_proteome'] if args.paralog_competition_scope == 'proteome'
                      else cand['target_id'])
        paralog_ev = [
            best_ev.get((qp, pid, tk)) if pid is not None and not pd.isna(pid) else None
            for qp, pid, tk in zip(cand['query_proteome'], paralog_ids, target_key)
        ]
        paralog_ev = pd.Series(paralog_ev, index=cand.index, dtype='float64')
        disqualified = paralog_ev.notna() & (paralog_ev < cand['evalue'])
        cand = cand[~disqualified]

    # Best (lowest) qualifying evalue per (protein_id, target_proteome).
    best_hit = (cand.groupby(['query_id', 'target_proteome'])['evalue'].min().to_dict()
               if not cand.empty else {})

    _write(args.output_matrix, candidates, context_ids, best_hit, as_presence=True)
    _write(args.output_evalues, candidates, context_ids, best_hit, as_presence=False)

    n_scored = sum(1 for (_, tp) in best_hit if tp in context_ids)
    print(f"Context presence: {len(candidates)} candidates x {len(context_ids)} context "
          f"proteomes, {n_scored} qualifying hits", file=sys.stderr)


def _write(path, candidates, context_ids, best_hit, as_presence):
    columns = ['protein_id', 'source_proteome'] + context_ids
    with open(path, 'w') as out:
        out.write('\t'.join(columns) + '\n')
        for pid, source in sorted(candidates.items()):
            vals = []
            for cid in context_ids:
                ev = best_hit.get((pid, cid))
                if as_presence:
                    vals.append('1' if ev is not None else '0')
                else:
                    vals.append(str(ev) if ev is not None else '')
            out.write('\t'.join([pid, source] + vals) + '\n')


if __name__ == '__main__':
    main()
