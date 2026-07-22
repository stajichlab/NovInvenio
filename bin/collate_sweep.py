#!/usr/bin/env python3
"""Collate a family-parameter sweep into a scored table and pick the knee (ADR-0002 Q8).

The sweep (issue #6) runs the family-profile pathway across a grid of family-definition
parameters — mmseqs `--min-seq-id` × `-c` (cov) × hmmsearch presence E-value — and scores
each grid point on: number of families, number of novelty candidates, BUSCO single-copy
recovery (bin/busco_family_recovery.py), and control recall / false-novelty rate
(bin/score_controls.py). This script reads the per-point metrics the harness gathered
(`bin/run_param_sweep.sh`) and turns the guessed default into an evidence-backed one.

Input: a wide TSV (`--metrics`) with one row per grid point and columns:
  min_seq_id, cov, hmm_evalue, n_families, n_novelties, busco_recovery, recall, fp_rate,
  tblastn_removed
(any missing numeric cell is treated as the worst value for ranking).

Knee selection is a transparent, documented rule — not a black box:
  1. Admissible = busco_recovery >= --min-busco-recovery AND fp_rate <= --max-fp
     (quality gates: a family definition that shreds single-copy orthologs or floods false
     novelties is disqualified regardless of recall).
  2. Rank admissible points by composite = recall + busco_recovery - fp_rate (all ~[0,1]).
  3. The knee: among near-top composites, looser identity thresholds keep inflating
     n_novelties without buying recall/recovery, so ties break to FEWER novelties (then
     the shipped default params). Pick that point's parameters as the recommended default.
If no point is admissible, fall back to the best composite overall and flag it.
"""
import argparse
import csv
import sys

PARAM_COLS = ('min_seq_id', 'cov', 'hmm_evalue')
# metric -> (higher_is_better, worst_value_when_missing)
METRICS = {
    'n_families': (None, 0),
    'n_novelties': (None, 0),
    'busco_recovery': (True, 0.0),
    'recall': (True, 0.0),
    'fp_rate': (False, 1.0),
    'tblastn_removed': (None, 0),
}
# Shipped default (nextflow.config) — the tie-break anchor when scores are equal.
DEFAULT_PARAMS = {'min_seq_id': 0.3, 'cov': 0.8, 'hmm_evalue': 1e-3}


def _num(value, fallback):
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def read_metrics(path):
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = []
        for raw in reader:
            point = {}
            for c in PARAM_COLS:
                point[c] = _num(raw.get(c), None)
            for m, (_hib, worst) in METRICS.items():
                point[m] = _num(raw.get(m), worst)
            rows.append(point)
    return rows


def score_points(rows, min_busco_recovery, max_fp):
    for p in rows:
        p['admissible'] = int(p['busco_recovery'] >= min_busco_recovery
                              and p['fp_rate'] <= max_fp)
        p['composite'] = round(p['recall'] + p['busco_recovery'] - p['fp_rate'], 6)
    return rows


def _params_match_default(p):
    return all(p.get(c) == DEFAULT_PARAMS[c] for c in PARAM_COLS)


def select_knee(rows):
    """Return (chosen_point, used_fallback). Ranking rule documented in the module docstring."""
    if not rows:
        return None, False
    admissible = [p for p in rows if p['admissible']]
    pool = admissible or rows
    used_fallback = not admissible
    # Best composite; tie-break: fewer novelties (the knee), then shipped-default params,
    # then higher recall — all deterministic.
    chosen = min(
        pool,
        key=lambda p: (
            -p['composite'],
            p['n_novelties'],
            0 if _params_match_default(p) else 1,
            -p['recall'],
        ),
    )
    return chosen, used_fallback


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--metrics', required=True, help='wide per-grid-point metrics TSV')
    ap.add_argument('--min-busco-recovery', type=float, default=0.9,
                    dest='min_busco_recovery', help='admissibility gate (default 0.9)')
    ap.add_argument('--max-fp', type=float, default=0.05, dest='max_fp',
                    help='max control false-novelty rate to be admissible (default 0.05)')
    ap.add_argument('--output', required=True, help='scored sweep table TSV (sorted best-first)')
    args = ap.parse_args()

    rows = read_metrics(args.metrics)
    if not rows:
        sys.exit(f'no grid points in {args.metrics}')
    score_points(rows, args.min_busco_recovery, args.max_fp)
    chosen, used_fallback = select_knee(rows)

    # Sort best-first: admissible first, then composite, then fewer novelties.
    ordered = sorted(rows, key=lambda p: (-p['admissible'], -p['composite'], p['n_novelties']))
    for p in ordered:
        p['chosen'] = int(p is chosen)

    cols = list(PARAM_COLS) + list(METRICS) + ['composite', 'admissible', 'chosen']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for p in ordered:
            w.writerow({c: p.get(c, '') for c in cols})

    if used_fallback:
        print('WARNING: no grid point met the admissibility gates '
              f'(busco_recovery >= {args.min_busco_recovery}, fp_rate <= {args.max_fp}); '
              'reporting the best composite overall.', file=sys.stderr)
    print(f"Recommended default: --family_min_seq_id {chosen['min_seq_id']} "
          f"--family_cov {chosen['cov']} --hmm_presence_evalue {chosen['hmm_evalue']:g} "
          f"(composite {chosen['composite']:.3f}, recall {chosen['recall']:.3f}, "
          f"busco {chosen['busco_recovery']:.3f}, fp {chosen['fp_rate']:.3f}, "
          f"novelties {int(chosen['n_novelties'])})", file=sys.stderr)


if __name__ == '__main__':
    main()
