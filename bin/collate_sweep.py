#!/usr/bin/env python3
"""Collate a family-parameter sweep into a scored table and pick the knee (ADR-0002 Q8;
extended 2026-09-03 to cover the presence-calling coverage floor).

The sweep (issue #6) runs the family-profile pathway across a grid of family-definition
parameters — mmseqs `--min-seq-id` × `-c` (cov) × hmmsearch presence E-value/coverage
fraction/min-residues — and scores each grid point on: number of families, number of
novelty candidates, BUSCO single-copy clustering recovery (bin/busco_family_recovery.py),
BUSCO PRESENCE recovery (bin/busco_presence_recovery.py -- the metric actually sensitive
to hmm_cov/hmm_residues; busco_recovery is not, see that script's docstring), and control
recall / false-novelty rate (bin/score_controls.py). This script reads the per-point
metrics the harness gathered (`bin/run_param_sweep.sh`) and turns the guessed default into
an evidence-backed one.

Input: a wide TSV (`--metrics`) with one row per grid point and columns:
  min_seq_id, cov, hmm_evalue, hmm_cov, hmm_residues, n_families, n_novelties,
  busco_recovery, presence_recovery, recall, fp_rate, tblastn_removed
(any missing numeric cell is treated as the worst value for ranking).

Knee selection is a transparent, documented rule — not a black box:
  1. Admissible = busco_recovery >= --min-busco-recovery AND
     presence_recovery >= --min-presence-recovery AND fp_rate <= --max-fp
     (quality gates: a family definition that shreds single-copy orthologs, a coverage
     floor that hides real orthologs from BUSCO-confirmed outgroup copies, or a
     combination that floods false novelties, is disqualified regardless of recall).
  2. Rank admissible points by
     composite = recall + busco_recovery + presence_recovery - fp_rate (all ~[0,1]).
  3. The knee: among near-top composites, looser identity thresholds keep inflating
     n_novelties without buying recall/recovery, so ties break to FEWER novelties (then
     the shipped default params). Pick that point's parameters as the recommended default.
If no point is admissible, fall back to the best composite overall and flag it.

NOTE: presence_recovery, recall, and fp_rate are each left blank when their optional
input is missing (BUSCO_OUTGROUP_TABLES, or a controls CSV) rather than genuinely
measured -- this script tracks which metrics were actually measured anywhere in the
sweep and skips that metric's admissibility gate and composite contribution entirely
when it wasn't, printing a NOTE to say so. Earlier versions defaulted a blank cell to
its "worst" value (fp_rate=1.0, recall=0.0) and gated on it like a real measurement,
which made "no BUSCO_OUTGROUP_TABLES / no controls csv" look identical to "every grid
point has a 100% false-novelty rate" -- always check the raw metrics TSV and this
script's stderr NOTEs before trusting an admissibility verdict.
"""
import argparse
import csv
import sys

PARAM_COLS = ('min_seq_id', 'cov', 'hmm_evalue', 'hmm_cov', 'hmm_residues')
# metric -> (higher_is_better, worst_value_when_missing)
METRICS = {
    'n_families': (None, 0),
    'n_novelties': (None, 0),
    'busco_recovery': (True, 0.0),
    'presence_recovery': (True, 0.0),
    'recall': (True, 0.0),
    'fp_rate': (False, 1.0),
    'tblastn_removed': (None, 0),
    # run_ok (1/0, default 1 for older metrics files without this column -- see
    # bin/run_param_sweep.sh) marks whether the grid point's pipeline run actually
    # completed; a failed run's other metrics are partial/unreliable and must never be
    # admissible regardless of what they happen to read.
    'run_ok': (None, 1),
}
# Shipped default (nextflow.config) — the tie-break anchor when scores are equal.
DEFAULT_PARAMS = {'min_seq_id': 0.3, 'cov': 0.8, 'hmm_evalue': 1e-3,
                  'hmm_cov': 0.5, 'hmm_residues': 100}


def _num(value, fallback):
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def read_metrics(path):
    """Return (rows, available_metrics). available_metrics is the set of METRICS names
    that had at least one non-empty cell anywhere in the file -- an all-missing column
    (e.g. presence_recovery with no BUSCO_OUTGROUP_TABLES configured) is distinguished
    from a genuinely-measured-and-worst column, so its admissibility gate can be skipped
    rather than silently failing every row (see module docstring's NOTE)."""
    available = set()
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = []
        for raw in reader:
            point = {}
            for c in PARAM_COLS:
                point[c] = _num(raw.get(c), None)
            for m, (_hib, worst) in METRICS.items():
                cell = raw.get(m)
                if cell not in (None, ''):
                    available.add(m)
                point[m] = _num(cell, worst)
            rows.append(point)
    return rows, available


def score_points(rows, min_busco_recovery, min_presence_recovery, max_fp, available_metrics):
    # Each of these is optional infrastructure (curated controls, BUSCO outgroup tables)
    # that may simply not exist for a given clade/run -- an unmeasured metric must never
    # silently default to its "worst" value and drive a false rejection (or, for recall,
    # a false zero-credit) the way a genuinely-measured-and-bad value would. Only
    # busco_recovery is unconditional: it needs nothing but the run's own cluster_tsv/
    # families.tsv, so a missing value there really does mean the run failed, not that
    # an optional input was omitted.
    check_presence_recovery = 'presence_recovery' in available_metrics
    check_recall = 'recall' in available_metrics
    check_fp = 'fp_rate' in available_metrics
    for p in rows:
        admissible = p['run_ok'] >= 1 and p['busco_recovery'] >= min_busco_recovery
        if check_fp:
            admissible = admissible and p['fp_rate'] <= max_fp
        if check_presence_recovery:
            admissible = admissible and p['presence_recovery'] >= min_presence_recovery
        p['admissible'] = int(admissible)
        p['composite'] = round(
            p['busco_recovery']
            + (p['recall'] if check_recall else 0.0)
            + (p['presence_recovery'] if check_presence_recovery else 0.0)
            - (p['fp_rate'] if check_fp else 0.0), 6)
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
                    dest='min_busco_recovery',
                    help='clustering-quality admissibility gate (default 0.9)')
    ap.add_argument('--min-presence-recovery', type=float, default=0.9,
                    dest='min_presence_recovery',
                    help='presence-calling-quality admissibility gate (default 0.9); '
                         'skipped automatically if presence_recovery was never measured '
                         '(no BUSCO_OUTGROUP_TABLES in the sweep run)')
    ap.add_argument('--max-fp', type=float, default=0.05, dest='max_fp',
                    help='max control false-novelty rate to be admissible (default 0.05)')
    ap.add_argument('--output', required=True, help='scored sweep table TSV (sorted best-first)')
    args = ap.parse_args()

    rows, available_metrics = read_metrics(args.metrics)
    if not rows:
        sys.exit(f'no grid points in {args.metrics}')
    for metric, why in (
        ('presence_recovery', 'no BUSCO_OUTGROUP_TABLES?'),
        ('recall', 'no CONTROLS csv?'),
        ('fp_rate', 'no CONTROLS csv?'),
    ):
        if metric not in available_metrics:
            print(f'NOTE: {metric} was not measured in any grid point ({why}) -- '
                  'skipping its admissibility gate and composite contribution.',
                  file=sys.stderr)
    if not {'presence_recovery', 'recall', 'fp_rate'} & available_metrics:
        print('WARNING: presence_recovery, recall, and fp_rate are ALL unmeasured -- '
              'ranking rests entirely on busco_recovery (clustering quality, insensitive '
              'to hmm_cov/hmm_residues) plus the fewer-novelties tie-break. That is a weak '
              'signal for tuning hmm_cov/hmm_residues specifically -- see collate_sweep.py '
              'module docstring and consider running with BUSCO_OUTGROUP_TABLES and/or a '
              'controls CSV before trusting this recommendation.', file=sys.stderr)
    score_points(rows, args.min_busco_recovery, args.min_presence_recovery, args.max_fp,
                available_metrics)
    chosen, used_fallback = select_knee(rows)

    # Sort best-first: admissible first, then composite, then fewer novelties.
    ordered = sorted(rows, key=lambda p: (-p['admissible'], -p['composite'], p['n_novelties']))
    for p in ordered:
        p['chosen'] = int(p is chosen)

    cols = list(PARAM_COLS) + list(METRICS) + ['composite', 'admissible', 'chosen']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for p in ordered:
            w.writerow({c: p.get(c, '') for c in cols})

    if used_fallback:
        gates = ['run_ok', f'busco_recovery >= {args.min_busco_recovery}']
        if 'presence_recovery' in available_metrics:
            gates.append(f'presence_recovery >= {args.min_presence_recovery}')
        if 'fp_rate' in available_metrics:
            gates.append(f'fp_rate <= {args.max_fp}')
        print(f"WARNING: no grid point met the admissibility gates ({', '.join(gates)}); "
              'reporting the best composite overall.', file=sys.stderr)
    def _fmt(metric):
        if metric not in available_metrics:
            return 'n/a'
        return f"{chosen[metric]:.3f}"

    print(f"Recommended default: --family_min_seq_id {chosen['min_seq_id']} "
          f"--family_cov {chosen['cov']} --hmm_presence_evalue {chosen['hmm_evalue']:g} "
          f"--hmm_presence_cov {chosen['hmm_cov']:g} "
          f"--hmm_presence_min_residues {int(chosen['hmm_residues'])} "
          f"(composite {chosen['composite']:.3f}, recall {_fmt('recall')}, "
          f"busco {chosen['busco_recovery']:.3f}, "
          f"presence {_fmt('presence_recovery')}, fp {_fmt('fp_rate')}, "
          f"novelties {int(chosen['n_novelties'])})", file=sys.stderr)


if __name__ == '__main__':
    main()
