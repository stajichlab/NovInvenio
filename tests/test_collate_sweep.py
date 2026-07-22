"""Unit tests for bin/collate_sweep.py (ADR-0002 Q8 sweep knee-selection, issue #6).

Covers the admissibility gates (BUSCO recovery + FP), the composite ranking, the
fewer-novelties knee tie-break, and the no-admissible-point fallback.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'collate_sweep.py'

_spec = importlib.util.spec_from_file_location('collate_sweep', BIN)
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)

HEADER = ('min_seq_id\tcov\thmm_evalue\tn_families\tn_novelties\t'
          'busco_recovery\trecall\tfp_rate\ttblastn_removed\n')


def _metrics(tmp_path, rows):
    p = tmp_path / 'metrics.tsv'
    p.write_text(HEADER + ''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    return p


def _score(tmp_path, rows, **kw):
    rows = cs.read_metrics(_metrics(tmp_path, rows))
    cs.score_points(rows, kw.get('min_busco_recovery', 0.9), kw.get('max_fp', 0.05))
    return rows


def test_highest_composite_admissible_point_wins(tmp_path):
    rows = _score(tmp_path, [
        # id, cov, e, n_fam, n_nov, busco, recall, fp, tb
        [0.3, 0.5, '1e-3', 500, 100, 0.95, 1.0, 0.0, 3],
        [0.2, 0.5, '1e-3', 700, 200, 0.96, 1.0, 0.0, 5],  # higher composite (+busco)
    ])
    chosen, fallback = cs.select_knee(rows)
    assert not fallback
    assert (chosen['min_seq_id'], chosen['cov']) == (0.2, 0.5)


def test_fewer_novelties_breaks_composite_ties(tmp_path):
    rows = _score(tmp_path, [
        [0.3, 0.5, '1e-3', 500, 100, 0.95, 1.0, 0.0, 3],   # same composite, fewer novelties
        [0.2, 0.5, '1e-3', 900, 400, 0.95, 1.0, 0.0, 9],   # same composite, more novelties
    ])
    chosen, _ = cs.select_knee(rows)
    assert chosen['n_novelties'] == 100
    assert chosen['min_seq_id'] == 0.3


def test_admissibility_gates_exclude_bad_quality(tmp_path):
    rows = _score(tmp_path, [
        [0.2, 0.5, '1e-3', 900, 500, 0.80, 1.0, 0.0, 1],   # busco below 0.9 gate
        [0.5, 0.7, '1e-5', 300, 40, 0.92, 0.8, 0.2, 0],    # fp above 0.05 gate
        [0.3, 0.5, '1e-3', 500, 100, 0.93, 0.9, 0.0, 3],   # the only admissible point
    ])
    assert [p['admissible'] for p in rows] == [0, 0, 1]
    chosen, fallback = cs.select_knee(rows)
    assert not fallback
    assert chosen['min_seq_id'] == 0.3


def test_fallback_when_no_point_is_admissible(tmp_path):
    rows = _score(tmp_path, [
        [0.2, 0.5, '1e-3', 900, 500, 0.60, 1.0, 0.3, 1],
        [0.5, 0.7, '1e-5', 300, 40, 0.70, 0.9, 0.2, 0],   # best composite of the two
    ])
    chosen, fallback = cs.select_knee(rows)
    assert fallback
    assert chosen['min_seq_id'] == 0.5


def test_missing_cells_treated_as_worst(tmp_path):
    p = tmp_path / 'm.tsv'
    # second point has blank busco_recovery/recall/fp → worst values, not admissible.
    p.write_text(HEADER +
                 "0.3\t0.5\t1e-3\t500\t100\t0.95\t1.0\t0.0\t3\n"
                 "0.4\t0.7\t1e-5\t300\t40\t\t\t\t0\n")
    rows = cs.read_metrics(p)
    cs.score_points(rows, 0.9, 0.05)
    assert rows[1]['busco_recovery'] == 0.0
    assert rows[1]['fp_rate'] == 1.0
    assert rows[1]['admissible'] == 0


def test_end_to_end_writes_scored_table(tmp_path):
    metrics = _metrics(tmp_path, [
        [0.3, 0.5, '1e-3', 500, 100, 0.95, 1.0, 0.0, 3],
        [0.2, 0.5, '1e-3', 900, 400, 0.95, 1.0, 0.0, 9],
        [0.5, 0.7, '1e-5', 300, 40, 0.80, 0.9, 0.0, 0],
    ])
    out = tmp_path / 'scores.tsv'
    r = subprocess.run(
        [sys.executable, str(BIN), '--metrics', str(metrics), '--output', str(out)],
        check=True, capture_output=True, text=True)
    lines = out.read_text().splitlines()
    header = lines[0].split('\t')
    assert 'composite' in header and 'chosen' in header
    chosen_col = header.index('chosen')
    chosen_rows = [ln for ln in lines[1:] if ln.split('\t')[chosen_col] == '1']
    assert len(chosen_rows) == 1
    # Recommendation is printed to stderr for the operator to wire into nextflow.config.
    assert 'Recommended default' in r.stderr
