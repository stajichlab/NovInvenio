"""Unit tests for lib/family_presence.py (2026-09-03 fixes: drop circular per-family
calibration from presence gating; fix coverage pooled across unrelated targets; fix
presence being gated by only the single best-E target's own coverage instead of "any
target clears both gates")."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from family_presence import family_presence_by_proteome, parse_domtblout  # noqa: E402


def _domtblout_line(target, query, full_e=1e-10, qlen=100, hmm_from=1, hmm_to=60):
    f = ['-'] * 23
    f[0] = target
    f[3] = query
    f[5] = str(qlen)
    f[6] = f'{full_e:g}'
    f[15] = str(hmm_from)
    f[16] = str(hmm_to)
    return ' '.join(f)


def _write_domtblout(path, rows):
    """rows: iterable of (target, query, full_e, qlen, hmm_from, hmm_to) with defaults."""
    with open(path, 'w') as fh:
        fh.write("# hmmsearch domtblout\n")
        for row in rows:
            fh.write(_domtblout_line(*row) + "\n")


def test_coverage_is_not_pooled_across_different_targets(tmp_path):
    # Two different, unrelated target proteins each hit by a weak, partial (30%-coverage)
    # domain. Pooling would sum to 60% coverage and clear a 0.5 threshold; per-target
    # coverage must not -- neither target qualifies on its own, so famX is absent.
    dom = tmp_path / 'x.domtblout'
    _write_domtblout(dom, [
        ('targetA', 'famX', 1e-3, 100, 1, 30),   # 30% coverage on targetA
        ('targetB', 'famX', 1e-3, 100, 1, 30),   # 30% coverage on targetB (different target)
    ])
    result = parse_domtblout(dom, default_evalue=1e-2, min_coverage=0.5)
    assert 'famX' not in result


def test_coverage_merges_overlapping_domains_on_the_same_target(tmp_path):
    dom = tmp_path / 'x.domtblout'
    _write_domtblout(dom, [
        ('targetA', 'famX', 1e-3, 100, 1, 40),
        ('targetA', 'famX', 1e-3, 100, 30, 70),  # overlaps the first domain
    ])
    # Merged span 1-70 = 70% coverage, not 40+41=81% double-counted -- qualifies at 0.70
    # but not at 0.80.
    assert 'famX' in parse_domtblout(dom, default_evalue=1e-2, min_coverage=0.70)
    assert 'famX' not in parse_domtblout(dom, default_evalue=1e-2, min_coverage=0.80)


def test_any_qualifying_target_counts_even_if_not_the_best_evalue_target(tmp_path):
    # Regression test for the 2026-09-03 fix: targetA has the best E-value (1e-50) but
    # only fragmentary coverage (20%); targetB has a weaker (but still significant)
    # E-value (1e-10) with full coverage (90%). Gating on the single best-E target's own
    # coverage would wrongly call famX absent (targetA's 20% fails 0.5) even though
    # targetB is a genuine, well-covered hit. Presence must be true (any target
    # qualifies), and the reported E-value must be targetB's own (1e-10), not mixed with
    # targetA's.
    dom = tmp_path / 'x.domtblout'
    _write_domtblout(dom, [
        ('targetA', 'famX', 1e-50, 100, 1, 20),   # best E-value, only 20% coverage
        ('targetB', 'famX', 1e-10, 100, 1, 90),   # weaker E-value, 90% coverage
    ])
    result = parse_domtblout(dom, default_evalue=1e-5, min_coverage=0.5)
    assert result['famX'] == 1e-10


def test_family_presence_by_proteome_uses_flat_default_not_calibration(tmp_path):
    # family_presence_by_proteome no longer accepts a calibrated-thresholds argument --
    # every proteome uses the same flat default_evalue.
    dom = tmp_path / 'D1.domtblout'
    _write_domtblout(dom, [('d1_x', 'famX', 3e-160, 500, 1, 500)])

    presence = family_presence_by_proteome([dom], default_evalue=1e-3, min_coverage=0.5)
    assert 'famX' in presence['D1']


def test_family_presence_by_proteome_respects_coverage_floor(tmp_path):
    dom = tmp_path / 'D1.domtblout'
    _write_domtblout(dom, [('d1_x', 'famX', 1e-10, 500, 1, 100)])  # 20% coverage

    presence = family_presence_by_proteome([dom], default_evalue=1e-3, min_coverage=0.5)
    assert 'famX' not in presence.get('D1', set())


def test_family_presence_by_proteome_any_target_qualifies(tmp_path):
    # Same any-target regression as above, through the family_presence_by_proteome
    # entry point bin/novelty_screen.py actually uses.
    dom = tmp_path / 'D1.domtblout'
    _write_domtblout(dom, [
        ('d1_fragment', 'famX', 1e-50, 100, 1, 20),
        ('d1_full', 'famX', 1e-10, 100, 1, 90),
    ])
    presence = family_presence_by_proteome([dom], default_evalue=1e-5, min_coverage=0.5)
    assert 'famX' in presence['D1']
