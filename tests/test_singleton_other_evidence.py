"""Issue #208: evidence record for the novelty_discovery singleton path."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))

from other_evidence import evidence_from_records  # noqa: E402
from singleton_presence import score_singleton_hits  # noqa: E402

# s1 hits OutA weakly (1e-10, qcov 20) but its paralog p1 hits the same target at 1e-100.
HITS = [
    ('s1', 't1', 1e-10, 'OutA', 20.0),
    ('p1', 't1', 1e-100, 'OutA', 95.0),
    ('s1', 't2', 1e-30, 'OutB', 80.0),
    ('s1', 't3', 1e-12, 'Tgt2', 60.0),
    ('s1', 't4', 1e-20, 'OutC', 20.0),
]


def score(**kw):
    ev = []
    pres, _ = score_singleton_hits(HITS, {'s1'}, {'s1': 'p1'}, 1e-5, 'target',
                                   rescue_evalue=1e-20, evidence_out=ev, **kw)
    return pres, ev


def test_evidence_records_filtered_and_kept_hits():
    pres, ev = score()
    by = {(r[3]): r for r in ev}
    assert by['OutA'][5] == 'paralog_filtered' and by['OutA'][6] == 'p1' and by['OutA'][7] == 1e-100
    assert by['OutB'][5] == 'kept'
    assert 'OutA' not in pres and 'OutB' in pres


def test_evidence_does_not_change_presence():
    pres_with, _ = score()
    pres_without, _ = score_singleton_hits(HITS, {'s1'}, {'s1': 'p1'}, 1e-5, 'target', rescue_evalue=1e-20)
    assert dict(pres_with) == dict(pres_without)


def test_coverage_floor_status_is_recorded():
    _, ev = score(coverage_floor_qcov=50, floor_shorts={'OutB', 'OutC'})
    assert [r for r in ev if r[3] == 'OutB'][0][5] == 'kept'          # qcov 80 >= 50
    row = [r for r in ev if r[3] == 'OutC'][0]
    assert row[5] == 'coverage_floor'


def test_aggregation_keeps_only_other_group_cells():
    _, ev = score()
    df = evidence_from_records(ev, {'s1': 'Tgt1'}, {'OutA', 'OutB'})
    assert set(df.other_proteome) == {'OutA', 'OutB'}
    a = df[df.other_proteome == 'OutA'].iloc[0]
    assert a.source_proteome == 'Tgt1' and a.status == 'paralog_filtered'
    assert abs(a.paralog_delta - 90) < 1e-6 and a.max_qcov == 20.0


def test_empty_records_give_a_header_only_frame():
    df = evidence_from_records([], {}, {'OutA'})
    assert len(df) == 0 and 'paralog_delta' in df.columns
