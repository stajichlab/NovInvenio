"""Issue #208: other-group evidence in the report payload ('oe' field)."""
import gzip
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))
sys.path.insert(0, str(REPO / 'tests'))

from test_report_data import payload_for, rows_by_id, run_dir, samples  # noqa: E402,F401

EV_HEAD = ('protein_id\tsource_proteome\tother_proteome\tevidence\tbest_target_id\tbest_evalue\t'
           'best_bitscore\tbest_qcov\tmax_qcov\tn_hits\tstatus\tparalog_id\tparalog_evalue\tparalog_delta\n')
COV_HEAD = 'protein_id\trep_id\tgenome\tn_hsps\tbest_evalue\tquery_span_cov\tspan_start\tspan_end\n'


def _gz(path, text):
    with gzip.open(path, 'wt') as fh:
        fh.write(text)
    return path


@pytest.fixture
def evidence(run_dir):
    ev = _gz(run_dir / 'ev.tsv.gz', EV_HEAD +
             'n1\tNcra\tSpom\tprotein\tt1\t1e-10\t60\t20\t20\t1\tparalog_filtered\tp1\t1e-100\t90\n')
    cov = _gz(run_dir / 'cov.tsv.gz', COV_HEAD + 'n2\tn2\tScer\t2\t1e-8\t0.800000\t5\t95\n')
    return ev, cov


def test_oe_carries_c_and_per_cell_evidence(run_dir, samples, evidence):
    ev, cov = evidence
    p = payload_for(run_dir, samples, other_evidence_path=ev, tblastn_coverage_path=cov)
    oe = p['fields'].index('oe')
    rows = rows_by_id(p)
    n1 = rows['n1'][oe]
    assert n1['c'] == 20
    assert n1['p'][0][:2] == ['Spom', 't1']
    assert n1['p'][0][5] == 'paralog_filtered'
    assert n1['p'][0][6] == 'p1'
    n2 = rows['n2'][oe]
    assert n2['c'] == 80.0
    assert n2['t'][0][0] == 'Scer' and n2['t'][0][1] == 0.8


def test_novelty_without_hits_gets_empty_evidence(run_dir, samples, evidence):
    ev, cov = evidence
    p = payload_for(run_dir, samples, other_evidence_path=ev, tblastn_coverage_path=cov)
    oe = p['fields'].index('oe')
    # 'shared' is not a novelty candidate -> not measured; so use a novelty with no rows.
    ev2 = _gz(run_dir / 'ev2.tsv.gz', EV_HEAD)
    p2 = payload_for(run_dir, samples, other_evidence_path=ev2)
    assert rows_by_id(p2)['n1'][oe] == {'c': None, 'p': [], 't': []}
    assert rows_by_id(p)['shared'][oe] is None


def test_flags_and_threshold(run_dir, samples, evidence):
    ev, cov = evidence
    p = payload_for(run_dir, samples, other_evidence_path=ev, tblastn_coverage_path=cov,
                    other_signal_qcov=65)
    assert p['has_other_evidence'] is True
    assert p['other_signal_qcov'] == 65


def test_no_sidecars_means_not_measured(run_dir, samples):
    p = payload_for(run_dir, samples)
    assert p['has_other_evidence'] is False
    oe = p['fields'].index('oe')
    assert all(r[oe] is None for r in p['rows'])
