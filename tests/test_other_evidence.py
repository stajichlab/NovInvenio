"""Issue #208: other-group evidence sidecar from bin/build_presence_matrix.py.

Evidence only: the presence matrix and candidate list must not change.
"""
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'build_presence_matrix.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,In one,,in1.pep.fa,,In1,X
IN,In two,,in2.pep.fa,,In2,X
OUT,Out one,,out1.pep.fa,,Out1,Y
OUT,Out two,,out2.pep.fa,,Out2,Y
"""
HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\tqcov\n'
PARALOG = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def run(tmp_path, hits, paralogs='', evidence=True, extra=()):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'hits.tsv').write_text(HEADER + hits)
    (tmp_path / 'par.tsv').write_text(PARALOG + paralogs)
    cmd = [sys.executable, str(SCRIPT), '--hits', str(tmp_path / 'hits.tsv'),
           '--paralog-cutoffs', str(tmp_path / 'par.tsv'),
           '--config', str(tmp_path / 'config.csv'),
           '--ingroup-min-frac', '1.0', '--other-max-frac', '0.0',
           '--output-matrix', str(tmp_path / 'm.tsv'),
           '--output-candidates', str(tmp_path / 'c.txt'), *extra]
    if evidence:
        cmd += ['--output-other-evidence', str(tmp_path / 'ev.tsv.gz')]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not evidence:
        return None
    with gzip.open(tmp_path / 'ev.tsv.gz', 'rt') as fh:
        return pd.read_csv(fh, sep='\t', keep_default_na=False, dtype=str)


# g1 hits In2 (so it is an ingroup-wide protein) and Out1 weakly (1e-10, qcov 20).
# Its paralog p1 hits the same Out1 target far better (1e-100): filter 2 removes g1's hit.
PODAN_LIKE = (
    'g1\tg1_in2\t1e-50\t200\tIn1\tIn2\t90\n'
    'g1\tt1\t1e-10\t60\tIn1\tOut1\t20\n'
    'p1\tt1\t1e-100\t300\tIn1\tOut1\t95\n'
)
PODAN_PAR = 'g1\tp1\t100\t1e-30\n'


def test_paralog_filtered_cell_records_paralog_evalue_and_delta(tmp_path):
    ev = run(tmp_path, PODAN_LIKE, PODAN_PAR)
    row = ev[(ev.protein_id == 'g1') & (ev.other_proteome == 'Out1')].iloc[0]
    assert row.status == 'paralog_filtered'
    assert row.evidence == 'protein'
    assert row.source_proteome == 'In1'
    assert row.best_target_id == 't1'
    assert float(row.best_evalue) == 1e-10
    assert float(row.best_qcov) == 20
    assert float(row.max_qcov) == 20
    assert row.n_hits == '1'
    assert row.paralog_id == 'p1'
    assert float(row.paralog_evalue) == 1e-100
    assert abs(float(row.paralog_delta) - 90) < 1e-6


def test_only_other_group_cells_are_recorded(tmp_path):
    ev = run(tmp_path, PODAN_LIKE, PODAN_PAR)
    assert set(ev.other_proteome) <= {'Out1', 'Out2'}


def test_cell_without_any_hit_has_no_row(tmp_path):
    ev = run(tmp_path, PODAN_LIKE, PODAN_PAR)
    assert not ((ev.protein_id == 'g1') & (ev.other_proteome == 'Out2')).any()


def test_kept_hit_is_recorded_as_kept_without_paralog_fields(tmp_path):
    hits = ('g2\tg2_in2\t1e-50\t200\tIn1\tIn2\t90\n'
            'g2\tt9\t1e-40\t150\tIn1\tOut2\t80\n')
    ev = run(tmp_path, hits)
    row = ev[(ev.protein_id == 'g2') & (ev.other_proteome == 'Out2')].iloc[0]
    assert row.status == 'kept'
    assert row.paralog_id == '' and row.paralog_evalue == '' and row.paralog_delta == ''


def test_cell_is_kept_when_any_hit_survives(tmp_path):
    hits = ('g3\tg3_in2\t1e-50\t200\tIn1\tIn2\t90\n'
            'g3\tt1\t1e-10\t60\tIn1\tOut1\t20\n'
            'g3\tt2\t1e-12\t70\tIn1\tOut1\t55\n'
            'p3\tt1\t1e-100\t300\tIn1\tOut1\t95\n')
    ev = run(tmp_path, hits, 'g3\tp3\t100\t1e-30\n')
    row = ev[(ev.protein_id == 'g3') & (ev.other_proteome == 'Out1')].iloc[0]
    assert row.status == 'kept'
    assert row.n_hits == '2'
    assert float(row.max_qcov) == 55
    assert row.best_target_id == 't2'   # lowest evalue is 1e-12 -> t2


def test_coverage_floor_status(tmp_path):
    hits = ('g4\tg4_in2\t1e-50\t200\tIn1\tIn2\t90\n'
            'g4\tt1\t1e-30\t100\tIn1\tOut1\t20\n')
    ev = run(tmp_path, hits, extra=['--other-coverage-floor-qcov', '50'])
    row = ev[(ev.protein_id == 'g4') & (ev.other_proteome == 'Out1')].iloc[0]
    assert row.status == 'coverage_floor'


def test_evidence_does_not_change_matrix_or_candidates(tmp_path):
    t1, t2 = tmp_path / 'a', tmp_path / 'b'
    t1.mkdir(); t2.mkdir()
    run(t1, PODAN_LIKE, PODAN_PAR, evidence=True)
    run(t2, PODAN_LIKE, PODAN_PAR, evidence=False)
    assert (t1 / 'm.tsv').read_bytes() == (t2 / 'm.tsv').read_bytes()
    assert (t1 / 'c.txt').read_bytes() == (t2 / 'c.txt').read_bytes()


def test_header_only_when_there_are_no_other_group_hits(tmp_path):
    hits = 'g5\tg5_in2\t1e-50\t200\tIn1\tIn2\t90\n'
    ev = run(tmp_path, hits)
    assert len(ev) == 0
    assert 'paralog_delta' in ev.columns


def test_payload_with_a_zero_paralog_evalue_is_valid_json():
    """paralog_delta is +-inf when an E-value is 0; JSON has no Infinity, so the payload carries +-400
    (2026-10-08: agaricomycetes_pairwise novelties.html was blank, JSON.parse failed on 'Infinity')."""
    import json
    from other_evidence import DELTA_CAP, other_evidence_payload
    rows = [{'other_proteome': 'Cneo', 'best_target_id': 'X', 'best_evalue': '7.9e-10', 'best_bitscore': '56.6',
             'best_qcov': '50.5', 'status': 'paralog_filtered', 'paralog_id': 'Y', 'paralog_evalue': '0.0',
             'paralog_delta': 'inf', 'max_qcov': '50.5'},
            {'other_proteome': 'Rtor', 'best_target_id': 'Z', 'best_evalue': '0.0', 'best_bitscore': '900',
             'best_qcov': '99', 'status': 'kept', 'paralog_id': 'Y', 'paralog_evalue': '1e-20',
             'paralog_delta': '-inf', 'max_qcov': '99'}]
    out = other_evidence_payload(rows, [])
    assert [r[8] for r in out['p']] == [DELTA_CAP, -DELTA_CAP]
    json.dumps(out, allow_nan=False)
