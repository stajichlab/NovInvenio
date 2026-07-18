import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'build_presence_matrix.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,In one,,in1.pep.fa,,In1,X
IN,In two,,in2.pep.fa,,In2,X
OUT,Out one,,out1.pep.fa,,Out1,Y
"""

HIT_HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'


@pytest.fixture
def run_dir(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    return tmp_path


def run(run_dir, hits_text, query_group='IN', min_frac='1.0'):
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + hits_text)
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    subprocess.run(
        [sys.executable, str(SCRIPT),
         '--hits', str(hits_path),
         '--config', str(run_dir / 'config.csv'),
         '--ingroup-min-frac', min_frac,
         '--query-group', query_group,
         '--output-matrix', str(matrix_out),
         '--output-candidates', str(candidates_out)],
        check=True, capture_output=True, text=True,
    )
    matrix = pd.read_csv(matrix_out, sep='\t')
    candidates = candidates_out.read_text().splitlines() if candidates_out.stat().st_size else []
    return matrix, candidates


def test_default_query_group_is_ingroup_novelty_direction(run_dir):
    # g1 (In1) hits In2 but never Out1 — a novelty candidate under the default
    # (unchanged) ingroup-as-query behaviour.
    matrix, candidates = run(run_dir, 'g1\tg1_in2\t1e-10\t100\tIn1\tIn2\n')
    assert candidates == ['In1::g1']
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['In1'] == 1 and row['In2'] == 1 and row['Out1'] == 0


def test_query_group_out_finds_a_loss_candidate(run_dir):
    # h1 (Out1) hits nothing in the ingroup — a candidate loss: conserved in
    # the (single-member) outgroup query set, absent from the ingroup.
    matrix, candidates = run(run_dir, 'h1\th1_self\t1e-10\t100\tOut1\tOut1\n',
                              query_group='OUT')
    assert candidates == ['Out1::h1']
    row = matrix[matrix['protein_id'] == 'h1'].iloc[0]
    assert row['Out1'] == 1 and row['In1'] == 0 and row['In2'] == 0


def test_query_group_out_excludes_hits_present_in_ingroup(run_dir):
    # h2 (Out1) also hits In1 — present in the ingroup, so not a loss.
    matrix, candidates = run(
        run_dir,
        'h2\th2_in1\t1e-10\t100\tOut1\tIn1\n',
        query_group='OUT',
    )
    assert candidates == []
    row = matrix[matrix['protein_id'] == 'h2'].iloc[0]
    assert row['Out1'] == 1 and row['In1'] == 1
