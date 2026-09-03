"""Unit tests for bin/context_presence.py (issue #48)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'bin' / 'context_presence.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,In one,,in1.pep.fa,,In1,X
NEAR_INGROUP,Near one,,near1.pep.fa,,Near1,X
BROAD_OUTGROUP,Broad one,,broad1.pep.fa,,Broad1,Y
"""

HIT_HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'
PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


@pytest.fixture
def run_dir(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    return tmp_path


def run(run_dir, hits_text, candidates_text, paralog_text=None, scope=None):
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + hits_text)
    candidates_path = run_dir / 'candidates.txt'
    candidates_path.write_text(candidates_text)
    matrix_out = run_dir / 'context_presence.tsv'
    evalues_out = run_dir / 'context_presence.evalues.tsv'
    cmd = [sys.executable, str(SCRIPT),
           '--hits', str(hits_path),
           '--candidates', str(candidates_path),
           '--config', str(run_dir / 'config.csv'),
           '--output-matrix', str(matrix_out),
           '--output-evalues', str(evalues_out)]
    if paralog_text is not None:
        paralog_path = run_dir / 'paralog_cutoffs.tsv'
        paralog_path.write_text(PARALOG_HEADER + paralog_text)
        cmd += ['--paralog-cutoffs', str(paralog_path)]
    if scope is not None:
        cmd += ['--paralog-competition-scope', scope]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    matrix = pd.read_csv(matrix_out, sep='\t')
    evalues = pd.read_csv(evalues_out, sep='\t', dtype=str, keep_default_na=False)
    return matrix, evalues


def test_candidate_present_in_context_proteome(run_dir):
    matrix, evalues = run(run_dir, 'g1\tg1_near\t1e-10\t100\tIn1\tNear1\n', 'In1::g1\n')
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['Near1'] == 1
    assert row['Broad1'] == 0
    ev_row = evalues[evalues['protein_id'] == 'g1'].iloc[0]
    assert ev_row['Near1'] == '1e-10'
    assert ev_row['Broad1'] == ''


def test_only_candidates_are_emitted_not_their_paralogs(run_dir):
    # eif1 is searched (as hexA's paralog, for the competition check) but never emitted --
    # it isn't itself in candidates.txt.
    matrix, _ = run(
        run_dir,
        'hexA\thex1\t1e-10\t100\tIn1\tNear1\neif1\teif2\t1e-30\t150\tIn1\tNear1\n',
        'In1::hexA\n',
    )
    assert list(matrix['protein_id']) == ['hexA']


def test_flat_default_evalue_rejects_a_weak_hit(run_dir):
    # Filter 1 is now a flat significance floor (--default-evalue, 1e-5), not a
    # per-query paralog-derived one -- a weak hit (1e-3) is rejected regardless of
    # whether paralog data is supplied at all.
    matrix, _ = run(run_dir, 'g1\tg1_near\t1e-3\t20\tIn1\tNear1\n', 'In1::g1\n',
                    paralog_text='g1\tg1p\t42\t1e-8\n')
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['Near1'] == 0


def test_supplied_paralog_cutoff_no_longer_rejects_a_strong_hit(run_dir):
    # Regression test for the 2026-09-03 fix: a strong hit (1e-10, well within the flat
    # default 1e-5) used to be rejected outright whenever the candidate's own in-genome
    # paralog e-value (here 1e-20) was tighter than the hit -- an absolute-magnitude
    # proxy, not an actual test of whether the paralog explains this hit. The
    # paralog-competition filter (filter 2) is the real test for that, and it doesn't
    # fire here because the paralog (g1p) was never itself searched against Near1.
    matrix, _ = run(run_dir, 'g1\tg1_near\t1e-10\t100\tIn1\tNear1\n', 'In1::g1\n',
                    paralog_text='g1\tg1p\t42\t1e-20\n')
    row = matrix[matrix['protein_id'] == 'g1'].iloc[0]
    assert row['Near1'] == 1


def test_competition_target_scope_keeps_hexa_like_ortholog(run_dir):
    # Mirrors tests/test_build_presence_matrix.py's HEX-1/eIF5A case: hexA's own hit to
    # hex1 (1e-69) beats its paralog eif1's hit to the *same* target (5e-12), so 'target'
    # scope keeps the call even though eif1 hits harder somewhere else in Near1 (1e-70).
    hits = (
        'hexA\thex1\t1e-69\t230\tIn1\tNear1\n'
        'eif1\teif2\t1e-70\t233\tIn1\tNear1\n'
        'eif1\thex1\t5e-12\t45\tIn1\tNear1\n'
    )
    matrix, _ = run(run_dir, hits, 'In1::hexA\n',
                    paralog_text='hexA\teif1\t42\t4.2e-11\n', scope='target')
    row = matrix[matrix['protein_id'] == 'hexA'].iloc[0]
    assert row['Near1'] == 1


def test_competition_proteome_scope_drops_hexa_like_ortholog(run_dir):
    hits = (
        'hexA\thex1\t1e-69\t230\tIn1\tNear1\n'
        'eif1\teif2\t1e-70\t233\tIn1\tNear1\n'
        'eif1\thex1\t5e-12\t45\tIn1\tNear1\n'
    )
    matrix, _ = run(run_dir, hits, 'In1::hexA\n',
                    paralog_text='hexA\teif1\t42\t4.2e-11\n', scope='proteome')
    row = matrix[matrix['protein_id'] == 'hexA'].iloc[0]
    assert row['Near1'] == 0


def test_no_context_proteomes_in_config_yields_empty_columns(tmp_path):
    (tmp_path / 'config.csv').write_text('GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n'
                                         'IN,In one,,in1.pep.fa,,In1,X\n')
    matrix, evalues = run(tmp_path, 'g1\tg1_x\t1e-10\t100\tIn1\tSomewhere\n', 'In1::g1\n')
    assert list(matrix.columns) == ['protein_id', 'source_proteome']
    assert list(evalues.columns) == ['protein_id', 'source_proteome']
