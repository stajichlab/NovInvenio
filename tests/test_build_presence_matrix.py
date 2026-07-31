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


PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def run(run_dir, hits_text, query_group='IN', min_frac='1.0', other_max_frac='0.0',
        paralog_text=None, competition_scope=None):
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + hits_text)
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    cmd = [sys.executable, str(SCRIPT),
           '--hits', str(hits_path),
           '--config', str(run_dir / 'config.csv'),
           '--ingroup-min-frac', min_frac,
           '--query-group', query_group,
           '--other-max-frac', other_max_frac,
           '--output-matrix', str(matrix_out),
           '--output-candidates', str(candidates_out)]
    if paralog_text is not None:
        paralog_path = run_dir / 'paralog_cutoffs.tsv'
        paralog_path.write_text(PARALOG_HEADER + paralog_text)
        cmd += ['--paralog-cutoffs', str(paralog_path)]
    if competition_scope is not None:
        cmd += ['--paralog-competition-scope', competition_scope]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
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
    # h2 (Out1) also hits In1 — present in 1/2 of the ingroup, so not a loss
    # under the default strict absence (--other-max-frac 0.0).
    matrix, candidates = run(
        run_dir,
        'h2\th2_in1\t1e-10\t100\tOut1\tIn1\n',
        query_group='OUT',
    )
    assert candidates == []
    row = matrix[matrix['protein_id'] == 'h2'].iloc[0]
    assert row['Out1'] == 1 and row['In1'] == 1


def test_query_group_out_allows_nearly_missing_with_other_max_frac(run_dir):
    # Same h2 present in 1/2 of the ingroup (frac 0.5): excluded at the default
    # 0.0 above, but kept once --other-max-frac allows up to half the ingroup —
    # a "nearly missing" (lost from most, not all, of the ingroup) candidate.
    matrix, candidates = run(
        run_dir,
        'h2\th2_in1\t1e-10\t100\tOut1\tIn1\n',
        query_group='OUT',
        other_max_frac='0.5',
    )
    assert candidates == ['Out1::h2']


# --- Filter 2 (paralog competition) scope -----------------------------------
#
# Mirrors the pezizo5_fungi hexA/hex-1 case (see docs/hexA_filtering.md): hexA
# (In1) is a true ortholog of hex1 (In2's HEX-1 gene) and matches it strongly
# (1e-69), but hexA is a derivative of eIF5A. In2's genome carries both hex1 and
# its own eIF5A (eif2). hexA's within-genome paralog (eif1) hits In2's eif2 even
# harder (1e-70) than hexA hits hex1 — but only on a *different* target gene.
HEXA_HITS = (
    'hexA\thex1\t1e-69\t230\tIn1\tIn2\n'   # hexA -> hex1: the real ortholog hit
    'eif1\teif2\t1e-70\t233\tIn1\tIn2\n'   # eIF5A paralog -> In2's eIF5A: wins proteome-wide
    'eif1\thex1\t5e-12\t45\tIn1\tIn2\n'    # eIF5A paralog -> hex1: loses on this target
)
# hexA's within-genome paralog is eif1 (cutoff loose enough that filter 1 keeps
# the 1e-69 hit); eif1's paralog is hexA.
HEXA_PARALOGS = 'hexA\teif1\t42\t4.2e-11\neif1\thexA\t42\t4.2e-11\n'


def test_competition_proteome_scope_drops_hexa_like_ortholog(run_dir):
    # Default 'proteome' scope: the eIF5A paralog out-scores hexA *anywhere* in
    # In2, so hexA's only qualifying hit is dropped. With no surviving cross-hit
    # hexA gets no matrix row at all -> never a candidate.
    matrix, candidates = run(run_dir, HEXA_HITS, paralog_text=HEXA_PARALOGS,
                             competition_scope='proteome')
    assert (matrix['protein_id'] == 'hexA').sum() == 0
    assert 'In1::hexA' not in candidates


def test_competition_target_scope_keeps_hexa_like_ortholog(run_dir):
    # 'target' scope: on the shared target gene hex1, hexA (1e-69) beats its
    # paralog eif1 (5e-12), so the call survives -> present in 2/2 ingroup -> candidate.
    matrix, candidates = run(run_dir, HEXA_HITS, paralog_text=HEXA_PARALOGS,
                             competition_scope='target')
    row = matrix[matrix['protein_id'] == 'hexA'].iloc[0]
    assert row['In1'] == 1 and row['In2'] == 1
    assert 'In1::hexA' in candidates


def test_output_evalues_sidecar_matches_presence_calls(run_dir):
    # issue #44: --output-evalues emits the qualifying hit's e-value alongside each
    # presence=1 cell, empty for absence and for the protein's own source proteome.
    hits_path = run_dir / 'hits.tsv'
    hits_path.write_text(HIT_HEADER + 'g1\tg1_in2\t1e-10\t100\tIn1\tIn2\n')
    matrix_out = run_dir / 'matrix.tsv'
    candidates_out = run_dir / 'candidates.txt'
    evalues_out = run_dir / 'evalues.tsv'
    subprocess.run([
        sys.executable, str(SCRIPT),
        '--hits', str(hits_path),
        '--config', str(run_dir / 'config.csv'),
        '--ingroup-min-frac', '1.0',
        '--query-group', 'IN',
        '--other-max-frac', '0.0',
        '--output-matrix', str(matrix_out),
        '--output-candidates', str(candidates_out),
        '--output-evalues', str(evalues_out),
    ], check=True, capture_output=True, text=True)

    evalues = pd.read_csv(evalues_out, sep='\t', dtype=str, keep_default_na=False)
    row = evalues[evalues['protein_id'] == 'g1'].iloc[0]
    assert row['In2'] == '1e-10'
    assert row['In1'] == ''   # source proteome: presence is definitional, not a hit
    assert row['Out1'] == ''  # absent: no qualifying hit
