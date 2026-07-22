"""Unit tests for bin/score_controls.py (ADR-0002 Q8 validation battery, issue #5).

Covers the pure scoring/resolution helpers directly and an end-to-end subprocess run over
a fixture controls CSV + presence matrix: positive-recall, negative-FP, the protein_id /
fasta / busco anchor types, and placeholder skipping.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'score_controls.py'

# Import the bin script as a module to exercise its pure functions directly.
_spec = importlib.util.spec_from_file_location('score_controls', BIN)
sc = importlib.util.module_from_spec(_spec)
sys.path.insert(0, str(REPO / 'lib'))
_spec.loader.exec_module(sc)


CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Species one,,In1.faa,In1.fna,In1,Pezizomycotina
IN,Species two,,In2.faa,In2.fna,In2,Pezizomycotina
OUT,Species three,,Out1.faa,Out1.fna,Out1,Saccharomycotina
"""

# Family A (pA1/pA2): present in both ingroup, absent from the outgroup → NOVEL.
# Family B (pB1/pB2): present everywhere → NOT novel (a negative control's "core").
MATRIX = """\
protein_id\tsource_proteome\tIn1\tIn2\tOut1
pA1\tIn1\t1\t1\t0
pA2\tIn2\t1\t1\t0
pB1\tIn1\t1\t1\t1
pB2\tIn2\t1\t1\t1
"""

CLUSTER = "pA1\tpA1\npA1\tpA2\npB1\tpB1\npB1\tpB2\npC1\tpC1\n"

FAMILIES = (
    "family_index\trepresentative_id\tn_members\n"
    "fam_000001\tpA1\t2\n"
    "fam_000002\tpB1\t2\n"
)


# --------------------------------------------------------------------- pure functions
def test_family_call_novel_vs_not():
    ingroup, outgroup = ['In1', 'In2'], ['Out1']
    novel = sc.family_call({'In1': 1, 'In2': 1, 'Out1': 0}, ingroup, outgroup, 0.75, 0.0)
    assert novel == 'novel'
    # present in the outgroup → disqualified
    core = sc.family_call({'In1': 1, 'In2': 1, 'Out1': 1}, ingroup, outgroup, 0.75, 0.0)
    assert core == 'not_novel'
    # below the ingroup fraction → not a candidate
    sparse = sc.family_call({'In1': 1, 'In2': 0, 'Out1': 0}, ingroup, outgroup, 0.75, 0.0)
    assert sparse == 'not_novel'


def test_family_presence_vector_ors_members(tmp_path):
    matrix = pd.read_csv(pd.io.common.StringIO(MATRIX), sep='\t')
    cols = ['In1', 'In2', 'Out1']
    vec = sc.family_presence_vector(matrix, ['pA1', 'pA2'], cols)
    assert vec == {'In1': 1, 'In2': 1, 'Out1': 0}
    assert sc.family_presence_vector(matrix, ['nope'], cols) is None


def test_is_placeholder():
    assert sc.is_placeholder({'control_id': 'EXAMPLE_POS01', 'anchor': 'pA1'})
    assert sc.is_placeholder({'control_id': 'POS01', 'anchor': '<protein_id_here>'})
    assert not sc.is_placeholder({'control_id': 'POS01', 'anchor': 'pA1'})


def test_best_family_from_domtblout(tmp_path):
    dom = tmp_path / 'anchor.domtblout'
    # two hits; the lower full-seq E-value (col idx 6) wins.
    row = lambda q, e: ' '.join(['t', '-', '-', q, '-', '100', e] + ['-'] * 16)  # noqa: E731
    dom.write_text("# header\n" + row('pB1', '1e-5') + "\n" + row('pA1', '1e-30') + "\n")
    assert sc.best_family_from_domtblout(str(dom)) == 'pA1'


def test_load_busco_map(tmp_path):
    m = tmp_path / 'busco.tsv'
    m.write_text("BUSCO123\tpB1\nBUSCO456\tpA2\n")
    assert sc.load_busco_map(str(m)) == {'BUSCO123': 'pB1', 'BUSCO456': 'pA2'}
    assert sc.load_busco_map(None) == {}


# ------------------------------------------------------------------------ end-to-end
def _setup(tmp_path, controls_csv):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    (tmp_path / 'families.tsv').write_text(FAMILIES)
    (tmp_path / 'controls.csv').write_text(controls_csv)


def _run(tmp_path, *extra):
    out = tmp_path / 'results.tsv'
    cmd = [
        sys.executable, str(BIN),
        '--controls', str(tmp_path / 'controls.csv'),
        '--matrix', str(tmp_path / 'matrix.tsv'),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--families', str(tmp_path / 'families.tsv'),
        '--config', str(tmp_path / 'config.csv'),
        '--output', str(out),
        *extra,
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    results = pd.read_csv(out, sep='\t').set_index('control_id')
    summary = pd.read_csv(tmp_path / 'results.summary.tsv', sep='\t').set_index('metric')['value']
    return results, summary


CONTROLS = """\
control_id,class,expected_call,anchor_type,anchor,proteome_short,gene_name,expected_origin,source,notes
POS01,positive,novel,protein_id,pA1,In1,,clade,prior_run,true novelty
NEG01,negative,core,protein_id,pB1,In1,rpb2,pan,housekeeping,conserved
NEG_FP,negative,core,protein_id,pA2,In2,,pan,synthetic,should trip an FP
EXAMPLE_POS99,positive,novel,protein_id,<placeholder>,In1,,clade,template,skip me
"""


def test_end_to_end_protein_id(tmp_path):
    _setup(tmp_path, CONTROLS)
    results, summary = _run(tmp_path)

    assert results.loc['POS01', 'outcome'] == 'hit'
    assert results.loc['POS01', 'resolved_family'] == 'pA1'
    assert results.loc['NEG01', 'outcome'] == 'tn'
    assert results.loc['NEG_FP', 'outcome'] == 'fp'  # family A is novel; expected core
    assert 'EXAMPLE_POS99' not in results.index  # placeholder skipped

    assert int(summary['n_placeholder_skipped']) == 1
    assert float(summary['recall']) == 1.0          # 1/1 positives recovered
    assert float(summary['fp_rate']) == 0.5         # 1 of 2 negatives tripped


BUSCO_CONTROLS = """\
control_id,class,expected_call,anchor_type,anchor,proteome_short,gene_name,expected_origin,source,notes
BUSCO01,negative,core,busco,12345at4751,,,pan-fungal,fungi_odb12,single-copy ortholog
"""


def test_busco_anchor_unresolved_without_map(tmp_path):
    _setup(tmp_path, BUSCO_CONTROLS)
    results, summary = _run(tmp_path)
    assert results.loc['BUSCO01', 'actual_call'] == 'unresolved'
    assert int(summary['n_unresolved']) == 1


def test_busco_anchor_resolves_with_map(tmp_path):
    _setup(tmp_path, BUSCO_CONTROLS)
    (tmp_path / 'busco.tsv').write_text("12345at4751\tpB1\n")
    results, _ = _run(tmp_path, '--busco-map', str(tmp_path / 'busco.tsv'))
    # BUSCO → pB1 → family B (present everywhere) → not novel → correct negative (tn).
    assert results.loc['BUSCO01', 'resolved_family'] == 'pB1'
    assert results.loc['BUSCO01', 'outcome'] == 'tn'


FASTA_CONTROLS = """\
control_id,class,expected_call,anchor_type,anchor,proteome_short,gene_name,expected_origin,source,notes
FASTA01,positive,novel,fasta,seqs/missing.faa,,,clade,prior_run,anchor by sequence
"""


def test_fasta_anchor_missing_file_is_unresolved(tmp_path):
    _setup(tmp_path, FASTA_CONTROLS)
    results, _ = _run(tmp_path)
    assert results.loc['FASTA01', 'actual_call'] == 'unresolved'
    assert 'not found' in results.loc['FASTA01', 'note']
