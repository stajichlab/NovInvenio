"""Gating of the optional Starship / SM-cluster overlap step (STARSHIP_OVERLAP, SM_OVERLAP).

Drives a real `nextflow run pangenome.nf -preview` (same fixture as
tests/test_bool_cli_params.py): no new process without the params, each part alone with
its own param, and the workflow-start table check stops a bad table before any task."""
import re
import shutil
import subprocess

import pytest

from tests.test_bool_cli_params import REPO, _dag

NEXTFLOW = shutil.which('nextflow')


@pytest.fixture
def minimal_study(tmp_path):
    """Same 3-strain study as tests/test_bool_cli_params.py's fixture (A, B = IN; C = OUT)."""
    data = tmp_path / 'data'
    for sub in ('pep', 'dna', 'gff3'):
        (data / sub).mkdir(parents=True)
    for s in 'ABC':
        (data / 'pep' / f'{s}.fa').touch()
        (data / 'dna' / f'{s}.fa').touch()
        (data / 'gff3' / f'{s}.gff3').touch()
    ss = tmp_path / 'ss.csv'
    ss.write_text('GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n'
                  'IN,Sp a,,A.fa,A.fa,A.gff3,A,\nIN,Sp a,,B.fa,B.fa,B.gff3,B,\n'
                  'OUT,Sp o,,C.fa,C.fa,C.gff3,C,\n')
    return tmp_path, ss, data
pytestmark = pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')

WF = (REPO / 'workflows' / 'pangenome_profile.nf').read_text()
ELS = 'element_id\tstrain\tcontig\tstart\tend\tlabel\tconfidence\nE1\tA\tc1\t1\t100\tShip\thigh\n'
SM = 'cluster_id\tstrain\tcontig\tstart\tend\ttype\nK1\tA\tc1\t1\t50\tT1PKS\n'


def _processes(dag: str) -> set[str]:
    return set(re.findall(r'\b([A-Z][A-Z0-9_]{2,})\b', dag)) & set(re.findall(r'process\s+([A-Z0-9_]+)', _all_nf()))


def _all_nf() -> str:
    return '\n'.join(p.read_text() for p in (REPO / 'modules').rglob('*.nf'))


def _run(minimal_study, *extra):
    d, ss, data = minimal_study
    return subprocess.run(
        [NEXTFLOW, '-q', 'run', str(REPO / 'pangenome.nf'), '-preview', '-with-dag', str(d / 'dag.mmd'),
         '--pangenome_samplesheet', str(ss), '--pangenome_data_dir', str(data),
         '--outdir', str(d / 'out'), *extra],
        cwd=d, capture_output=True, text=True, timeout=600)


def _table(minimal_study, name, text):
    p = minimal_study[0] / name
    p.write_text(text)
    return str(p)


def test_both_calls_are_gated_in_the_workflow():
    for call, param in (('STARSHIP_OVERLAP(', 'params.pangenome_starship_annotation'),
                        ('SM_OVERLAP(', 'params.pangenome_sm_clusters')):
        before = WF.split(call, 1)[0]
        assert before.rstrip().endswith(f'if ({param}) {{'), call


def test_default_graph_has_no_overlap_process(minimal_study):
    dag = _dag(minimal_study)
    assert 'STARSHIP_OVERLAP' not in dag and 'SM_OVERLAP' not in dag
    assert 'BUILD_ISLANDS' in dag


def test_starship_only(minimal_study):
    dag = _dag(minimal_study, '--pangenome_starship_annotation', _table(minimal_study, 'e.tsv', ELS))
    assert 'STARSHIP_OVERLAP' in dag and 'SM_OVERLAP' not in dag


def test_sm_only(minimal_study):
    dag = _dag(minimal_study, '--pangenome_sm_clusters', _table(minimal_study, 's.tsv', SM))
    assert 'SM_OVERLAP' in dag and 'STARSHIP_OVERLAP' not in dag


def test_both_and_no_other_new_process(minimal_study):
    base = _processes(_dag(minimal_study))
    both = _processes(_dag(minimal_study, '--pangenome_starship_annotation', _table(minimal_study, 'e.tsv', ELS),
                           '--pangenome_sm_clusters', _table(minimal_study, 's.tsv', SM)))
    assert both - base == {'STARSHIP_OVERLAP', 'SM_OVERLAP'}
    assert base - both == set()


@pytest.mark.parametrize('param, text, msg', [
    ('--pangenome_starship_annotation', ELS.replace('\tA\tc1', '\tZZ\tc1'), "strain 'ZZ' is not a Short ID"),
    ('--pangenome_starship_annotation', ELS.replace('\t1\t100\t', '\t200\t100\t'), 'start 200 > end 100'),
    ('--pangenome_starship_annotation', ELS.replace('\t1\t100\t', '\tx\t100\t'), "start 'x' is not an integer"),
    ('--pangenome_sm_clusters', SM.replace('\ttype', '\tkind'), "missing column(s) [type]"),
    ('--pangenome_sm_clusters', SM + 'K1\tB\tc1\t1\t5\tNRPS\n', "duplicate cluster_id 'K1'"),
])
def test_bad_table_stops_before_any_task(minimal_study, param, text, msg):
    proc = _run(minimal_study, param, _table(minimal_study, 'bad.tsv', text))
    assert proc.returncode != 0
    assert msg in proc.stdout + proc.stderr
