"""Issue #191: under Nextflow 26 a command-line `--flag false` reaches the script
as the String "false", which is truthy. Every boolean param is read through
Helpers.asBool(); this drives a real `nextflow run -preview` and checks the
rescue branch really disappears from the workflow graph."""
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
NEXTFLOW = shutil.which('nextflow')

pytestmark = pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')


@pytest.fixture
def minimal_study(tmp_path):
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


def _dag(minimal_study, *extra):
    d, ss, data = minimal_study
    dag = d / 'dag.mmd'
    proc = subprocess.run(
        [NEXTFLOW, '-q', 'run', str(REPO / 'pangenome.nf'), '-preview', '-with-dag', str(dag),
         '--pangenome_samplesheet', str(ss), '--pangenome_data_dir', str(data),
         '--outdir', str(d / 'out'), *extra],
        cwd=d, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return dag.read_text()


def test_rescue_is_on_by_default(minimal_study):
    assert 'TBLASTN_PER_STRAIN' in _dag(minimal_study)


@pytest.mark.parametrize('flag', [['--pangenome_rescue_enable', 'false'],
                                  ['--pangenome_rescue_enable=false']])
def test_cli_false_turns_rescue_off(minimal_study, flag):
    dag = _dag(minimal_study, *flag)
    assert 'TBLASTN_PER_STRAIN' not in dag
    assert 'EXTRACT_ABSENT_QUERIES' not in dag
