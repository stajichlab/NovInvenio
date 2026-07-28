"""End-to-end integration test for --cluster_tool novelty_discovery (issue #29).

Runs the real Nextflow pipeline (NOVELTY_DISCOVERY -> NOVELTY_SCREEN -> ANNOTATE ->
REPORT) against a small, deterministic, purpose-built fixture -- not tests/data/, which
this repo's own nextflow.config `test` profile references but which is not actually
checked into git (a pre-existing gap, out of scope here). Skipped when `nextflow` is not
on PATH; the pixi-managed tools it needs (mmseqs2, famsa, hmmer, blast+) are activated
per-task by nextflow.config's beforeScript, not required on the test runner's own PATH.
"""
import random
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(
    shutil.which('nextflow') is None,
    reason='nextflow not on PATH -- integration test requires a real pipeline run',
)


def _seq(seed, length=110):
    """A fixed, reproducible pseudo-random amino acid sequence."""
    rng = random.Random(seed)
    return ''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=length))


def _write_fasta(path, records):
    with open(path, 'w') as fh:
        for name, seq in records.items():
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')


def _write_dna(path, n_records=2, length=300):
    records = {}
    for i in range(n_records):
        rng = random.Random(2000 + i)
        records[f'scaffold{i}'] = ''.join(rng.choices('ACGT', k=length))
    _write_fasta(path, records)


@pytest.fixture
def fixture_dir(tmp_path):
    """Build DISCOVERY_TARGET x2 / DISCOVERY_OUT x1 / NEAR_INGROUP x1 / BROAD_OUTGROUP x1
    proteomes + genomes.

    NOVEL1 is an identical sequence planted in both DISCOVERY_TARGET proteomes and in
    NEAR_INGROUP, but absent from DISCOVERY_OUT and BROAD_OUTGROUP -- a deterministic setup
    that should survive phase 1 (present in the targets, absent from DISCOVERY_OUT) and land
    in the 'clade_specific' category in phase 2 (found in NEAR_INGROUP, not in
    BROAD_OUTGROUP). DECOY* sequences are unrelated pseudo-random sequences that should
    never cluster with NOVEL1 at the default --family_min_seq_id 0.3.
    """
    data_dir = tmp_path / 'data'
    data_dir.mkdir()

    novel1 = _seq(1)
    decoy_t1 = _seq(11)
    decoy_t2 = _seq(12)
    decoy_d1 = _seq(21)
    decoy_n1 = _seq(31)
    decoy_b1 = _seq(41)

    _write_fasta(data_dir / 'target1.pep.fa', {'NOVEL1_T1': novel1, 'DECOY_T1': decoy_t1})
    _write_fasta(data_dir / 'target2.pep.fa', {'NOVEL1_T2': novel1, 'DECOY_T2': decoy_t2})
    _write_fasta(data_dir / 'discout1.pep.fa', {'DECOY_D1': decoy_d1})
    _write_fasta(data_dir / 'nearin1.pep.fa', {'NOVEL1_N1': novel1, 'DECOY_N1': decoy_n1})
    _write_fasta(data_dir / 'broadout1.pep.fa', {'DECOY_B1': decoy_b1})

    for name in ('target1', 'target2', 'discout1', 'nearin1', 'broadout1'):
        _write_dna(data_dir / f'{name}.dna.fa')

    config = tmp_path / 'config.csv'
    config.write_text(
        "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
        "DISCOVERY_TARGET,Target species one,,target1.pep.fa,target1.dna.fa,Targ1,TestClade\n"
        "DISCOVERY_TARGET,Target species two,,target2.pep.fa,target2.dna.fa,Targ2,TestClade\n"
        "DISCOVERY_OUT,Discovery outgroup one,,discout1.pep.fa,discout1.dna.fa,Disc1,OtherClade\n"
        "NEAR_INGROUP,Near ingroup one,,nearin1.pep.fa,nearin1.dna.fa,Near1,TestClade\n"
        "BROAD_OUTGROUP,Broad outgroup one,,broadout1.pep.fa,broadout1.dna.fa,Broad1,FarClade\n"
    )
    return tmp_path, data_dir, config


def test_novelty_discovery_and_screen_end_to_end(fixture_dir):
    tmp_path, data_dir, config = fixture_dir
    project = 'it_novelty_discovery'

    result = subprocess.run(
        [
            'nextflow', 'run', str(REPO / 'main.nf'),
            '--config', str(config),
            '--data_dir', str(data_dir),
            '--run_tool', 'phmmer',
            '--cluster_tool', 'novelty_discovery',
            '--max_cpus', '2',
            '--project', project,
        ],
        cwd=tmp_path, capture_output=True, text=True, timeout=900,
    )
    assert result.returncode == 0, (
        f"nextflow run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    results_dir = tmp_path / 'results' / project
    view_dir = tmp_path / 'view' / project

    for name in ('screened_presence_matrix.tsv', 'screened_candidates.txt',
                'novelties.html', 'core.html', 'losses.html'):
        assert (results_dir / name).exists(), f'missing results/{project}/{name}'
    for name in ('novelties.html', 'core.html', 'losses.html', 'report.html', 'summary.pdf'):
        assert (view_dir / name).exists(), f'missing view/{project}/{name}'

    header = (results_dir / 'screened_presence_matrix.tsv').read_text().splitlines()[0]
    fields = header.split('\t')
    assert 'novelty_category' in fields
    assert 'Near1' in fields and 'Broad1' in fields

    # NOVEL1's family: present in both DISCOVERY_TARGET genomes and NEAR_INGROUP, absent
    # from DISCOVERY_OUT/BROAD_OUTGROUP -- should survive phase 1 and land in
    # 'clade_specific'.
    cat_idx = fields.index('novelty_category')
    pid_idx = fields.index('protein_id')
    rows = [line.split('\t') for line in
           (results_dir / 'screened_presence_matrix.tsv').read_text().splitlines()[1:]]
    novel_rows = [r for r in rows if r[pid_idx].startswith('NOVEL1')]
    assert novel_rows, 'expected at least one NOVEL1_* row in the screened matrix'
    assert any(r[cat_idx] == 'clade_specific' for r in novel_rows)

    screened_candidates = (results_dir / 'screened_candidates.txt').read_text()
    assert any('NOVEL1' in line for line in screened_candidates.splitlines())
