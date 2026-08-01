"""Unit tests for bin/extend_singleton_query.py (issue #52)."""
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / 'bin' / 'extend_singleton_query.py'

PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def _fasta(records):
    return ''.join(f'>{sid}\n{seq}\n' for sid, seq in records.items())


def _read_ids(path):
    return {line[1:].strip() for line in path.read_text().splitlines() if line.startswith('>')}


def run(tmp_path, singletons, seed, paralog_text=None):
    singletons_fa = tmp_path / 'singletons.fa'
    singletons_fa.write_text(_fasta(singletons))
    seed_fa = tmp_path / 'seed.fa'
    seed_fa.write_text(_fasta(seed))
    output = tmp_path / 'singleton_query.fa'
    cmd = [sys.executable, str(SCRIPT),
           '--singletons-fa', str(singletons_fa),
           '--seed-fasta', str(seed_fa),
           '--output', str(output)]
    if paralog_text is not None:
        paralog_path = tmp_path / 'paralog_cutoffs.tsv'
        paralog_path.write_text(PARALOG_HEADER + paralog_text)
        cmd += ['--paralog-cutoffs', str(paralog_path)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return _read_ids(output)


def test_singleton_with_no_paralog_passes_through_unchanged(tmp_path):
    ids = run(tmp_path, {'pS1': 'MKV'}, {'pS1': 'MKV', 'pOther': 'AAA'})
    assert ids == {'pS1'}


def test_singleton_paralog_is_added_from_seed_fasta(tmp_path):
    ids = run(tmp_path, {'pS1': 'MKV'}, {'pS1': 'MKV', 'pEif1': 'QQQ'},
             paralog_text='pS1\tpEif1\t42\t4.2e-11\n')
    assert ids == {'pS1', 'pEif1'}


def test_paralog_already_a_singleton_is_not_duplicated(tmp_path):
    # pEif1 is itself also a singleton (in singletons.fa already) -- must appear once.
    ids = run(tmp_path, {'pS1': 'MKV', 'pEif1': 'QQQ'}, {'pS1': 'MKV', 'pEif1': 'QQQ'},
             paralog_text='pS1\tpEif1\t42\t4.2e-11\n')
    assert ids == {'pS1', 'pEif1'}


def test_no_singletons_produces_empty_output(tmp_path):
    ids = run(tmp_path, {}, {'pOther': 'AAA'})
    assert ids == set()
