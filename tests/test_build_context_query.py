"""Unit tests for bin/build_context_query.py (issue #48)."""
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / 'bin' / 'build_context_query.py'

PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def run(tmp_path, candidates_text, paralog_text=None):
    candidates_path = tmp_path / 'candidates.txt'
    candidates_path.write_text(candidates_text)
    cmd = [sys.executable, str(SCRIPT),
           '--candidates', str(candidates_path),
           '--output', str(tmp_path / 'context_query.txt')]
    if paralog_text is not None:
        paralog_path = tmp_path / 'paralog_cutoffs.tsv'
        paralog_path.write_text(PARALOG_HEADER + paralog_text)
        cmd += ['--paralog-cutoffs', str(paralog_path)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    out = tmp_path / 'context_query.txt'
    return [line for line in out.read_text().splitlines() if line] if out.stat().st_size else []


def test_candidate_with_no_paralog_passes_through_unchanged(tmp_path):
    entries = run(tmp_path, 'In1::g1\n')
    assert entries == ['In1::g1']


def test_candidate_paralog_is_added_from_same_proteome(tmp_path):
    entries = run(tmp_path, 'In1::hexA\n',
                  paralog_text='hexA\teif1\t42\t4.2e-11\n')
    assert set(entries) == {'In1::hexA', 'In1::eif1'}


def test_duplicate_entries_are_deduplicated(tmp_path):
    # eif1 is itself a candidate too, and also hexA's paralog -- must appear once.
    entries = run(tmp_path, 'In1::hexA\nIn1::eif1\n',
                  paralog_text='hexA\teif1\t42\t4.2e-11\n')
    assert sorted(entries) == ['In1::eif1', 'In1::hexA']


def test_empty_candidates_produces_empty_output(tmp_path):
    entries = run(tmp_path, '')
    assert entries == []
