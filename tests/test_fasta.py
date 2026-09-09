"""Unit tests for lib/fasta.py's read_fasta() (2026-09-03 fix: tolerate duplicate IDs
instead of crashing, a real occurrence when concatenating independently-sourced
proteomes -- e.g. deep_broad_1kfg's PROFILE_LOSS_SEARCH family building hit a real
cross-genome ID collision, 'Pchr|PCH_Pc18g05710.1', that used to abort the whole run)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from fasta import mmseqs_id, read_fasta  # noqa: E402


def test_no_duplicates_reads_normally(tmp_path):
    fa = tmp_path / 'x.fasta'
    fa.write_text(">a\nMAAA\n>b\nMBBB\n")
    records = read_fasta(fa)
    assert set(records) == {'a', 'b'}
    assert str(records['a'].seq) == 'MAAA'


def test_duplicate_id_keeps_first_and_warns(tmp_path, capsys):
    fa = tmp_path / 'x.fasta'
    fa.write_text(">shared\nMFIRST\n>b\nMBBB\n>shared\nMSECOND\n")
    records = read_fasta(fa)
    assert set(records) == {'shared', 'b'}
    assert str(records['shared'].seq) == 'MFIRST'  # first occurrence kept
    captured = capsys.readouterr()
    assert 'duplicate' in captured.err.lower()
    assert '1 duplicate' in captured.err


def test_mmseqs_id_extracts_uniprot_accession():
    # Empirically confirmed against a real mmseqs2 easy-cluster run (issue #85):
    # both "sp|" (reviewed) and "tr|" (unreviewed) UniProt deflines report the
    # middle field as the sequence's cluster.tsv id.
    assert mmseqs_id('sp|O74225|YCF1_SCHPO') == 'O74225'
    assert mmseqs_id('tr|A8PCG9|A8PCG9_ASPNG') == 'A8PCG9'


def test_mmseqs_id_is_a_no_op_for_non_uniprot_headers():
    # BFD/FungiDB-style bare protein IDs (no recognized prefix) pass through
    # unchanged -- this is why the pairwise (--cluster_tool pairwise) path
    # never surfaced this mismatch.
    assert mmseqs_id('Ncra_00001') == 'Ncra_00001'
    assert mmseqs_id('Pchr|PCH_Pc18g05710.1') == 'Pchr|PCH_Pc18g05710.1'
