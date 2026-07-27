"""Unit tests for bin/extract_singletons.py."""
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / 'bin' / 'extract_singletons.py'

CLUSTER_TSV = "repA\trepA\nrepB\tm1\nrepB\tm2\nrepC\trepC\n"
FASTA = ">repA\nMKL\n>m1\nMNP\n>m2\nMQR\n>repC\nMST\n"


def _run(tmp_path):
    ctsv = tmp_path / 'cluster.tsv'
    ctsv.write_text(CLUSTER_TSV)
    fa = tmp_path / 'seed.fa'
    fa.write_text(FASTA)
    out = tmp_path / 'singletons.fa'
    subprocess.run(
        [sys.executable, str(BIN),
         '--cluster-tsv', str(ctsv),
         '--fasta', str(fa),
         '--output', str(out)],
        check=True,
    )
    return out.read_text()


def test_extracts_singleton_clusters(tmp_path):
    result = _run(tmp_path)
    # repA and repC are singletons (rep == member, only 1 member)
    assert '>repA' in result
    assert '>repC' in result
    # repB is a multi-member family → NOT a singleton
    assert '>repB' not in result
    # m1 and m2 are members of repB's family, not singletons
    assert '>m1' not in result
    assert '>m2' not in result


def test_singleton_sequences_match_fasta(tmp_path):
    result = _run(tmp_path)
    lines = result.strip().split('\n')
    # repA should have its sequence MKL
    repa_idx = lines.index('>repA')
    assert lines[repa_idx + 1] == 'MKL'
