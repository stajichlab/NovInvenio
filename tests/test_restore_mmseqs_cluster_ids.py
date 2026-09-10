"""Unit tests for bin/restore_mmseqs_cluster_ids.py.

mmseqs easy-cluster silently collapses certain FASTA header conventions
(sp|ACC|NAME, tr|ACC|NAME, and others -- see lib/fasta.py::mmseqs_id()) to a
shorter field in its own *_cluster.tsv output, while every other pipeline
artifact keeps the full header token as protein_id. This script corrects
*_cluster.tsv back to the full-header form immediately after mmseqs
produces it, using the same FASTA it was run on to reconstruct the mapping.
"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'restore_mmseqs_cluster_ids.py'


def _run(fasta_path, cluster_tsv_path):
    return subprocess.run(
        [sys.executable, str(BIN),
         '--input-fasta', str(fasta_path),
         '--cluster-tsv', str(cluster_tsv_path)],
        capture_output=True, text=True,
    )


def test_restores_uniprot_style_collapsed_ids(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(
        ">sp|O74225|YCF1_SCHPO Uncharacterized OS=S. pombe\nMKVLA\n"
        ">tr|A8PCG9|A8PCG9_ASPNG Putative OS=A. niger\nMKVLB\n"
    )
    cluster_tsv = tmp_path / 'cluster.tsv'
    # mmseqs's own raw output shape: bare accessions.
    cluster_tsv.write_text("O74225\tO74225\nO74225\tA8PCG9\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr
    assert cluster_tsv.read_text() == (
        "sp|O74225|YCF1_SCHPO\tsp|O74225|YCF1_SCHPO\n"
        "sp|O74225|YCF1_SCHPO\ttr|A8PCG9|A8PCG9_ASPNG\n"
    )


def test_idempotent_no_op_for_non_uniprot_headers(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">UHM102bin47__k141_1_1 desc\nMKVLA\n>UHM102bin47__k141_2_1 desc\nMKVLB\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    original = "UHM102bin47__k141_1_1\tUHM102bin47__k141_1_1\nUHM102bin47__k141_1_1\tUHM102bin47__k141_2_1\n"
    cluster_tsv.write_text(original)

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr
    assert cluster_tsv.read_text() == original  # unchanged -- already full-token


def test_running_twice_is_idempotent(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("O74225\tO74225\n")

    r1 = _run(fasta, cluster_tsv)
    assert r1.returncode == 0, r1.stderr
    once = cluster_tsv.read_text()
    r2 = _run(fasta, cluster_tsv)
    assert r2.returncode == 0, r2.stderr
    assert cluster_tsv.read_text() == once


def test_fails_loud_on_ambiguous_bare_id_collision(tmp_path):
    fasta = tmp_path / 'seed.faa'
    # Two distinct full headers that mmseqs would both collapse to "O74225".
    fasta.write_text(
        ">sp|O74225|NAME_A desc one\nMKVLA\n"
        ">tr|O74225|NAME_B desc two\nMKVLB\n"
    )
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("O74225\tO74225\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'O74225' in r.stderr


def test_fails_loud_on_unreconcilable_id(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    # An id that matches neither the full header nor its computed bare form.
    cluster_tsv.write_text("SOMETHING_ELSE\tSOMETHING_ELSE\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'SOMETHING_ELSE' in r.stderr
