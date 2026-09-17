"""Unit tests for bin/verify_diamond_cluster_ids.py.

Unlike mmseqs easy-cluster, `diamond cluster` is expected to keep the full
FASTA header token verbatim in its own *_cluster.tsv output -- no collapsing,
so no restoration step is needed the way bin/restore_mmseqs_cluster_ids.py
is needed for mmseqs. This script is the fail-loud safety net that verifies
that expectation actually holds for a given run, instead of assuming it
(see docs/adr/0003-diamond-tier1-clustering-backend.md).
"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'verify_diamond_cluster_ids.py'


def _run(fasta_path, cluster_tsv_path):
    return subprocess.run(
        [sys.executable, str(BIN),
         '--input-fasta', str(fasta_path),
         '--cluster-tsv', str(cluster_tsv_path)],
        capture_output=True, text=True,
    )


def test_passes_when_ids_match_full_headers_verbatim(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">Ncra|geneA desc\nMKVLA\n>Ncra|geneB desc\nMKVLB\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("Ncra|geneA\tNcra|geneA\nNcra|geneA\tNcra|geneB\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr


def test_does_not_modify_cluster_tsv(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">Ncra|geneA desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    original = "Ncra|geneA\tNcra|geneA\n"
    cluster_tsv.write_text(original)

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr
    assert cluster_tsv.read_text() == original  # verify-only, never rewritten


def test_passes_with_short_prefixed_multipipe_uniprot_style_header(tmp_path):
    """The adversarial case flagged in the ADR as never having been tested:
    a pangenome Short-prefixed header whose original id itself contained
    pipes (Short|sp|ACC|NAME). diamond is expected to preserve the whole
    thing as one atomic token, unlike mmseqs which would try to collapse it."""
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">Afum|sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("Afum|sp|O74225|YCF1_SCHPO\tAfum|sp|O74225|YCF1_SCHPO\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode == 0, r.stderr


def test_fails_loud_on_unknown_id(tmp_path):
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">Ncra|geneA desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("Ncra|geneA\tNcra|geneZZZ\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'Ncra|geneZZZ' in r.stderr


def test_fails_loud_when_id_looks_mmseqs_collapsed(tmp_path):
    """If a cluster.tsv id only matches after stripping an sp|/tr|-style
    wrapper -- i.e. it looks like diamond collapsed a header the way mmseqs
    does -- that must be a distinct, actionable failure, not a silent
    mismatch or a generic 'unknown id' message."""
    fasta = tmp_path / 'seed.faa'
    fasta.write_text(">sp|O74225|YCF1_SCHPO desc\nMKVLA\n")
    cluster_tsv = tmp_path / 'cluster.tsv'
    cluster_tsv.write_text("O74225\tO74225\n")

    r = _run(fasta, cluster_tsv)
    assert r.returncode != 0
    assert 'O74225' in r.stderr
    assert 'collapsed' in r.stderr.lower()
