"""Compression support for the large pangenome TSV intermediates (issue #111).

`lib/compressed_io.py` has provided `open_maybe_compressed_write()` since it
was added, but nothing called it -- every pangenome output except
RESCUE_PASS's tblastn published as plain text, ~4.1 GB per 529-strain run
across `pair_classification.tsv`, `cooccurring_pairs.tsv`,
`family_positions.tsv` and `gene_positions.tsv`.

These tests pin both halves of the contract: a producer writes real
compressed output when its `--output` is named `.zst`, and every consumer of
those four files reads a compressed file exactly as it reads a plain one.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).parent.parent / "bin"
LIB = Path(__file__).parent.parent / "lib"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(LIB))

from compressed_io import open_maybe_compressed, open_maybe_compressed_write  # noqa: E402

needs_zstd = pytest.mark.skipif(
    shutil.which("zstd") is None, reason="zstd CLI not available"
)

FAMILY_POSITIONS = (
    "Short\tfamily\tcontig\trank\n"
    "s1\tfamA\tc1\t1\n"
    "s1\tfamB\tc1\t2\n"
    "s2\tfamA\tc9\t5\n"
)

PAIR_CLASSIFICATION = (
    "family_a\tfamily_b\tclassification\n"
    "famA\tfamB\tunexplained_physical\n"
    "famA\tfamC\ttrans\n"
)

GENE_POSITIONS = (
    "Short\tprotein_id\tcontig\tstart\tend\n"
    "s1\tp1\tc1\t100\t400\n"
    "s1\tp2\tc1\t500\t900\n"
    "s2\tp3\tc9\t10\t60\n"
)


def write_zst(path: Path, text: str) -> Path:
    """Write `text` to `path` as a real zstd file, via the CLI rather than
    the helper under test -- a round-trip through one implementation would
    not catch that implementation being wrong in both directions."""
    plain = path.with_suffix(path.suffix + ".plain")
    plain.write_text(text)
    subprocess.run(["zstd", "-q", "-f", str(plain), "-o", str(path)], check=True)
    return path


# --- lib/compressed_io.py itself (previously untested) -----------------

@pytest.mark.parametrize("suffix", ["", ".gz", ".zst"])
def test_compressed_io_write_then_read_round_trips(tmp_path, suffix):
    if suffix == ".zst" and shutil.which("zstd") is None:
        pytest.skip("zstd CLI not available")
    path = tmp_path / f"data.tsv{suffix}"
    payload = "col\n" + "".join(f"row{i}\n" for i in range(500))
    with open_maybe_compressed_write(path) as fh:
        fh.write(payload)
    with open_maybe_compressed(path) as fh:
        assert fh.read() == payload


@needs_zstd
def test_compressed_io_write_zst_produces_a_real_zstd_file(tmp_path):
    # Guards against the suffix being honoured in name only -- a plain-text
    # file called .zst would still round-trip through a buggy reader.
    path = tmp_path / "data.tsv.zst"
    with open_maybe_compressed_write(path) as fh:
        fh.write("col\nvalue\n")
    assert path.read_bytes()[:4] == b"\x28\xb5\x2f\xfd"  # zstd magic number


# --- consumers of pair_classification.tsv / family_positions.tsv -------

@needs_zstd
def test_build_islands_reads_zst_family_positions(tmp_path):
    from pangenome_build_islands import load_strain_gene_orders

    plain = tmp_path / "family_positions.tsv"
    plain.write_text(FAMILY_POSITIONS)
    zst = write_zst(tmp_path / "family_positions.tsv.zst", FAMILY_POSITIONS)
    assert load_strain_gene_orders(str(zst)) == load_strain_gene_orders(str(plain))
    assert load_strain_gene_orders(str(zst))["s1"][0][0] == "famA"


@needs_zstd
def test_build_islands_reads_zst_pair_classification(tmp_path):
    from pangenome_build_islands import load_significant_physical_pairs

    plain = tmp_path / "pair_classification.tsv"
    plain.write_text(PAIR_CLASSIFICATION)
    zst = write_zst(tmp_path / "pair_classification.tsv.zst", PAIR_CLASSIFICATION)
    result = load_significant_physical_pairs(str(zst))
    assert result == load_significant_physical_pairs(str(plain))
    assert result == {frozenset({"famA", "famB"}): "unexplained_physical"}


@needs_zstd
def test_report_tables_reads_zst_pair_classification(tmp_path):
    from pangenome_report_tables import classification_counts

    plain = tmp_path / "pair_classification.tsv"
    plain.write_text(PAIR_CLASSIFICATION)
    zst = write_zst(tmp_path / "pair_classification.tsv.zst", PAIR_CLASSIFICATION)
    assert classification_counts(str(zst)) == classification_counts(str(plain))
    assert classification_counts(str(zst)) == {"unexplained_physical": 1, "trans": 1}


# --- producer round-trip, across two scripts ---------------------------

@needs_zstd
def test_build_family_positions_reads_zst_input_and_writes_zst_output(tmp_path):
    # The real chain: BUILD_GENE_POSITIONS writes gene_positions.tsv.zst,
    # BUILD_FAMILY_POSITIONS consumes it and writes family_positions.tsv.zst,
    # BUILD_ISLANDS consumes that. Exercise all three links at once.
    from pangenome_build_islands import load_strain_gene_orders

    gene_positions = write_zst(tmp_path / "gene_positions.tsv.zst", GENE_POSITIONS)
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("famA\ts1|p1\nfamA\ts2|p3\nfamB\ts1|p2\n")
    out = tmp_path / "family_positions.tsv.zst"

    subprocess.run(
        [sys.executable, str(BIN / "pangenome_build_family_positions.py"),
         "--gene_positions", str(gene_positions),
         "--cluster_tsv", str(cluster_tsv),
         "--output", str(out)],
        check=True,
    )

    assert out.read_bytes()[:4] == b"\x28\xb5\x2f\xfd"
    by_strain = load_strain_gene_orders(str(out))
    assert set(by_strain) == {"s1", "s2"}
    assert {f for rows in by_strain.values() for f, *_ in rows} == {"famA", "famB"}


# --- defect fix: missing or corrupt compressed files should fail loudly ------

def test_open_maybe_compressed_missing_plain_raises_filenotfound(tmp_path):
    """Missing plain file raises FileNotFoundError (existing correct behavior)."""
    path = tmp_path / "nonexistent.tsv"
    with pytest.raises(FileNotFoundError):
        with open_maybe_compressed(str(path)) as fh:
            fh.read()


@needs_zstd
def test_open_maybe_compressed_missing_zst_raises_filenotfound(tmp_path):
    """Missing .zst file raises FileNotFoundError, not silent empty stream."""
    path = tmp_path / "nonexistent.tsv.zst"
    with pytest.raises(FileNotFoundError):
        with open_maybe_compressed(str(path)) as fh:
            fh.read()


@needs_zstd
def test_open_maybe_compressed_corrupt_zst_raises(tmp_path):
    """Corrupt .zst file (invalid header) raises, not silent empty stream."""
    import os
    path = tmp_path / "corrupt.tsv.zst"
    path.write_bytes(os.urandom(200))  # Random garbage, not valid zstd
    with pytest.raises(RuntimeError) as exc_info:
        with open_maybe_compressed(str(path)) as fh:
            fh.read()
    assert "zstd" in str(exc_info.value).lower()


def test_open_maybe_compressed_missing_gz_raises_filenotfound(tmp_path):
    """.gz file missing raises FileNotFoundError, not silent empty stream."""
    path = tmp_path / "nonexistent.tsv.gz"
    with pytest.raises(FileNotFoundError):
        with open_maybe_compressed(str(path)) as fh:
            fh.read()


@needs_zstd
def test_open_maybe_compressed_valid_zst_still_reads_correctly(tmp_path):
    """Valid .zst files still decompress and read correctly."""
    path = write_zst(tmp_path / "valid.tsv.zst", "col\nvalue1\nvalue2\n")
    with open_maybe_compressed(str(path)) as fh:
        lines = fh.readlines()
    assert lines == ["col\n", "value1\n", "value2\n"]


@needs_zstd
def test_open_maybe_compressed_symlink_to_zst_still_works(tmp_path):
    """Symlink to valid .zst file (Nextflow staging) still works correctly."""
    real_zst = write_zst(tmp_path / "real.tsv.zst", "a\tb\n1\t2\n")
    link = tmp_path / "link.tsv.zst"
    link.symlink_to(real_zst)
    with open_maybe_compressed(str(link)) as fh:
        content = fh.read()
    assert content == "a\tb\n1\t2\n"
