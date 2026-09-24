import hashlib
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
LIB = REPO / "tests" / "data" / "uniprot_library"
sys.path.insert(0, str(REPO / "lib"))
from uniprot_index import FORMAT_VERSION, IndexFormatError, UniProtIndex  # noqa: E402


def build_fixture_index(tmp_path):
    """Plan, parse and build the fixture library into tmp_path/index."""
    work = tmp_path / "work"

    def run(*a):
        subprocess.run([sys.executable, *map(str, a)], check=True, capture_output=True)

    run(REPO / "bin" / "uniprot_plan_chunks.py", "--library", LIB,
        "--library-csv", "proteomes.csv", "--chunk-gb", "8", "--output", tmp_path / "chunks.tsv")
    run(REPO / "bin" / "uniprot_parse_dat.py", "--library", LIB, "--chunks",
        tmp_path / "chunks.tsv", "--chunk-id", "chunk_0000", "--outdir", work)
    run(REPO / "bin" / "uniprot_build_index.py", "--records-dir", work / "records",
        "--stats-dir", work / "stats", "--library", LIB, "--library-csv", "proteomes.csv",
        "--output-dir", tmp_path / "index")
    return tmp_path / "index"


def test_manifest_and_tables(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    m = json.loads((idx_dir / "manifest.json").read_text())
    assert m["format_version"] == FORMAT_VERSION
    assert m["release"] == "2099_01"
    assert m["n_proteomes"] == 2 and m["n_records"] == 3
    assert m["missing"] == ["UP000000003"]
    con = sqlite3.connect(idx_dir / "seq_index.sqlite")
    assert con.execute("select count(*) from seq").fetchone()[0] == 3
    assert con.execute(
        "select accession from refseq where refseq_id='XP_960565.1'").fetchone()[0] == "Q7S6W2"


def test_reader_lookups(tmp_path):
    idx = UniProtIndex(build_fixture_index(tmp_path))
    assert idx.by_accession("Q7S6W2")["proteome_id"] == "UP000000001"
    assert idx.by_accession("NOPE") is None
    assert idx.by_refseq("XP_960565.1") == ["Q7S6W2"]
    hits = idx.by_md5(hashlib.md5(b"MKVLLAQ").hexdigest())
    assert sorted(h["accession"] for h in hits) == ["P00002", "Q7S6W2"]
    assert idx.taxid_for_species("Neurospora crassa") == 367110
    assert idx.taxid_for_species("Neurospora crassa OR74A") == 367110
    assert idx.taxid_for_species("Nonexistent fungus") is None
    assert idx.taxids_for_species("Neurospora crassa") == {367110}
    assert idx.taxids_for_species("Nonexistent fungus") == set()
    assert "Saccharomyces cerevisiae" in idx.species_name("UP000000002")
    recs = idx.records("UP000000001", {"Q7S6W2"})
    assert list(recs) == ["Q7S6W2"] and recs["Q7S6W2"]["alphafold_id"] == "Q7S6W2"


def test_incomplete_index_is_a_hard_error(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    (idx_dir / "manifest.json").unlink()
    with pytest.raises(IndexFormatError, match="manifest.json"):
        UniProtIndex(idx_dir)


def test_wrong_format_version_is_a_hard_error(tmp_path):
    idx_dir = build_fixture_index(tmp_path)
    m = json.loads((idx_dir / "manifest.json").read_text())
    m["format_version"] = 99
    (idx_dir / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(IndexFormatError, match="format_version 99"):
        UniProtIndex(idx_dir)
