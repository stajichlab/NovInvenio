import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LIB = REPO / "tests" / "data" / "uniprot_library"
sys.path.insert(0, str(REPO / "bin"))
sys.path.insert(0, str(REPO / "lib"))
import uniprot_plan_chunks as upc  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


def _rows():
    with open(LIB / "proteomes.csv") as fh:
        return list(csv.DictReader(fh))


def test_plan_chunks_groups_by_size_and_skips_missing():
    chunks = upc.plan_chunks(_rows(), LIB, chunk_bytes=1)   # 1 byte -> one proteome per chunk
    ids = [[r["proteome_id"] for r in c] for c in chunks]
    assert ids == [["UP000000001"], ["UP000000002"]]
    one = upc.plan_chunks(_rows(), LIB, chunk_bytes=10**9)
    assert len(one) == 1 and len(one[0]) == 2


def test_cli_writes_tsv_and_reports_missing(tmp_path):
    out = tmp_path / "chunks.tsv"
    p = subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_plan_chunks.py"),
                        "--library", str(LIB), "--library-csv", "proteomes.csv",
                        "--chunk-gb", "8", "--output", str(out)],
                       capture_output=True, text=True, check=True)
    assert "UP000000003" in p.stderr
    lines = out.read_text().splitlines()
    assert lines[0] == "chunk_id\tproteome_id\tfile_prefix"
    assert lines[1].startswith("chunk_0000\tUP000000001\t")


def test_parse_chunk_writes_records_and_stats(tmp_path):
    chunks = tmp_path / "chunks.tsv"
    chunks.write_text("chunk_id\tproteome_id\tfile_prefix\n"
                      "chunk_0000\tUP000000001\tUP000000001_367110\n")
    subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_parse_dat.py"),
                    "--library", str(LIB), "--chunks", str(chunks),
                    "--chunk-id", "chunk_0000", "--outdir", str(tmp_path)], check=True)
    with open_maybe_compressed(tmp_path / "records" / "UP000000001.records.tsv.zst") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["accession"] for r in rows] == ["Q7S6W2", "P00001"]
    st = json.loads((tmp_path / "stats" / "UP000000001.json").read_text())
    assert st["n_records"] == 2 and len(st["dat_sha256"]) == 64


def test_parse_chunk_unknown_chunk_id_is_an_error(tmp_path):
    chunks = tmp_path / "chunks.tsv"
    chunks.write_text("chunk_id\tproteome_id\tfile_prefix\n")
    p = subprocess.run([sys.executable, str(REPO / "bin" / "uniprot_parse_dat.py"),
                        "--library", str(LIB), "--chunks", str(chunks),
                        "--chunk-id", "chunk_0009", "--outdir", str(tmp_path)],
                       capture_output=True, text=True)
    assert p.returncode != 0 and "chunk_0009" in p.stderr
