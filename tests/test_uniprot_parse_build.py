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
