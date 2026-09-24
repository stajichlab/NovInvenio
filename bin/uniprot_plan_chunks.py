#!/usr/bin/env python3
"""Group a UniProt library's proteomes into parse batches of about --chunk-gb of
compressed .dat.gz each (toward 1-1.5 h per UNIPROT_PARSE_CHUNK job; the parser
reads ~2.8 MB/s of .dat.gz with all record fields, measured 2026-09-23 on
UP000001805). Proteomes listed in the CSV with no .dat.gz are reported on stderr
and left out. See docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md."""
import argparse
import csv
import sys
from pathlib import Path


def plan_chunks(rows, library, chunk_bytes):
    """rows: library-CSV dicts (proteome_id, file_prefix, ...) -> list of batches,
    each a list of rows whose .dat.gz sizes sum to about chunk_bytes."""
    chunks, cur, cur_size = [], [], 0
    for r in rows:
        f = Path(library) / "data" / f"{r['file_prefix']}.dat.gz"
        if not f.exists():
            continue
        size = f.stat().st_size
        if cur and cur_size + size > chunk_bytes:
            chunks.append(cur)
            cur, cur_size = [], 0
        cur.append(r)
        cur_size += size
    if cur:
        chunks.append(cur)
    return chunks


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--library-csv", required=True, help="CSV file name inside --library")
    ap.add_argument("--chunk-gb", type=float, default=8.0)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    with open(a.library / a.library_csv) as fh:
        rows = list(csv.DictReader(fh))
    missing = [r["proteome_id"] for r in rows
               if not (a.library / "data" / f"{r['file_prefix']}.dat.gz").exists()]
    for pid in missing:
        print(f"WARNING: no .dat.gz for {pid}; skipped", file=sys.stderr)
    chunks = plan_chunks(rows, a.library, int(a.chunk_gb * 1e9))
    with open(a.output, "w") as out:
        out.write("chunk_id\tproteome_id\tfile_prefix\n")
        for i, c in enumerate(chunks):
            for r in c:
                out.write(f"chunk_{i:04d}\t{r['proteome_id']}\t{r['file_prefix']}\n")
    print(f"{len(chunks)} chunks, {sum(map(len, chunks))} proteomes, "
          f"{len(missing)} missing", file=sys.stderr)


if __name__ == "__main__":
    main()
