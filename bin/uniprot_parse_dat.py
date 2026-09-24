#!/usr/bin/env python3
"""Parse one batch (--chunk-id) of a UniProt library's .dat.gz files into
records/<proteome_id>.records.tsv.zst + stats/<proteome_id>.json. See
docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed_write  # noqa: E402
from uniprot_dat import RECORD_COLUMNS, parse_dat_gz  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--chunks", required=True, type=Path)
    ap.add_argument("--chunk-id", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args()
    with open(a.chunks) as fh:
        todo = [r for r in csv.DictReader(fh, delimiter="\t") if r["chunk_id"] == a.chunk_id]
    if not todo:
        sys.exit(f"ERROR: no proteomes for {a.chunk_id} in {a.chunks}")
    (a.outdir / "records").mkdir(parents=True, exist_ok=True)
    (a.outdir / "stats").mkdir(parents=True, exist_ok=True)
    for r in todo:
        dat = a.library / "data" / f"{r['file_prefix']}.dat.gz"
        stats = {}
        out = a.outdir / "records" / f"{r['proteome_id']}.records.tsv.zst"
        with open_maybe_compressed_write(out) as fh:
            w = csv.DictWriter(fh, fieldnames=RECORD_COLUMNS, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            for rec in parse_dat_gz(dat, stats=stats):
                w.writerow(rec)
        stats.update(proteome_id=r["proteome_id"], file_prefix=r["file_prefix"],
                     dat_sha256=sha256(dat))
        (a.outdir / "stats" / f"{r['proteome_id']}.json").write_text(json.dumps(stats))
        print(f"{r['proteome_id']}: {stats['n_records']} records "
              f"({stats['n_skipped_no_ac']} no AC, {stats['n_skipped_no_seq']} no SQ)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
