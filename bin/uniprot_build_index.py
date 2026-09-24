#!/usr/bin/env python3
"""Merge records/<UP>.records.tsv.zst + stats/<UP>.json (bin/uniprot_parse_dat.py)
into a UniProt library index: seq_index.sqlite + manifest.json + records/
(format lib/uniprot_index.FORMAT_VERSION). manifest.json is written last, so an
interrupted build fails lib/uniprot_index.UniProtIndex's completeness check."""
import argparse
import csv
import datetime
import hashlib
import json
import re
import shutil
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from uniprot_index import FORMAT_VERSION  # noqa: E402

RELEASE_RE = re.compile(r"Release\s+(\d{4}_\d{2})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--records-dir", required=True, type=Path)
    ap.add_argument("--stats-dir", required=True, type=Path)
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--library-csv", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = a.output_dir / "manifest.json"
    if manifest_path.exists():
        manifest_path.unlink()
    csv_path = a.library / a.library_csv
    with open(csv_path) as fh:
        lib_rows = list(csv.DictReader(fh))
    readme = a.library / "README"
    m = RELEASE_RE.search(readme.read_text(errors="replace")) if readme.exists() else None

    db = a.output_dir / "seq_index.sqlite"
    if db.exists():
        db.unlink()
    con = sqlite3.connect(db)
    con.executescript("""
        create table seq(md5 text, accession text, proteome_id text, taxid integer, reviewed integer);
        create table acc(accession text primary key, proteome_id text, taxid integer);
        create table refseq(refseq_id text, accession text);
        create table proteome(proteome_id text primary key, taxid integer, species_name text,
                              file_prefix text);
    """)
    stats = {p.stem: json.loads(p.read_text()) for p in a.stats_dir.glob("*.json")}
    n_records = 0
    for row in lib_rows:
        pid = row["proteome_id"]
        if pid not in stats:
            continue
        con.execute("insert into proteome values (?,?,?,?)",
                    (pid, int(row["tax_id"]), row["species_name"], row["file_prefix"]))
        seq_rows, acc_rows, ref_rows = [], [], []
        with open_maybe_compressed(a.records_dir / f"{pid}.records.tsv.zst") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tax = int(r["taxon_id"] or row["tax_id"])
                seq_rows.append((r["seq_md5"], r["accession"], pid, tax, int(r["reviewed"])))
                acc_rows.append((r["accession"], pid, tax))
                for x in r["xrefs"].split("|"):
                    if x.startswith("RefSeq:"):
                        ref_rows.append((x.split(":", 1)[1], r["accession"]))
        con.executemany("insert into seq values (?,?,?,?,?)", seq_rows)
        con.executemany("insert or ignore into acc values (?,?,?)", acc_rows)
        con.executemany("insert into refseq values (?,?)", ref_rows)
        con.commit()
        n_records += len(seq_rows)
    con.executescript("create index seq_md5 on seq(md5); "
                      "create index refseq_id on refseq(refseq_id);")
    con.commit()
    con.close()

    out_records = a.output_dir / "records"
    if a.records_dir.resolve() != out_records.resolve():
        out_records.mkdir(exist_ok=True)
        for p in a.records_dir.glob("*.records.tsv.zst"):
            shutil.copyfile(p, out_records / p.name)

    missing = [r["proteome_id"] for r in lib_rows if r["proteome_id"] not in stats]
    manifest = {
        "format_version": FORMAT_VERSION,
        "library": str(a.library),
        "release": m.group(1) if m else "",
        "library_csv": a.library_csv,
        "library_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "n_proteomes": len(stats),
        "n_records": n_records,
        "missing": missing,
        "proteomes": [stats[k] for k in sorted(stats)],
        "build_date": datetime.date.today().isoformat(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=1))
    print(f"index: {len(stats)} proteomes, {n_records} records, {len(missing)} missing",
          file=sys.stderr)


if __name__ == "__main__":
    main()
