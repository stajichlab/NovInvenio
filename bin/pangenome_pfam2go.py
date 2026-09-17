#!/usr/bin/env python3
"""GO-term annotation of Pfam domains already found by this pipeline's own
hmmscan (island_pfam_enrichment.tsv) via a standard pfam2go mapping file
(http://current.geneontology.org/ontology/external2go/pfam2go). Ported from
studies/fungi/Afumigatus_pangenome/bin/map_pfam_to_go.py -- annotates only
domains already discovered upstream, does not run a fresh InterProScan.

Usage:
  pangenome_pfam2go.py --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --pfam2go pfam2go --output island_pfam_enrichment.go.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import bare_pfam_accession  # noqa: E402

PFAM2GO_LINE_RE = re.compile(r"^Pfam:(PF\d+)\s.*?;\s*(GO:\d+)\s*$")


def parse_pfam2go(path: str) -> dict[str, list[str]]:
    """{bare_pfam_accession: [go_id, ...]} from a standard pfam2go file.
    One line per Pfam-accession/GO-term pair -- a domain with N GO terms
    appears on N separate lines."""
    mapping: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            m = PFAM2GO_LINE_RE.match(line.strip())
            if not m:
                continue
            accession, go_id = m.group(1), m.group(2)
            mapping.setdefault(accession, []).append(go_id)
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--pfam2go", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    pfam2go = parse_pfam2go(args.pfam2go)
    print(f"pangenome_pfam2go: {len(pfam2go)} Pfam accessions with >=1 GO term in {args.pfam2go}",
          file=sys.stderr)

    with open(args.island_pfam_enrichment, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or []) + ["n_go_terms", "go_terms"]

    n_annotated = 0
    for row in rows:
        go_ids = pfam2go.get(bare_pfam_accession(row.get("pfam_accession", "-")), [])
        row["n_go_terms"] = str(len(go_ids))
        row["go_terms"] = ";".join(go_ids) if go_ids else "-"
        if go_ids:
            n_annotated += 1

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"pangenome_pfam2go: {n_annotated}/{len(rows)} domains got >=1 GO term",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
