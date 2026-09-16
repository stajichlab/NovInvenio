#!/usr/bin/env python3
"""Filter tier1_rep_seq.fasta down to the families actually eligible for
Pfam-domain enrichment testing: shell+cloud bins only, matching
pangenome_cooccurrence.py's own co-occurrence-eligibility selection (never
core/soft_core/singleton -- core was never eligible for co-occurrence
testing in the first place, and a whole-genome background would spuriously
enrich for "accessory-genome-typical" domains regardless of which specific
island is tested).

Usage:
  pangenome_select_background_reps.py --rep_fasta tier1_rep_seq.fasta \
      --frequency_table frequency_table.tsv --output background_reps.fa
"""
from __future__ import annotations

import argparse
import csv
import sys


def select_background_families(frequency_table_path: str) -> set[str]:
    """All families in the shell or cloud bin."""
    background = set()
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["bin"] in ("shell", "cloud"):
                background.add(row["family"])
    return background


def write_background_fasta(rep_fasta_path: str, background_families: set[str], out_path: str) -> int:
    """Writes only the FASTA records whose header ID is in
    `background_families`. Returns the number of records written."""
    n_written = 0
    writing = False
    with open(rep_fasta_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                family_id = line[1:].split(None, 1)[0].rstrip()
                writing = family_id in background_families
                if writing:
                    n_written += 1
            if writing:
                fout.write(line)
    return n_written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep_fasta", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    background = select_background_families(args.frequency_table)
    n = write_background_fasta(args.rep_fasta, background, args.output)
    print(f"pangenome_select_background_reps: {n} shell+cloud family reps written to {args.output}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
