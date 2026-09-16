#!/usr/bin/env python3
"""Fill a samplesheet's empty `TaxonGroup` cells from
pangenome_assign_clades.py's `clade_assignments.tsv` (Mash+PCoA+k-means
clade labels).

Ported/newly-written for the pangenome-profiling subworkflow: the
originating study's TaxonGroup fill (DAPC-clade > a published-table match >
Mash-clade fallback) was a manual, study-specific research step, never
automated in any script (NEXTFLOW_MIGRATION_NOTES.md section B flags this
explicitly). pangenome_cooccurrence.py reads TaxonGroup directly off the
samplesheet, not off clade_assignments.tsv -- so without SOME automated
fill step, a fresh study's co-occurrence permutation null would silently
run completely unstratified (every strain in one "unknown" clade) even
after ASSIGN_CLADES ran successfully. This script provides the generic
default the migration notes' "Open questions" section left unresolved: fill
only currently-empty TaxonGroup cells with the Mash-derived label, and
NEVER overwrite an existing (e.g. hand-curated DAPC) value. A study that
wants a richer priority scheme (DAPC > published-table match > Mash
fallback, as the originating study did) still applies its own fill before
this subworkflow runs -- this script only fills what a study left blank.

Preserves every other column and the original column order exactly (a
csv.DictReader/DictWriter round-trip), so this is safe to use as this
subworkflow's samplesheet for every downstream step without any schema
drift.

Usage:
  pangenome_fill_taxon_group.py --config config.csv \\
      --clade_assignments clade_assignments.tsv --output config.with_clades.csv
"""
from __future__ import annotations

import argparse
import csv
import sys


def read_clade_labels(path: str) -> dict[str, str]:
    """{Short: mash_clade} from pangenome_assign_clades.py's output."""
    labels: dict[str, str] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            short = row.get("Short", "").strip()
            clade = row.get("mash_clade", "").strip()
            if short and clade:
                labels[short] = clade
    return labels


def fill_taxon_group(config_path: str, clade_labels: dict[str, str], output_path: str) -> tuple[int, int]:
    """Returns (n_filled, n_total_rows)."""
    with open(config_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if fieldnames is None or "TaxonGroup" not in fieldnames or "Short" not in fieldnames:
        print(
            f"ERROR: {config_path} has no 'Short'/'TaxonGroup' column -- "
            f"got {fieldnames}", file=sys.stderr,
        )
        sys.exit(1)

    n_filled = 0
    for row in rows:
        if not (row.get("TaxonGroup") or "").strip():
            clade = clade_labels.get(row["Short"])
            if clade:
                row["TaxonGroup"] = clade
                n_filled += 1

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return n_filled, len(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--clade_assignments", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    clade_labels = read_clade_labels(args.clade_assignments)
    n_filled, n_total = fill_taxon_group(args.config, clade_labels, args.output)
    print(
        f"pangenome_fill_taxon_group: filled {n_filled}/{n_total} rows' empty "
        "TaxonGroup from clade_assignments.tsv (existing values kept as-is)",
        file=sys.stderr,
    )
    if n_filled == 0 and n_total > 0:
        print(
            "WARNING: pangenome_fill_taxon_group filled ZERO rows -- this "
            "usually means something upstream is broken (e.g. "
            "clade_assignments.tsv is empty/malformed, every samplesheet row "
            "already had a TaxonGroup, or the Short IDs in "
            "clade_assignments.tsv don't match the samplesheet's Short "
            "column). Verify this is expected before trusting downstream "
            "frequency/co-occurrence results that key on TaxonGroup.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
