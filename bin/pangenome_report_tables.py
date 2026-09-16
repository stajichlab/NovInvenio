#!/usr/bin/env python3
"""Tidy, matplotlib-free aggregation tables for the pangenome island+Pfam
report -- joins significant_islands.tsv to island_pfam_enrichment.tsv,
computes size distribution, classification counts, and per-strain gene/bin
summaries. Deliberately kept free of any plotting dependency so these
tables are independently useful (e.g. to this repo's docs/ publishing
pipeline) and independently testable. Generalizes
studies/fungi/Afumigatus_pangenome/bin/annotate_islands_with_enrichment.py's
join logic, plus a new per_strain_summary (Opus review: this is how that
study's own accessory-count outlier strain was found, by hand -- generating
it automatically here).

Usage:
  pangenome_report_tables.py --significant_islands significant_islands.tsv \\
      --island_pfam_enrichment island_pfam_enrichment.tsv \\
      --pair_classification pair_classification.tsv \\
      --presence_matrix presence_matrix.tsv \\
      --frequency_table frequency_table.tsv \\
      --domtblout background_reps_vs_pfam.domtblout \\
      --out_dir report_tables/
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pangenome_domain_enrichment import parse_domtblout  # noqa: E402


def annotate_islands_with_domains(islands_rows: list[dict], family_domains: dict[str, set[str]]) -> list[dict]:
    """Adds a `pfam_domains` key (comma-joined sorted domain names, or `-`
    if none) to each island row, from the union of its member families'
    domains."""
    out = []
    for row in islands_rows:
        row = dict(row)
        members = row["member_families"].split(",")
        domains: set[str] = set()
        for m in members:
            domains |= family_domains.get(m, set())
        row["pfam_domains"] = ",".join(sorted(domains)) if domains else "-"
        out.append(row)
    return out


def island_size_distribution(islands_rows: list[dict]) -> dict[int, int]:
    """{island_size: count}."""
    dist: dict[int, int] = {}
    for row in islands_rows:
        size = int(row["island_size"])
        dist[size] = dist.get(size, 0) + 1
    return dist


def classification_counts(pair_classification_path: str) -> dict[str, int]:
    """{classification: count} across all rows in pair_classification.tsv."""
    counts: dict[str, int] = {}
    with open(pair_classification_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c = row["classification"]
            counts[c] = counts.get(c, 0) + 1
    return counts


def per_strain_summary(presence_matrix_path: str, family_bin: dict[str, str]) -> list[dict]:
    """One row per strain: total genes present, plus counts broken down by
    frequency bin (core/soft_core/shell/cloud/singleton). `genome_only`
    counts as present (matches lib/pangenome_matrix.PresenceMatrix.is_present's
    semantics) -- rescued genome-only calls are real presence evidence, not
    a weaker state."""
    with open(presence_matrix_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        strains = header[1:]
        totals = {
            s: {"Short": s, "n_genes": 0, "core": 0, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0}
            for s in strains
        }
        for row in reader:
            family = row[0]
            b = family_bin.get(family)
            for strain, call in zip(strains, row[1:]):
                if call != "absent":
                    totals[strain]["n_genes"] += 1
                    if b in totals[strain]:
                        totals[strain][b] += 1
    return list(totals.values())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_islands", required=True)
    ap.add_argument("--island_pfam_enrichment", required=True)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--presence_matrix", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--domtblout", required=True, action="append")
    ap.add_argument("--domain_evalue", type=float, default=1e-3,
                     help="Domain-level i-Evalue cutoff for Pfam domain hits (default: 1e-3).")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.significant_islands, newline="") as fh:
        islands_rows = list(csv.DictReader(fh, delimiter="\t"))

    family_domains = parse_domtblout(args.domtblout, max_ievalue=args.domain_evalue)
    annotated = annotate_islands_with_domains(islands_rows, family_domains)
    with open(out_dir / "islands_with_domains.tsv", "w", newline="") as out:
        fieldnames = list(annotated[0].keys()) if annotated else [
            "n_strains", "example_strain", "island_size", "member_families",
            "n_supporting_pairs", "classifications", "pfam_domains",
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(annotated, key=lambda r: -int(r["island_size"])):
            writer.writerow(row)

    dist = island_size_distribution(islands_rows)
    with open(out_dir / "island_size_distribution.tsv", "w") as out:
        out.write("island_size\tcount\n")
        for size, count in sorted(dist.items()):
            out.write(f"{size}\t{count}\n")

    counts = classification_counts(args.pair_classification)
    with open(out_dir / "classification_counts.tsv", "w") as out:
        out.write("classification\tcount\n")
        for classification, count in sorted(counts.items(), key=lambda kv: -kv[1]):
            out.write(f"{classification}\t{count}\n")

    family_bin: dict[str, str] = {}
    with open(args.frequency_table, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            family_bin[row["family"]] = row["bin"]
    strain_rows = per_strain_summary(args.presence_matrix, family_bin)
    with open(out_dir / "per_strain_summary.tsv", "w", newline="") as out:
        fieldnames = ["Short", "n_genes", "core", "soft_core", "shell", "cloud", "singleton"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(strain_rows, key=lambda r: r["n_genes"]):
            writer.writerow(row)

    print(f"pangenome_report_tables: wrote islands_with_domains.tsv "
          f"({len(annotated)} islands), island_size_distribution.tsv, "
          f"classification_counts.tsv, per_strain_summary.tsv "
          f"({len(strain_rows)} strains) to {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
