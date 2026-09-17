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
from pangenome_build_presence_matrix import split_member_id  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv  # noqa: E402


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


def add_island_locus(
    islands_rows: list[dict],
    member_to_rep: dict[str, str],
    gene_positions: dict[tuple[str, str], dict],
    id_sep: str = "|",
) -> list[dict]:
    """Adds locus_id/locus_contig/locus_start/locus_end/
    n_members_with_coordinates/n_contigs_in_locus to each island row, by
    resolving member_families (rep-protein IDs, comma-joined -- always,
    regardless of `id_sep`, matching bin/pangenome_build_islands.py's
    hardcoded ','.join(entry['members'])) down to the island's own
    example_strain's actual proteins via member_to_rep (rep -> member
    inverted, restricted to that strain), then spanning gene_positions.
    `id_sep` is used only to split each individual resolved member's
    `Short<id_sep>protein_id` prefix (see
    bin/pangenome_build_presence_matrix.py's split_member_id) -- member_to_rep
    values are always Short-prefixed like that in real pipeline data, while
    gene_positions is keyed on the bare protein_id. Members with no
    resolvable coordinate (e.g. a rescue-pass genome_only call with no
    annotated protein_id) are excluded from the span and counted, not
    treated as an error. A same-family member belonging to a strain other
    than this island's own example_strain is not a candidate for this
    island's locus. n_contigs_in_locus > 1 is reported, not silently
    collapsed (probable paralog-copy pull-in)."""
    # Invert member_to_rep (member -> rep) into rep -> [members on this strain].
    rep_to_members: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        rep_to_members.setdefault(rep, []).append(member)

    out = []
    for row in islands_rows:
        strain = row["example_strain"]
        starts, ends, contigs = [], [], set()
        n_resolved = 0
        for family in row["member_families"].split(","):
            resolved_pos = None
            for candidate in rep_to_members.get(family, [family]):
                cand_strain, cand_protein = split_member_id(candidate, id_sep)
                if cand_strain != strain:
                    continue
                if (strain, cand_protein) in gene_positions:
                    resolved_pos = gene_positions[(strain, cand_protein)]
                    break
            if resolved_pos is None:
                continue
            starts.append(resolved_pos["start"])
            ends.append(resolved_pos["end"])
            contigs.add(resolved_pos["contig"])
            n_resolved += 1

        new_row = dict(row)
        if n_resolved == 0:
            new_row.update({
                "locus_id": "-", "locus_contig": "-", "locus_start": "-",
                "locus_end": "-", "n_members_with_coordinates": 0,
                "n_contigs_in_locus": 0,
            })
        else:
            contig = sorted(contigs)[0]
            new_row.update({
                "locus_id": f"{strain}:{contig}:{min(starts)}-{max(ends)}",
                "locus_contig": contig, "locus_start": min(starts),
                "locus_end": max(ends), "n_members_with_coordinates": n_resolved,
                "n_contigs_in_locus": len(contigs),
            })
        out.append(new_row)
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


def marker_summary(islands_rows: list[dict], fieldnames: list[str]) -> list[dict]:
    """One row per marker name found in the `has_<marker>` columns of
    significant_islands.tsv: {marker_name, n_islands_with_marker,
    n_islands_total, pct_islands_with_marker}. `fieldnames` is
    significant_islands.tsv's own header (not `islands_rows[0].keys()`, which
    would be empty/unavailable when there are zero island rows) -- passed
    in explicitly so a header-only marker_summary.tsv can still be written
    when no marker search was run for this study (zero `has_<marker>`
    columns present) without erroring on an empty `islands_rows`."""
    marker_names = sorted(fn[len("has_"):] for fn in fieldnames if fn.startswith("has_"))
    n_total = len(islands_rows)
    rows = []
    for name in marker_names:
        n_with = sum(1 for r in islands_rows if r.get(f"has_{name}") == "Y")
        pct = round(100 * n_with / n_total, 1) if n_total else 0.0
        rows.append({
            "marker_name": name,
            "n_islands_with_marker": n_with,
            "n_islands_total": n_total,
            "pct_islands_with_marker": pct,
        })
    return rows


def per_strain_summary(presence_matrix_path: str, family_bin: dict[str, str]) -> list[dict]:
    """One row per strain: total families present (`n_families` -- paralogs
    collapse to one family in the presence matrix, so this is a family
    count, not a raw gene count), plus counts broken down by frequency bin
    (core/soft_core/shell/cloud/singleton). `genome_only` counts as present
    (matches lib/pangenome_matrix.PresenceMatrix.is_present's semantics) --
    rescued genome-only calls are real presence evidence, not a weaker
    state."""
    with open(presence_matrix_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        strains = header[1:]
        totals = {
            s: {"Short": s, "n_families": 0, "core": 0, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0}
            for s in strains
        }
        for row in reader:
            family = row[0]
            b = family_bin.get(family)
            for strain, call in zip(strains, row[1:]):
                if call != "absent":
                    totals[strain]["n_families"] += 1
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
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.significant_islands, newline="") as fh:
        significant_islands_reader = csv.DictReader(fh, delimiter="\t")
        islands_rows = list(significant_islands_reader)
        significant_islands_fieldnames = significant_islands_reader.fieldnames or []

    family_domains = parse_domtblout(args.domtblout, max_ievalue=args.domain_evalue)
    annotated = annotate_islands_with_domains(islands_rows, family_domains)

    member_to_rep = read_cluster_tsv(args.cluster_tsv)
    gene_positions: dict[tuple[str, str], dict] = {}
    with open(args.gene_positions, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene_positions[(row["Short"], row["protein_id"])] = {
                "contig": row["contig"], "start": int(row["start"]), "end": int(row["end"]),
            }
    annotated = add_island_locus(annotated, member_to_rep, gene_positions, id_sep=args.id_sep)

    with open(out_dir / "islands_with_domains.tsv", "w", newline="") as out:
        fieldnames = list(annotated[0].keys()) if annotated else [
            "n_strains", "example_strain", "island_size", "member_families",
            "n_supporting_pairs", "classifications", "pfam_domains",
            "locus_id", "locus_contig", "locus_start", "locus_end",
            "n_members_with_coordinates", "n_contigs_in_locus",
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

    markers = marker_summary(islands_rows, significant_islands_fieldnames)
    with open(out_dir / "marker_summary.tsv", "w", newline="") as out:
        fieldnames = ["marker_name", "n_islands_with_marker", "n_islands_total", "pct_islands_with_marker"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in markers:
            writer.writerow(row)

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
        fieldnames = ["Short", "n_families", "core", "soft_core", "shell", "cloud", "singleton"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(strain_rows, key=lambda r: r["n_families"]):
            writer.writerow(row)

    print(f"pangenome_report_tables: wrote islands_with_domains.tsv "
          f"({len(annotated)} islands), island_size_distribution.tsv, "
          f"classification_counts.tsv, marker_summary.tsv ({len(markers)} markers), "
          f"per_strain_summary.tsv ({len(strain_rows)} strains) to {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
