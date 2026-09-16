#!/usr/bin/env python3
"""Accessory-island construction + statistical-significance gate, generalized
from studies/fungi/Afumigatus_pangenome/bin/synteny_windows.py +
find_accessory_islands.py (NovInvenio_Investigations repo) into a permanent
nf_NovInvenio pipeline step -- see
notes/superpowers/specs/2026-09-16-pangenome-island-pfam-enrichment-design.md.

Consumes family_positions.tsv's already-computed per-strain `rank` column
directly (no GFF3 re-parsing needed -- accessory_islands() only needs
(gene_id, contig) pairs in traversal order, and rank already IS that order).
"Significant" island membership is gated entirely by pair_classification.tsv's
`classification` column (starship_explained/unexplained_physical/
ambiguous_linkage) -- that column already reflects the upstream FDR-filtered,
exact-permutation-tested co-occurrence result; this script does not re-derive
or recompute any significance test.

Usage:
  pangenome_build_islands.py --family_positions family_positions.tsv \\
      --frequency_table frequency_table.tsv \\
      --pair_classification pair_classification.tsv \\
      --cluster_tsv tier1_cluster.tsv \\
      [--marker_tblout captain=captain_vs_study.tblout] \\
      [--marker_tblout sm_backbone=SM_backbone_vs_study.tblout] \\
      --output significant_islands.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv, iter_tblout_family_hits  # noqa: E402


PHYSICAL_CLASSIFICATIONS = frozenset(
    {"starship_explained", "unexplained_physical", "ambiguous_linkage"}
)


def accessory_islands(
    gene_order: list[tuple], is_core: dict[str, bool]
) -> list[list[tuple]]:
    """Maximal runs of consecutive non-core genes, never crossing a contig
    boundary. A gene_id absent from `is_core` defaults to core (conservative:
    an unrecognized gene more likely reflects a GFF3-ID-vs-presence-matrix-ID
    mismatch than a genuine accessory gene) -- warns if >10% of genes seen are
    missing from `is_core`, since that's the real failure signature of such a
    mismatch (every island silently collapsing to nothing, no exception
    raised)."""
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    missing = 0
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
        if gene_id not in is_core:
            missing += 1
        core = is_core.get(gene_id, True)
        if contig != prev_contig and current:
            islands.append(current)
            current = []
        if not core:
            current.append(gene)
        elif current:
            islands.append(current)
            current = []
        prev_contig = contig
    if current:
        islands.append(current)
    if gene_order and missing / len(gene_order) > 0.1:
        warnings.warn(
            f"accessory_islands: {missing}/{len(gene_order)} genes "
            f"({missing / len(gene_order):.0%}) were absent from is_core and "
            "defaulted to core -- check for a GFF3 ID vs presence-matrix "
            "protein ID format mismatch before trusting these islands.",
            stacklevel=2,
        )
    return islands


def load_is_core(frequency_table_path: str) -> dict[str, bool]:
    """{family: True if core/soft_core, False if shell/cloud/singleton}."""
    is_core: dict[str, bool] = {}
    with open(frequency_table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            is_core[row["family"]] = row["bin"] in ("core", "soft_core")
    return is_core


def load_strain_gene_orders(family_positions_path: str) -> dict[str, list[tuple]]:
    """{Short: [(family, contig, rank, rank), ...]}, sorted by rank."""
    by_strain: dict[str, list[tuple]] = {}
    with open(family_positions_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            short, family, contig, rank = row["Short"], row["family"], row["contig"], int(row["rank"])
            by_strain.setdefault(short, []).append((family, contig, rank, rank))
    for rows in by_strain.values():
        rows.sort(key=lambda r: r[2])
    return by_strain


def load_significant_physical_pairs(pair_classification_path: str) -> dict[frozenset, str]:
    """{frozenset({family_a, family_b}): classification} for every physically-linked pair."""
    pairs: dict[frozenset, str] = {}
    with open(pair_classification_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["classification"] in PHYSICAL_CLASSIFICATIONS:
                pairs[frozenset({row["family_a"], row["family_b"]})] = row["classification"]
    return pairs


def load_hit_families(tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|") -> set[str]:
    """Families with >=1 hmmsearch/hmmscan hit anywhere in the cohort --
    generic across any tblout (captain, SM-backbone, or any future named
    marker search). Parsing itself lives in
    lib/pangenome_matrix.iter_tblout_family_hits (shared with
    pangenome_pair_classification.py's load_captain_families)."""
    families: set[str] = set()
    for _short, family in iter_tblout_family_hits(tblout_path, member_to_rep, id_sep=id_sep):
        families.add(family)
    return families


def build_pair_index(significant_pairs: dict[frozenset, str]) -> dict[str, list[tuple[frozenset, str]]]:
    """{family: [(pair, classification), ...]} -- avoids an O(islands x
    total_pairs) full scan; indexed lookup is O(islands x island_size x
    avg_pairs_per_family), the difference between finishing and not (a naive
    full-scan version was killed after 12+ minutes with no output on the
    Afumigatus real-data run)."""
    index: dict[str, list[tuple[frozenset, str]]] = {}
    for pair, classification in significant_pairs.items():
        for family in pair:
            index.setdefault(family, []).append((pair, classification))
    return index


def find_significant_islands(
    strain_gene_orders: dict[str, list[tuple]],
    is_core: dict[str, bool],
    significant_pairs: dict[frozenset, str],
) -> list[dict]:
    """Per strain: build accessory islands, keep only islands containing >=1
    significant physically-linked pair (both members present in that SAME
    island in that strain)."""
    pairs_by_family = build_pair_index(significant_pairs)
    results: list[dict] = []
    for short, gene_order in strain_gene_orders.items():
        islands = accessory_islands(gene_order, is_core)
        for island in islands:
            if len(island) < 2:
                continue
            members = [g[0] for g in island]
            member_set = set(members)
            candidates: dict[frozenset, str] = {}
            for family in member_set:
                for pair, classification in pairs_by_family.get(family, ()):
                    if pair <= member_set:
                        candidates[pair] = classification
            if not candidates:
                continue
            results.append({
                "strain": short,
                "contig": island[0][1],
                "size": len(island),
                "members": members,
                "supporting_pairs": list(candidates.items()),
            })
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--frequency_table", required=True)
    ap.add_argument("--pair_classification", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--min_island_size", type=int, default=2)
    ap.add_argument(
        "--marker_tblout", action="append", default=[], metavar="NAME=PATH",
        help="Repeatable. e.g. --marker_tblout captain=captain.tblout "
             "--marker_tblout sm_backbone=sm_backbone.tblout",
    )
    ap.add_argument("--id_sep", default="|",
                     help="Short-prefix separator in clustering-input FASTA headers "
                          "(default: '|'); must match --pangenome_id_sep used "
                          "everywhere else in this pipeline run.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    markers: dict[str, str] = {}
    for entry in args.marker_tblout:
        if "=" not in entry:
            sys.exit(f"ERROR: --marker_tblout must be NAME=PATH, got: {entry!r}")
        name, path = entry.split("=", 1)
        markers[name] = path

    is_core = load_is_core(args.frequency_table)
    strain_gene_orders = load_strain_gene_orders(args.family_positions)
    significant_pairs = load_significant_physical_pairs(args.pair_classification)
    print(f"pangenome_build_islands: {len(significant_pairs)} significant physically-linked "
          f"pairs, {len(strain_gene_orders)} strains", file=sys.stderr)

    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    islands = [isl for isl in islands if isl["size"] >= args.min_island_size]
    print(f"pangenome_build_islands: {len(islands)} significant islands found across all strains",
          file=sys.stderr)

    marker_families: dict[str, set[str]] = {name: set() for name in markers}
    if markers:
        member_to_rep = read_cluster_tsv(args.cluster_tsv)
        for name, path in markers.items():
            marker_families[name] = load_hit_families(path, member_to_rep, id_sep=args.id_sep)

    marker_names = sorted(markers)
    with open(args.output, "w") as out:
        header = (
            "n_strains\texample_strain\tisland_size\tmember_families\t"
            "n_supporting_pairs\tclassifications"
        )
        for name in marker_names:
            header += f"\thas_{name}"
        out.write(header + "\n")

        by_key: dict[tuple, dict] = {}
        for island in islands:
            member_set = frozenset(island["members"])
            classifications = frozenset(c for _, c in island["supporting_pairs"])
            key = (member_set, classifications)
            entry = by_key.setdefault(key, {
                "strains": [], "size": island["size"], "members": island["members"],
                "classifications": classifications,
                "n_supporting_pairs": len(island["supporting_pairs"]),
            })
            entry["strains"].append(island["strain"])

        for (member_set, classifications), entry in sorted(
            by_key.items(), key=lambda kv: -len(kv[1]["strains"])
        ):
            row = (
                f"{len(entry['strains'])}\t{entry['strains'][0]}\t{entry['size']}\t"
                f"{','.join(entry['members'])}\t{entry['n_supporting_pairs']}\t"
                f"{','.join(sorted(classifications))}"
            )
            for name in marker_names:
                has_marker = bool(member_set & marker_families[name])
                row += f"\t{'Y' if has_marker else 'N'}"
            out.write(row + "\n")

    print(f"pangenome_build_islands: {len(by_key)} distinct significant islands "
          f"(deduped across strains) written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
