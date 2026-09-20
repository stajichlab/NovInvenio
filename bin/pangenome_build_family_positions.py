#!/usr/bin/env python3
"""Join gene_positions.tsv (protein_id -> genomic position, per strain)
with a tier-1 cluster TSV (protein_id -> family_id) into per-strain,
per-family gene-order RANK positions -- the input lib/pangenome_synteny.py's
`linkage_fraction` and pangenome_pair_classification.py actually consume.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/build_family_positions.py).

Rank, not raw bp coordinate: `linkage_fraction`'s `k` parameter is a gene-
COUNT window ("within k genes of each other") -- a protein's rank is its
index (0-based) among all proteins on the SAME contig, sorted by start
position. A family with multiple copies in one strain gets a LIST of
(contig, rank) positions per strain, not one.

Usage:
  pangenome_build_family_positions.py --gene_positions gene_positions.tsv \\
      --cluster_tsv tier1_cluster.tsv --output family_positions.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import (  # noqa: E402
    open_maybe_compressed,
    open_maybe_compressed_write,
)
from pangenome_matrix import read_cluster_tsv  # noqa: E402


def build_family_positions(
    gene_position_rows: list[tuple[str, str, str, int, int]],
    member_to_rep: dict[str, str],
    rescue_position_rows: list[tuple[str, str, str, int]] | None = None,
    id_sep: str = "|",
) -> dict[str, dict[str, list[tuple[str, int]]]]:
    """`gene_position_rows`: [(Short, protein_id, contig, start, end), ...],
    resolved to a family via `member_to_rep` (GFF3-annotated genes).
    `rescue_position_rows`: [(Short, family, contig, start), ...] -- the
    family is already known directly (a genome-level tblastn rescue hit has
    no protein_id/gene model to look up; see
    pangenome_extract_rescue_positions.py), so these skip the
    `member_to_rep` join entirely. Returns
    {Short: {family_id: [(contig, rank), ...]}}.

    A protein_id with no cluster membership (e.g. a rare gene filtered out
    before clustering) is silently skipped -- not every protein needs to be
    in a family for this to work, only ones that are.

    Merging rescue positions into the SAME per-strain ordered list (not
    appending them separately after rank assignment) is the point -- a
    rescued family's rank must be interleaved with real annotated genes by
    genomic coordinate, not tacked on at the end, or `linkage_fraction`'s
    "within k genes" adjacency check would compare ranks that don't reflect
    real physical distance."""
    # First pass: per-strain ordered (contig, start, protein_id_or_None,
    # family_or_None) list, sorted by contig then start, so a running index
    # gives contig-relative rank (a jump in the running index at a contig
    # boundary is harmless -- linkage_fraction only ever compares positions
    # when contig_a == contig_b). Exactly one of protein_id/family is set
    # per entry -- protein_id defers its family lookup to the second pass
    # (so ranks fall out inclusive of skipped/unclustered proteins in their
    # true genomic order), family is already resolved for a rescue entry.
    per_strain_ordered: dict[str, list[tuple[str, int, str | None, str | None]]] = {}
    for short, protein_id, contig, start, _end in gene_position_rows:
        per_strain_ordered.setdefault(short, []).append((contig, start, protein_id, None))
    for short, family, contig, start in (rescue_position_rows or []):
        per_strain_ordered.setdefault(short, []).append((contig, start, None, family))

    result: dict[str, dict[str, list[tuple[str, int]]]] = {}
    for short, ordered in per_strain_ordered.items():
        ordered.sort(key=lambda row: (row[0], row[1]))
        family_positions: dict[str, list[tuple[str, int]]] = {}
        for rank, (contig, _start, protein_id, family) in enumerate(ordered):
            if family is None:
                member_id = f"{short}{id_sep}{protein_id}"
                family = member_to_rep.get(member_id)
            if family is None:
                continue
            family_positions.setdefault(family, []).append((contig, rank))
        result[short] = family_positions
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene_positions", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument(
        "--rescue_positions",
        help="optional Short/family/contig/start TSV from "
        "pangenome_extract_rescue_positions.py -- genomic positions for "
        "GENOME_ONLY (rescue-pass) presence calls, which have no "
        "GFF3-annotated protein_id to resolve via --gene_positions/--cluster_tsv",
    )
    ap.add_argument("--id_sep", default="|")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    member_to_rep = read_cluster_tsv(args.cluster_tsv)

    rows = []
    with open_maybe_compressed(args.gene_positions) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            short, protein_id, contig, start, end = parts
            rows.append((short, protein_id, contig, int(start), int(end)))

    rescue_rows = []
    if args.rescue_positions:
        with open_maybe_compressed(args.rescue_positions) as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                short, family, contig, start = parts
                rescue_rows.append((short, family, contig, int(start)))

    positions = build_family_positions(
        rows, member_to_rep, rescue_position_rows=rescue_rows, id_sep=args.id_sep,
    )

    n_rows = 0
    with open_maybe_compressed_write(args.output) as out:
        out.write("Short\tfamily\tcontig\trank\n")
        for short, family_positions in positions.items():
            for family, entries in family_positions.items():
                for contig, rank in entries:
                    out.write(f"{short}\t{family}\t{contig}\t{rank}\n")
                    n_rows += 1

    print(f"build_family_positions: {n_rows} (strain, family, copy) positions written", file=sys.stderr)


if __name__ == "__main__":
    main()
