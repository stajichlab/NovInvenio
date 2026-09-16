#!/usr/bin/env python3
"""Genome-level rescue pass: fold per-strain tblastn hits (family
representative vs. that strain's own genome) into the presence matrix,
upgrading ABSENT calls to GENOME_ONLY wherever a qualifying hit exists --
catches coverage-failed protein-model absences (fragmented/split gene
models, draft-assembly annotation gaps) that protein-level clustering alone
would miss.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/rescue_pass.py). This script itself is unchanged by that study's
2026-09-15 per-strain-database redesign (see
pangenome_extract_absent_family_queries.py's docstring) -- it already
accepted repeated `--tblastn_tsv` regardless of how the search was
partitioned; only the process that PRODUCES those files (searching each
strain's absent-family queries against a per-strain genome database,
never a shared combined database) needed to change.

Usage:
  pangenome_rescue_pass.py --matrix presence_matrix.tsv \\
      --tblastn_tsv strain1.tblastn.tsv.zst --tblastn_tsv strain2.tblastn.tsv.zst \\
      --output presence_matrix.rescued.tsv

`--tblastn_tsv` accepts a plain, `.gz`, or `.zst` file, and may be repeated
to read several per-strain outputs without concatenating them first.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import ABSENT, GENOME_ONLY, PresenceMatrix  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


def parse_tblastn_hits(
    lines: list[str], min_pident: float = 90.0, min_qcov: float = 80.0
) -> set[tuple[str, str]]:
    """Parse tblastn outfmt6 lines (with a trailing qcovs column) into the
    set of (family_rep, strain) pairs with a qualifying genomic hit.
    Subject IDs are expected as '<Short>|<contig>' (via renamed genome FASTA
    headers before makeblastdb, so the strain is recoverable from the hit)."""
    hits: set[tuple[str, str]] = set()
    for line in lines:
        line = line.rstrip("\n")
        # Skip blank lines and comment lines
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 13:  # Minimum columns for outfmt6 with qcovs
            continue
        try:
            family, subject, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[-1])
        except (ValueError, IndexError):
            continue
        if pident >= min_pident and qcovs >= min_qcov:
            strain = subject.split("|", 1)[0]
            hits.add((family, strain))
    return hits


def apply_rescue(matrix: PresenceMatrix, rescue_hits: set[tuple[str, str]]) -> tuple[int, int]:
    """Upgrade ABSENT calls to GENOME_ONLY wherever a qualifying rescue
    hit exists. Never downgrades an existing PRESENT call.

    Returns: (num_applied, num_skipped) where skipped = hits with
    unrecognized strain or family."""
    applied = 0
    skipped = 0
    for family, strain in rescue_hits:
        if strain not in matrix.strains or family not in matrix.families:
            skipped += 1
            continue
        if matrix.call(family, strain) == ABSENT:
            matrix.set_call(family, strain, GENOME_ONLY)
            applied += 1
    return applied, skipped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument(
        "--tblastn_tsv", action="append", default=[],
        help="plain, .gz, or .zst tblastn outfmt6+qcovs file; repeat to read "
        "several per-strain outputs without concatenating them first. "
        "May be omitted entirely (e.g. zero ABSENT calls upstream, so the "
        "workflow feeds an empty collected channel) -- in that case no "
        "rescue is applied and the matrix passes through unchanged.",
    )
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--min_qcov", type=float, default=80.0)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    hits: set[tuple[str, str]] = set()
    for tblastn_path in args.tblastn_tsv:
        with open_maybe_compressed(tblastn_path) as fh:
            hits |= parse_tblastn_hits(fh.readlines(), args.min_pident, args.min_qcov)

    applied, skipped = apply_rescue(matrix, hits)
    print(f"Rescue: {applied} ABSENT→GENOME_ONLY, {skipped} skipped (unrecognized strain/family)", file=sys.stderr)

    if hits and skipped == len(hits):
        print(f"ERROR: All {len(hits)} parsed tblastn hits were skipped (likely wrong genome-DB naming)", file=sys.stderr)
        sys.exit(1)

    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
