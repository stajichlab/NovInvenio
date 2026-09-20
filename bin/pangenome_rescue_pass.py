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
    lines: list[str],
    min_pident: float = 90.0,
    min_qcov: float = 80.0,
    stats: dict | None = None,
) -> set[tuple[str, str]]:
    """Parse tblastn outfmt6 lines (with a trailing qcovs column) into the
    set of (family_rep, strain) pairs with a qualifying genomic hit.
    Subject IDs are expected as '<Short>|<contig>' (via renamed genome FASTA
    headers before makeblastdb, so the strain is recoverable from the hit).

    If `stats` is given, it is updated in place with two funnel counters
    (issue #126): `rows_parsed` -- every non-blank, non-comment data row
    seen, regardless of whether it is well-formed or passes threshold (this
    is the "did we actually read real data" count) -- and
    `rows_passed_threshold` -- rows that met `min_pident`/`min_qcov`, before
    deduping into the returned set."""
    if stats is None:
        stats = {}
    stats.setdefault("rows_parsed", 0)
    stats.setdefault("rows_passed_threshold", 0)
    hits: set[tuple[str, str]] = set()
    for line in lines:
        line = line.rstrip("\n")
        # Skip blank lines and comment lines
        if not line or line.startswith("#"):
            continue
        stats["rows_parsed"] += 1
        parts = line.split("\t")
        if len(parts) < 13:  # Minimum columns for outfmt6 with qcovs
            continue
        try:
            family, subject, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[-1])
        except (ValueError, IndexError):
            continue
        if pident >= min_pident and qcovs >= min_qcov:
            stats["rows_passed_threshold"] += 1
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
    ap.add_argument(
        "--allow_zero_rescue", action="store_true", default=False,
        help="Do not fail when zero hit rows are parsed across all "
        "--tblastn_tsv inputs. Use only when a study genuinely expects no "
        "rescue hits -- otherwise zero rows parsed almost always means the "
        "inputs were unreadable, mis-staged, or otherwise empty by mistake "
        "(issue #126). Has no effect when --tblastn_tsv is omitted entirely, "
        "which is the separate, always-legitimate 'zero ABSENT calls "
        "upstream' case.",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)

    hits: set[tuple[str, str]] = set()
    total_rows_parsed = 0
    total_rows_passed = 0
    for tblastn_path in args.tblastn_tsv:
        file_stats: dict = {}
        with open_maybe_compressed(tblastn_path) as fh:
            hits |= parse_tblastn_hits(fh.readlines(), args.min_pident, args.min_qcov, stats=file_stats)
        total_rows_parsed += file_stats["rows_parsed"]
        total_rows_passed += file_stats["rows_passed_threshold"]
        print(
            f"  {tblastn_path}: {file_stats['rows_parsed']} rows parsed, "
            f"{file_stats['rows_passed_threshold']} passed identity/coverage thresholds",
            file=sys.stderr,
        )

    files_read = len(args.tblastn_tsv)
    print(
        f"Rescue funnel: {files_read} input file(s), {total_rows_parsed} hit rows parsed, "
        f"{total_rows_passed} rows passed identity/coverage thresholds, "
        f"{len(hits)} unique (family, strain) hits",
        file=sys.stderr,
    )

    # Issue #126: files WERE provided but produced not a single parsed row --
    # the exact shape of the symlink-staging incident (530 files present,
    # all read empty). Zero *files* (files_read == 0) is a separate,
    # always-legitimate case (zero ABSENT calls upstream) and is not gated
    # here -- see the --tblastn_tsv help text and the existing
    # test_main_succeeds_with_zero_tblastn_tsv_args regression test.
    if files_read > 0 and total_rows_parsed == 0:
        message = (
            f"zero hit rows parsed across {files_read} --tblastn_tsv input "
            "file(s). This almost never means 'no rescue hits exist' -- it "
            "means the inputs were unreadable or mis-staged (e.g. Nextflow "
            "staging files as symlinks that a decompressor refused), a "
            "decompression step failed or silently produced nothing, or an "
            "upstream TBLASTN step produced empty output. Pass "
            "--allow_zero_rescue if this run genuinely expects zero rescue "
            "hits."
        )
        if args.allow_zero_rescue:
            print(f"WARNING: {message}", file=sys.stderr)
        else:
            print(f"ERROR: {message}", file=sys.stderr)
            sys.exit(1)

    applied, skipped = apply_rescue(matrix, hits)
    print(f"Rescue: {applied} ABSENT→GENOME_ONLY, {skipped} skipped (unrecognized strain/family)", file=sys.stderr)

    if hits and skipped == len(hits):
        print(f"ERROR: All {len(hits)} parsed tblastn hits were skipped (likely wrong genome-DB naming)", file=sys.stderr)
        sys.exit(1)

    if total_rows_parsed > 0 and applied == 0:
        print(
            "WARNING: hit rows were parsed but zero cells were applied to "
            "the matrix -- the identity/coverage thresholds "
            "(--min_pident/--min_qcov) or strain/family naming may have "
            "rejected every hit. This can be a legitimate outcome, but "
            "verify it is expected.",
            file=sys.stderr,
        )

    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
