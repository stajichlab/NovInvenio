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
from pangenome_matrix import (  # noqa: E402
    ABSENT,
    GENOME_ONLY,
    PresenceMatrix,
    read_cluster_tsv,
    read_gene_positions_tsv,
)
from compressed_io import open_maybe_compressed  # noqa: E402


def read_fasta_lengths(path: str) -> dict[str, int]:
    """{record_id: sequence length} for a FASTA file, keyed by the first
    whitespace-delimited token after '>' -- matches the "family ID = tier-1
    cluster representative ID, verbatim" convention (see
    pangenome_extract_absent_family_queries.py's read_fasta_records), so a
    presence-matrix family name looks up its representative's length
    directly."""
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split(None, 1)[0].rstrip("\n")
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            lengths[current_id] = current_len
    return lengths


class StructuralFilter:
    """Issue #133: rescue-as-configured is dominated by hits that land on a
    predicted gene already assigned to a DIFFERENT family in that strain
    (77.6% of all rescuable cells, exhaustively measured) -- the same locus
    counted twice, not a genuine annotation dropout. Thresholds (identity/
    coverage) cannot separate this from a real hit: the class mix stays flat
    across the whole grid because these hits are 100%-identical real DNA.

    Three structural signals, checked in this order (each measured to
    remove a specific share of the false-positive rescuable cells):
      1. `overlapping_other_family` -- the hit's span overlaps an existing
         predicted gene (from `gene_positions.tsv`) already assigned to a
         different family. Removes 77.6% of cells outright.
      2. `rep_too_short` -- the query family representative is under
         `min_rep_length` aa (default 150; 83.2% of the remaining
         intergenic cells come from reps under this, median real family rep
         is 201 aa).
      3. Repeat-hotspot windows (applied in `filter_candidates`, not a
         single-hit method, since it needs every candidate on a
         strain/contig at once) -- a `hotspot_window_bp` window hit by
         `hotspot_min_families` or more distinct families is a repeat/
         mobile-element region, not that many independent real gene losses
         (47.5% of the remaining intergenic cells).
    """

    def __init__(
        self,
        gene_positions: list[tuple[str, str, str, int, int]],
        member_to_rep: dict[str, str],
        family_lengths: dict[str, int],
        id_sep: str = "|",
        min_rep_length: int = 150,
        hotspot_window_bp: int = 2000,
        hotspot_min_families: int = 3,
    ) -> None:
        self.family_lengths = family_lengths
        self.min_rep_length = min_rep_length
        self.hotspot_window_bp = hotspot_window_bp
        self.hotspot_min_families = hotspot_min_families
        # (strain, contig) -> [(start, end, family), ...] sorted by start,
        # so overlap queries can stop scanning once a gene starts past the
        # query's end -- a real interval tree only matters at a scale this
        # pipeline doesn't reach (a strain's own contig's own gene count).
        self._genes: dict[tuple[str, str], list[tuple[int, int, str]]] = {}
        for short, protein_id, contig, start, end in gene_positions:
            family = member_to_rep.get(f"{short}{id_sep}{protein_id}")
            if family is None:
                continue
            self._genes.setdefault((short, contig), []).append((start, end, family))
        for genes in self._genes.values():
            genes.sort()

    def overlapping_other_family(self, strain: str, contig: str, start: int, end: int, family: str) -> bool:
        """True if [start, end] overlaps a predicted gene on this
        strain/contig that is assigned to a DIFFERENT family than `family`."""
        for g_start, g_end, g_family in self._genes.get((strain, contig), ()):
            if g_start > end:
                break
            if g_end >= start and g_family != family:
                return True
        return False

    def rep_too_short(self, family: str) -> bool:
        """True if `family`'s representative sequence is under
        `min_rep_length` aa. A family with no recorded length (rep not
        found in the FASTA given) is never rejected on this ground alone --
        absence of data isn't evidence of a short rep."""
        length = self.family_lengths.get(family)
        return length is not None and length < self.min_rep_length

    def filter_candidates(
        self,
        candidates: list[tuple[str, str, str, int, int]],
        stats: dict,
    ) -> set[tuple[str, str]]:
        """Apply all three criteria to one file's (already identity/coverage
        -qualifying) candidate hits -- `candidates` is
        [(family, strain, contig, start, end), ...]. The hotspot criterion
        is windowed per (strain, contig) over what survives criteria 1-2,
        which is correct scope here since each tblastn file is one strain's
        hits only (see modules/pangenome/rescue.nf's TBLASTN_PER_STRAIN)."""
        pass1: list[tuple[str, str, str, int, int]] = []
        for family, strain, contig, start, end in candidates:
            if self.overlapping_other_family(strain, contig, start, end, family):
                stats["rows_rejected_overlap"] += 1
                continue
            if self.rep_too_short(family):
                stats["rows_rejected_short_rep"] += 1
                continue
            pass1.append((family, strain, contig, start, end))

        windows: dict[tuple[str, str, int], set[str]] = {}
        for family, strain, contig, start, _end in pass1:
            w = start // self.hotspot_window_bp
            windows.setdefault((strain, contig, w), set()).add(family)

        accepted: set[tuple[str, str]] = set()
        for family, strain, contig, start, _end in pass1:
            w = start // self.hotspot_window_bp
            if len(windows[(strain, contig, w)]) >= self.hotspot_min_families:
                stats["rows_rejected_hotspot"] += 1
                continue
            accepted.add((family, strain))
        return accepted


def parse_tblastn_hits(
    lines: list[str],
    min_pident: float = 90.0,
    min_qcov: float = 80.0,
    stats: dict | None = None,
    structural: "StructuralFilter | None" = None,
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
    deduping into the returned set.

    `structural`, when given a `StructuralFilter` (issue #133), additionally
    rejects a threshold-qualifying row whose best-HSP span lands on an
    existing gene of a different family, whose query family rep is
    implausibly short, or that falls in a repeat-hotspot window -- see
    `StructuralFilter`'s docstring. When `structural` is None (the default),
    behavior is unchanged from before issue #133."""
    if stats is None:
        stats = {}
    stats.setdefault("rows_parsed", 0)
    stats.setdefault("rows_passed_threshold", 0)
    if structural is not None:
        stats.setdefault("rows_rejected_overlap", 0)
        stats.setdefault("rows_rejected_short_rep", 0)
        stats.setdefault("rows_rejected_hotspot", 0)

    candidates: list[tuple[str, str, str, int, int]] = []
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
        if pident < min_pident or qcovs < min_qcov:
            continue
        stats["rows_passed_threshold"] += 1
        strain, _, contig = subject.partition("|")
        if structural is None:
            candidates.append((family, strain, contig, 0, 0))
            continue
        try:
            sstart, send = int(parts[8]), int(parts[9])
        except (ValueError, IndexError):
            continue
        start, end = (sstart, send) if sstart <= send else (send, sstart)
        candidates.append((family, strain, contig, start, end))

    if structural is None:
        return {(family, strain) for family, strain, _contig, _start, _end in candidates}
    return structural.filter_candidates(candidates, stats)


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
    ap.add_argument(
        "--gene_positions",
        help="gene_positions.tsv (plain/.gz/.zst) -- enables the issue #133 "
        "structural rescue filter. Must be given together with "
        "--cluster_tsv and --rep_fasta, or not at all.",
    )
    ap.add_argument(
        "--cluster_tsv",
        help="tier1_cluster.tsv (rep\\tmember) -- resolves a gene_positions "
        "protein_id to its family, for the structural filter's overlap "
        "check. Must be given together with --gene_positions/--rep_fasta.",
    )
    ap.add_argument(
        "--rep_fasta",
        help="tier1_rep_seq.fasta (family representative sequences) -- "
        "supplies rep length for the structural filter's min-length check. "
        "Must be given together with --gene_positions/--cluster_tsv.",
    )
    ap.add_argument("--id_sep", default="|", help="Short-prefix separator in cluster_tsv member IDs")
    ap.add_argument(
        "--rescue_min_rep_length", type=int, default=150,
        help="structural filter: reject a hit whose query family rep is under "
        "this many aa (default 150; median real family rep is ~201 aa, "
        "83.2%% of intergenic false-positive rescues come from reps under "
        "this length -- issue #133).",
    )
    ap.add_argument(
        "--rescue_hotspot_window", type=int, default=2000,
        help="structural filter: genomic window size (bp) for the repeat-"
        "hotspot check (issue #133).",
    )
    ap.add_argument(
        "--rescue_hotspot_min_families", type=int, default=3,
        help="structural filter: reject a hit if its --rescue_hotspot_window "
        "window is hit by this many or more distinct families -- a repeat/"
        "mobile-element region, not that many independent real gene losses "
        "(issue #133).",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    structural_inputs = {
        "--gene_positions": args.gene_positions,
        "--cluster_tsv": args.cluster_tsv,
        "--rep_fasta": args.rep_fasta,
    }
    given = [name for name, val in structural_inputs.items() if val]
    missing = [name for name, val in structural_inputs.items() if not val]
    if given and missing:
        ap.error(
            "the structural rescue filter (issue #133) requires "
            "--gene_positions, --cluster_tsv, and --rep_fasta together; "
            f"got {given} but not {missing}"
        )

    structural: StructuralFilter | None = None
    if given:
        gene_positions = read_gene_positions_tsv(args.gene_positions)
        member_to_rep = read_cluster_tsv(args.cluster_tsv)
        family_lengths = read_fasta_lengths(args.rep_fasta)
        structural = StructuralFilter(
            gene_positions, member_to_rep, family_lengths,
            id_sep=args.id_sep,
            min_rep_length=args.rescue_min_rep_length,
            hotspot_window_bp=args.rescue_hotspot_window,
            hotspot_min_families=args.rescue_hotspot_min_families,
        )

    matrix = PresenceMatrix.from_tsv(args.matrix)

    hits: set[tuple[str, str]] = set()
    total_rows_parsed = 0
    total_rows_passed = 0
    total_rejected_overlap = 0
    total_rejected_short_rep = 0
    total_rejected_hotspot = 0
    for tblastn_path in args.tblastn_tsv:
        file_stats: dict = {}
        with open_maybe_compressed(tblastn_path) as fh:
            hits |= parse_tblastn_hits(
                fh.readlines(), args.min_pident, args.min_qcov,
                stats=file_stats, structural=structural,
            )
        total_rows_parsed += file_stats["rows_parsed"]
        total_rows_passed += file_stats["rows_passed_threshold"]
        total_rejected_overlap += file_stats.get("rows_rejected_overlap", 0)
        total_rejected_short_rep += file_stats.get("rows_rejected_short_rep", 0)
        total_rejected_hotspot += file_stats.get("rows_rejected_hotspot", 0)
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
    if structural is not None:
        print(
            f"  structural filter (issue #133): {total_rejected_overlap} rejected: "
            f"overlaps a different family's gene, {total_rejected_short_rep} rejected: "
            f"query rep under {args.rescue_min_rep_length} aa, {total_rejected_hotspot} "
            f"rejected: repeat-hotspot window (>= {args.rescue_hotspot_min_families} "
            "families)",
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
        structural_note = (
            " or the structural filter (--gene_positions/--cluster_tsv/--rep_fasta)"
            if structural is not None else ""
        )
        print(
            "WARNING: hit rows were parsed but zero cells were applied to "
            "the matrix -- the identity/coverage thresholds "
            f"(--min_pident/--min_qcov){structural_note} or strain/family "
            "naming may have rejected every hit. This can be a legitimate "
            "outcome, but verify it is expected.",
            file=sys.stderr,
        )

    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
