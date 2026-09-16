#!/usr/bin/env python3
"""Genomic positions for rescue-pass (GENOME_ONLY) presence calls -- closes
a real gap: pangenome_build_family_positions.py only knows positions for
GFF3-ANNOTATED proteins (pangenome_build_gene_positions.py), so a family
present ONLY via the genome-level tblastn rescue (a real hit, but no
annotated gene model in that strain) would otherwise have no resolvable
position at all, making it invisible to pangenome_pair_classification.py's
physical-linkage test.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/extract_rescue_positions.py). A tblastn hit already carries a genomic
position directly (subject start/end on a specific contig) -- no
protein_id/gene-model crosswalk needed, unlike the GFF3 path. This
re-parses the SAME per-strain tblastn output files pangenome_rescue_pass.py
already consumed, keeping only hits for (family, strain) pairs the final
rescued matrix actually calls GENOME_ONLY. When a GENOME_ONLY (family,
strain) pair has multiple qualifying hits (e.g. real multi-copy
duplication), only the single highest-bitscore hit's position is kept.

`--processes` parses independent tblastn files in a process pool; the
original study found (2026-09-15) that the real bottleneck at 295-strain
scale was NOT file parsing but an O(N_families) linear membership scan on a
plain list, fixed here by using sets for `matrix.families`/`matrix.strains`
-- worth keeping in mind before reaching for `--processes` on a similarly
shaped bottleneck elsewhere: profile before parallelizing.

Usage:
  pangenome_extract_rescue_positions.py --matrix presence_matrix.rescued.tsv \\
      --tblastn_tsv strain1.tblastn.tsv.zst --tblastn_tsv strain2.tblastn.tsv.zst \\
      --processes 8 --output rescue_positions.tsv
  # Short<TAB>family<TAB>contig<TAB>start -- feed into
  # pangenome_build_family_positions.py's --rescue_positions
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import GENOME_ONLY, PresenceMatrix  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402


def parse_tblastn_best_hit_positions(
    lines: list[str], min_pident: float = 90.0, min_qcov: float = 80.0,
) -> dict[tuple[str, str], tuple[str, int, float]]:
    """Parse tblastn outfmt6+qcovs lines into
    {(family, strain): (contig, start, bitscore)} for the single
    highest-bitscore qualifying hit per (family, strain) pair. `start` is
    min(sstart, send) -- tblastn reports subject coordinates in hit
    orientation, so a minus-strand hit has sstart > send; normalizing to
    the smaller value matches pangenome_build_gene_positions.py's
    GFF3-derived convention (min(start), max(end)) so ranks computed later
    are ordered consistently regardless of strand."""
    best: dict[tuple[str, str], tuple[str, int, float]] = {}
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        try:
            family, subject = parts[0], parts[1]
            pident, qcovs = float(parts[2]), float(parts[-1])
            sstart, send, bitscore = int(parts[8]), int(parts[9]), float(parts[11])
        except (ValueError, IndexError):
            continue
        if pident < min_pident or qcovs < min_qcov:
            continue
        strain, _, contig = subject.partition("|")
        if not contig:
            continue
        key = (family, strain)
        start = min(sstart, send)
        existing = best.get(key)
        if existing is None or bitscore > existing[2]:
            best[key] = (contig, start, bitscore)
    return best


def _parse_one_file(
    args: tuple[str, float, float],
) -> dict[tuple[str, str], tuple[str, int, float]]:
    """Pool worker: parse a single tblastn file in isolation. A top-level
    (not nested) function so it's picklable for `multiprocessing.Pool`."""
    tblastn_path, min_pident, min_qcov = args
    with open_maybe_compressed(tblastn_path) as fh:
        return parse_tblastn_best_hit_positions(fh.readlines(), min_pident, min_qcov)


def merge_best_hits(
    per_file_hits: list[dict[tuple[str, str], tuple[str, int, float]]],
) -> dict[tuple[str, str], tuple[str, int, float]]:
    """Reduce several files' {(family, strain): (contig, start, bitscore)}
    dicts into one, keeping the highest-bitscore hit per key -- a simple
    associative/commutative max reduction, so the result doesn't depend on
    what order the per-file dicts arrive in (parallel-safe)."""
    merged: dict[tuple[str, str], tuple[str, int, float]] = {}
    for hits in per_file_hits:
        for key, (contig, start, bitscore) in hits.items():
            existing = merged.get(key)
            if existing is None or bitscore > existing[2]:
                merged[key] = (contig, start, bitscore)
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument(
        "--tblastn_tsv", action="append", default=[],
        help="plain, .gz, or .zst tblastn outfmt6+qcovs file; repeat for "
        "several per-strain files. May be omitted entirely (e.g. zero "
        "ABSENT calls upstream, so the workflow feeds an empty collected "
        "channel) -- in that case an output file with only the header is "
        "written and nothing is rescued.",
    )
    ap.add_argument("--min_pident", type=float, default=90.0)
    ap.add_argument("--min_qcov", type=float, default=80.0)
    ap.add_argument(
        "--processes", type=int, default=1,
        help="parse this many tblastn files in parallel (each file's parsing "
        "is independent -- see the module docstring). 1 (default) parses "
        "sequentially in-process, identical to the original behavior.",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)

    tasks = [(path, args.min_pident, args.min_qcov) for path in args.tblastn_tsv]
    if args.processes <= 1:
        per_file_hits = [_parse_one_file(task) for task in tasks]
    else:
        with mp.Pool(processes=args.processes) as pool:
            per_file_hits = pool.map(_parse_one_file, tasks)
    best = merge_best_hits(per_file_hits)

    # Sets, not the underlying lists: `family in matrix.families` on a large
    # LIST is an O(F) linear scan, and this membership check runs once per
    # entry in `best` (potentially hundreds of thousands to millions of
    # entries at real scale) -- the actual dominant cost at real study scale,
    # dwarfing per-file parsing. Converting to sets makes each check O(1)
    # without changing PresenceMatrix's own list-based public fields.
    strain_set = set(matrix.strains)
    family_set = set(matrix.families)

    n_written, n_not_genome_only = 0, 0
    with open(args.output, "w") as out:
        out.write("Short\tfamily\tcontig\tstart\n")
        for (family, strain), (contig, start, _bitscore) in sorted(best.items()):
            if strain not in strain_set or family not in family_set:
                continue
            if matrix.call(family, strain) != GENOME_ONLY:
                n_not_genome_only += 1
                continue
            out.write(f"{strain}\t{family}\t{contig}\t{start}\n")
            n_written += 1

    print(
        f"extract_rescue_positions: {n_written} GENOME_ONLY positions written, "
        f"{n_not_genome_only} qualifying hits skipped (matrix call was not "
        "GENOME_ONLY -- e.g. already PRESENT at the protein level)", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
