"""Build the island synteny view's JSON payload.

For one accessory island: rows = strains, columns = member families in locus
order, filled cell = family present in that strain. Deletion breakpoints show
up as vertical edges.

No I/O at import time; every function here is pure and takes already-parsed
rows, so the whole payload is unit-testable without touching a pipeline run.
"""
from __future__ import annotations


def _as_int(value, default: int = 0) -> int:
    """Parse a TSV cell to int, tolerating '' and the '-' missing sentinel."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def select_islands(rows: list[dict], min_strains: int = 2,
                   top_n: int = 50) -> list[dict]:
    """Islands worth drawing, largest first.

    Two filters, in this order:

    1. A locus must be resolved (`locus_id != '-'`), since the view is
       locus-ordered and an unplaced island has no column order.
    2. The island must be carried by at least `min_strains` strains. On the
       real genus_vs_ureesii run 17,218 of 27,836 located islands (62%) are
       single-strain, and all 20 of the largest were -- an unfiltered size
       ranking is therefore entirely strain-private content, and a
       single-strain island's grid is one filled row with no breakpoint in
       it. This is the same filter PR #110 added to report.md's top-islands
       table, for the same reason.

    Truncation to `top_n` happens AFTER filtering, so the filter never costs
    slots in the output.
    """
    located = [r for r in rows if r.get("locus_id", "-") not in ("-", "", None)]
    kept = [r for r in located if _as_int(r.get("n_strains")) >= min_strains]
    kept.sort(key=lambda r: -_as_int(r.get("island_size")))
    return kept[:top_n]


def collapse_haplotypes(presence_rows: dict[str, list[bool]]) -> list[dict]:
    """Collapse identical strain rows into distinct presence patterns.

    Returns one entry per distinct pattern: {"pattern", "count", "strains"},
    sorted lexicographically on the 0/1 pattern string. That sort is
    equivalent to a Hamming-distance leaf ordering for this data and is what
    makes deletion breakpoints align into visible vertical edges.

    Row collapsing is the load-bearing decision of this view: without it, a
    530-strain study draws 530 rows for a 20-56 column island and nothing is
    legible. With it, the row count is the number of haplotypes, typically
    tens.
    """
    by_pattern: dict[str, list[str]] = {}
    for strain, vector in presence_rows.items():
        pattern = "".join("1" if v else "0" for v in vector)
        by_pattern.setdefault(pattern, []).append(strain)
    return [
        {"pattern": p, "count": len(s), "strains": sorted(s)}
        for p, s in sorted(by_pattern.items())
    ]
