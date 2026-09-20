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
