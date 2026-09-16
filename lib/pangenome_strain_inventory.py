"""Reader for `strain_inventory.tsv` (bin/dereplicate_strains.py's output).

The plan's global constraint is that frequency counts use DEREPLICATED strains,
not the raw 293 -- public strain sets can carry the same isolate twice under
different names, which inflates every presence fraction. The inventory names
one representative per mash dedup group (`is_representative == 1`); this module
is how the downstream steps (frequency_bins.py, cooccurrence.py) consume it.
"""
from __future__ import annotations

from pathlib import Path


def read_representative_shorts(path: str | Path) -> list[str]:
    """Return the `Short` values with `is_representative == 1`, in file order.

    Raises:
        ValueError: if the file has no `Short`/`is_representative` columns, or
            names no representative at all (a silently-empty strain set would
            make every downstream frequency 0).
    """
    reps: list[str] = []
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")
        if "Short" not in header or "is_representative" not in header:
            raise ValueError(
                f"{path}: not a strain inventory -- expected 'Short' and "
                f"'is_representative' columns, got {header}"
            )
        short_i = header.index("Short")
        rep_i = header.index("is_representative")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(short_i, rep_i):
                continue
            if parts[rep_i].strip() == "1":
                reps.append(parts[short_i].strip())
    if not reps:
        raise ValueError(f"{path}: no rows with is_representative == 1")
    return reps
