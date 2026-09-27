"""Island locus view: exemplar-anchored columns, per-strain cell states,
row classes and shared indel breakpoints.

Spec: docs/superpowers/specs/2026-09-24-island-locus-view-design.md.

Pure functions, no I/O at import time. A position is a (contig, rank) pair
from family_positions.tsv. A rank is a strain-wide running index
(bin/pangenome_build_family_positions.py sorts every strain's genes by
(contig, start) and enumerates them), so two ranks compare only on the same
contig, and a contig's first and last gene are its minimum and maximum rank.

"Ruling Rn" refers to the Rulings section of
docs/superpowers/plans/2026-09-26-island-locus-view-clinker.md.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass

# ---- cell state codes (one character per column in a row's `codes`) -------
ABSENT = "0"
IN_PLACE = "1"
ELSEWHERE = "2"
RESCUE_IN_PLACE = "3"
RESCUE_ELSEWHERE = "4"
CONTIG_BREAK = "5"
STATE_LABELS = {
    ABSENT: "absent",
    IN_PLACE: "in place",
    ELSEWHERE: "elsewhere",
    RESCUE_IN_PLACE: "rescue (in place)",
    RESCUE_ELSEWHERE: "rescue (elsewhere)",
    CONTIG_BREAK: "contig break",
}
IN_PLACE_CODES = frozenset({IN_PLACE, RESCUE_IN_PLACE})

ROW_CLASSES = ("full", "partial", "empty", "uninformative")

DEFAULT_FLANK = 5
DEFAULT_FLANK_MIN = 3
DEFAULT_K = 10
DEFAULT_EMPTY_FRAC = 0.8
DEFAULT_CONTAINMENT = 0.5
INFORMATIVE_MIN_EMPTY = 10
INFORMATIVE_MIN_FULL = 2
MAX_DETAIL_ROWS = 300

Positions = dict  # {(strain, family): [(contig, rank), ...]}
ContigSpans = dict  # {(strain, contig): (min_rank, max_rank)}


def _as_int(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def island_members(row: dict) -> list[str]:
    return [f for f in row.get("member_families", "").split(",") if f]


# ---- 1. loci ---------------------------------------------------------------
@dataclass
class Locus:
    """A group of located islands: the largest (`root`) plus every island
    that joined it (spec section 1). `variants[0]` is the root."""
    root: dict
    variants: list

    @property
    def members(self) -> list[str]:
        return self.root["members"]

    @property
    def locus_id(self) -> str:
        return self.root.get("locus_id", "-")

    @property
    def size(self) -> int:
        return len(self.root["members"])

    def carrier_proxy(self, n_strains_total: int) -> int:
        """Sum of the variants' n_strains, capped at the strain count. A
        cheap pre-filter only: one strain can carry two variants."""
        return min(n_strains_total, sum(_as_int(v.get("n_strains")) for v in self.variants))


def group_loci(island_rows: list[dict], containment: float = DEFAULT_CONTAINMENT) -> list[Locus]:
    """Group located islands into loci (spec section 1, the M4 criterion).

    Islands are visited largest first (stable on input order, as the
    feasibility script does). An island joins the locus of a strictly larger,
    earlier island when >= `containment` of its families are in that island.
    With several such islands it joins the one sharing the most families;
    a tie goes to the earlier (larger) island. It inherits that island's
    locus, so grouping is transitive to the root.
    """
    rows = []
    for r in island_rows:
        if r.get("locus_id", "-") in ("-", "", None):
            continue
        members = island_members(r)
        if members:
            rows.append(dict(r, members=members))
    by_size = sorted(rows, key=lambda r: -len(r["members"]))
    sets = [set(r["members"]) for r in by_size]
    fam_index: dict[str, list[int]] = collections.defaultdict(list)
    locus_of: list[Locus] = []
    loci: list[Locus] = []
    for i, s in enumerate(sets):
        overlap = collections.Counter(
            j for f in s for j in fam_index[f] if len(sets[j]) > len(s))
        best = None
        for j, c in overlap.items():
            if c / len(s) < containment:
                continue
            if best is None or c > overlap[best] or (c == overlap[best] and j < best):
                best = j
        if best is None:
            loc = Locus(root=by_size[i], variants=[by_size[i]])
            loci.append(loc)
        else:
            loc = locus_of[best]
            loc.variants.append(by_size[i])
        locus_of.append(loc)
        for f in s:
            fam_index[f].append(i)
    return loci


def candidate_loci(loci: list[Locus], n_strains_total: int, rank_by: str,
                   n_candidates: int, min_strains: int = 2) -> list[Locus]:
    """Loci whose states are computed. `rank_by == "size"` pre-selects by root
    size; any other rank pre-selects by carrier_proxy()."""
    kept = [loc for loc in loci if loc.carrier_proxy(n_strains_total) >= min_strains]
    if rank_by == "size":
        kept.sort(key=lambda loc: (-loc.size, -loc.carrier_proxy(n_strains_total), loc.locus_id))
    else:
        kept.sort(key=lambda loc: (-loc.carrier_proxy(n_strains_total), -loc.size, loc.locus_id))
    return kept[:n_candidates]
