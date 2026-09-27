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


# ---- 2. exemplar -----------------------------------------------------------
@dataclass(frozen=True)
class Placement:
    """Where one strain carries a locus's root variant."""
    strain: str
    contig: str
    lo: int
    hi: int
    left_avail: int
    right_avail: int
    contig_genes: int


def _min_cover_window(copies: list[tuple[int, str]], need: int) -> tuple[int, int] | None:
    """Smallest rank window over sorted (rank, family) copies that holds
    every one of `need` distinct families. Ties: lowest start."""
    best = None
    count: dict[str, int] = collections.Counter()
    have = 0
    left = 0
    for right, (r_rank, r_fam) in enumerate(copies):
        count[r_fam] += 1
        if count[r_fam] == 1:
            have += 1
        while have == need:
            lo, hi = copies[left][0], r_rank
            if best is None or hi - lo < best[1] - best[0]:
                best = (lo, hi)
            l_fam = copies[left][1]
            count[l_fam] -= 1
            if count[l_fam] == 0:
                have -= 1
            left += 1
    return best


def carrier_placements(members: list[str], positions: Positions, strains: list[str],
                       spans: ContigSpans, k: int = DEFAULT_K) -> list[Placement]:
    """Strains that carry the variant `members`, one Placement each.

    A strain carries it when one contig holds a copy of every member family
    within a span of at most len(set(members)) - 1 + k ranks (Ruling R3).
    Per strain the tightest window wins; ties go to the contig name.
    """
    fams = sorted(set(members))
    need = len(fams)
    out = []
    for s in strains:
        by_contig: dict[str, list[tuple[int, str]]] = collections.defaultdict(list)
        missing = False
        for f in fams:
            copies = positions.get((s, f))
            if not copies:
                missing = True
                break
            for contig, rank in copies:
                by_contig[contig].append((rank, f))
        if missing:
            continue
        best = None
        for contig in sorted(by_contig):
            copies = sorted(by_contig[contig])
            win = _min_cover_window(copies, need)
            if win is None or win[1] - win[0] > need - 1 + k:
                continue
            if best is None or win[1] - win[0] < best[2] - best[1]:
                best = (contig, win[0], win[1])
        if best is None:
            continue
        contig, lo, hi = best
        cmin, cmax = spans.get((s, contig), (lo, hi))
        out.append(Placement(s, contig, lo, hi, lo - cmin, cmax - hi, cmax - cmin + 1))
    return out


def quality_key(strain: str, n50: dict[str, int], contig_genes: int) -> tuple:
    """Sort key, best first. Strains with an N50 rank before strains without
    one; with an N50: N50 descending, then name (spec section 2). Without an
    N50 (assembly_quality_vs_content.tsv covers only the ingroup): gene count
    of the locus contig descending, then name (Ruling R4)."""
    if strain in n50:
        return (0, -n50[strain], strain)
    return (1, -contig_genes, strain)


def choose_exemplar(placements: list[Placement], n50: dict[str, int],
                    flank: int = DEFAULT_FLANK,
                    flank_min: int = DEFAULT_FLANK_MIN) -> tuple[Placement, str] | None:
    """(exemplar, tier). Tier "full": >= `flank` genes on both sides.
    "short_flanks": >= `flank_min` on both sides. "contig_end": the strain(s)
    with the most flank genes in total. None when no strain carries it."""
    if not placements:
        return None

    def best(pool):
        return min(pool, key=lambda p: quality_key(p.strain, n50, p.contig_genes))

    tier = [p for p in placements if p.left_avail >= flank and p.right_avail >= flank]
    if tier:
        return best(tier), "full"
    tier = [p for p in placements if p.left_avail >= flank_min and p.right_avail >= flank_min]
    if tier:
        return best(tier), "short_flanks"
    top = max(p.left_avail + p.right_avail for p in placements)
    return best([p for p in placements if p.left_avail + p.right_avail == top]), "contig_end"


# ---- 3. columns ------------------------------------------------------------
@dataclass(frozen=True)
class Columns:
    left: tuple
    locus: tuple
    right: tuple
    left_avail: int
    right_avail: int

    @property
    def families(self) -> list[str]:
        return list(self.left) + list(self.locus) + list(self.right)


def locus_columns(order: list[tuple[int, str]], members, flank: int,
                  block: tuple[int, int] | None = None) -> Columns | None:
    """Columns from the exemplar's contig order [(rank, family), ...].

    `block` (lo, hi ranks) is the exemplar's Placement. Without it the block
    is every copy of a member family on the contig, which is the feasibility
    script's rule (used by the regression check). The locus block is every
    gene between the block ends, in rank order; the flanks are up to `flank`
    genes on each side. None when no member sits on the contig.
    """
    seq = sorted(order)
    mem = set(members)
    if block is None:
        idx = [i for i, (_, f) in enumerate(seq) if f in mem]
    else:
        idx = [i for i, (r, _) in enumerate(seq) if block[0] <= r <= block[1]]
    if not idx:
        return None
    lo, hi = min(idx), max(idx)
    return Columns(
        left=tuple(f for _, f in seq[max(0, lo - flank):lo]),
        locus=tuple(f for _, f in seq[lo:hi + 1]),
        right=tuple(f for _, f in seq[hi + 1:hi + 1 + flank]),
        left_avail=lo,
        right_avail=len(seq) - 1 - hi,
    )
