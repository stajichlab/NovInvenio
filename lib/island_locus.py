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

from island_labels import column_labels

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
                   n_candidates: int, min_strains: int = 2,
                   strata: tuple = ()) -> list[Locus]:
    """Loci whose states are computed. `rank_by == "size"` pre-selects by root
    size; any other rank pre-selects by carrier_proxy().

    `strata` is a sorted tuple of carrier-fraction cut points, for example
    (0.05, 0.2, 0.5, 0.8). The loci are cut into len(strata) + 1 bins by
    carrier_proxy / n_strains_total. Each bin supplies the same number of
    candidates, largest locus first, and a bin that is short passes its share
    to the others. With no strata, the highest-carrier loci are kept (a set
    that is biased toward near-fixed loci)."""
    kept = [loc for loc in loci if loc.carrier_proxy(n_strains_total) >= min_strains]
    if rank_by == "size":
        kept.sort(key=lambda loc: (-loc.size, -loc.carrier_proxy(n_strains_total), loc.locus_id))
    else:
        kept.sort(key=lambda loc: (-loc.carrier_proxy(n_strains_total), -loc.size, loc.locus_id))
    if not strata:
        return kept[:n_candidates]
    bins: list[list[Locus]] = [[] for _ in range(len(strata) + 1)]
    for loc in kept:
        frac = loc.carrier_proxy(n_strains_total) / max(n_strains_total, 1)
        bins[sum(frac >= cut for cut in strata)].append(loc)
    for b in bins:
        b.sort(key=lambda loc: (-loc.size, loc.locus_id))
    take = [0] * len(bins)
    left = n_candidates
    while left > 0 and any(take[i] < len(b) for i, b in enumerate(bins)):
        open_bins = [i for i, b in enumerate(bins) if take[i] < len(b)]
        share = max(left // len(open_bins), 1)
        for i in open_bins:
            n = min(share, len(bins[i]) - take[i], left)
            take[i] += n
            left -= n
            if left == 0:
                break
    return [loc for i, b in enumerate(bins) for loc in b[:take[i]]]


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
                    flank_min: int = DEFAULT_FLANK_MIN,
                    outgroup: frozenset | set = frozenset(),
                    preferred: list[str] | tuple = ()) -> tuple[Placement, str] | None:
    """(exemplar, tier). Tier "full": >= `flank` genes on both sides.
    "short_flanks": >= `flank_min` on both sides. "contig_end": the strain(s)
    with the most flank genes in total. None when no strain carries it.

    `outgroup` names outgroup strains. Within a tier, an ingroup carrier is
    preferred to any outgroup carrier, whatever their assembly quality: the
    exemplar names the locus's coordinates and column labels, and the page is
    about the ingroup. An outgroup strain is chosen only when it is the sole carrier
    in the best available tier.

    `preferred` (ordered) lists reference strains: within a tier, the first listed carrier wins
    over any other, whatever its quality. The tier is still chosen first, so a reference strain
    with short flanks does not displace a strain with full flanks."""
    if not placements:
        return None
    pref = {s: i for i, s in enumerate(preferred)}

    def best(pool):
        return min(pool, key=lambda p: (pref.get(p.strain, len(pref)), p.strain in outgroup,
                                        quality_key(p.strain, n50, p.contig_genes)))

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


# ---- 4. cell states ----------------------------------------------------------
@dataclass
class StrainCells:
    """One strain's view of one locus."""
    base: list  # ABSENT / IN_PLACE / ELSEWHERE only (the feasibility rule)
    codes: list  # display codes, adds rescue and contig break
    detail: list  # per column: None or (contig, rank, neighbour_col, delta)
    in_place_copies: list  # [(col, contig, rank)] for every in-place copy


def _nearest_other_family(copies: list[tuple[int, int]], idx: int, k: int,
                          columns: list[str]) -> tuple[int, int] | None:
    """(column, signed rank delta) of the nearest copy of a DIFFERENT FAMILY
    within `k` ranks of copies[idx], in a rank-sorted list of (rank, column).
    A family that fills two columns (a tandem paralog in the exemplar) is
    never its own neighbour. Ties: the left one."""
    rank, ci = copies[idx]
    fam = columns[ci]
    left = right = None
    j = idx - 1
    while j >= 0 and rank - copies[j][0] <= k:
        if columns[copies[j][1]] != fam:
            left = (copies[j][1], copies[j][0] - rank)
            break
        j -= 1
    j = idx + 1
    while j < len(copies) and copies[j][0] - rank <= k:
        if columns[copies[j][1]] != fam:
            right = (copies[j][1], copies[j][0] - rank)
            break
        j += 1
    if left is None:
        return right
    if right is None or -left[1] <= right[1]:
        return left
    return right


def strain_cells(strain: str, columns: list[str], positions: Positions,
                 spans: ContigSpans, k: int = DEFAULT_K,
                 genome_only: frozenset = frozenset(),
                 present_unplaced: frozenset = frozenset()) -> StrainCells:
    """Cell states for one strain (spec section 4).

    A copy is in place when a copy of a different column's family sits on
    the same contig within `k` ranks (Ruling R9). `genome_only` names the
    families this strain carries only as a TBLASTN rescue hit (matrix call
    genome_only).
    `present_unplaced` names families the matrix calls present with no
    position row; they are "elsewhere" (Ruling R6).
    """
    n = len(columns)
    by_contig: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    has_copy = [False] * n
    for ci, fam in enumerate(columns):
        for contig, rank in positions.get((strain, fam), ()):
            by_contig[contig].append((rank, ci))
            has_copy[ci] = True
    best_detail: list = [None] * n
    in_place_copies = []
    for contig, copies in by_contig.items():
        copies.sort()
        for idx, (rank, ci) in enumerate(copies):
            nbr = _nearest_other_family(copies, idx, k, columns)
            if nbr is None:
                continue
            in_place_copies.append((ci, contig, rank))
            cur = best_detail[ci]
            cand = (contig, rank, nbr[0], nbr[1])
            if cur is None or (abs(cand[3]), cand[0], cand[1]) < (abs(cur[3]), cur[0], cur[1]):
                best_detail[ci] = cand
    base = []
    codes = []
    detail = []
    for ci, fam in enumerate(columns):
        if best_detail[ci] is not None:
            b = IN_PLACE
            d = best_detail[ci]
        elif has_copy[ci]:
            b = ELSEWHERE
            contig, rank = min(positions[(strain, fam)])
            d = (contig, rank, None, None)
        elif fam in present_unplaced:
            b = ELSEWHERE
            d = None
        else:
            b = ABSENT
            d = None
        base.append(b if (has_copy[ci] or b == ABSENT) else ABSENT)
        rescued = fam in genome_only and b != ABSENT
        if rescued:
            codes.append(RESCUE_IN_PLACE if b == IN_PLACE else RESCUE_ELSEWHERE)
        else:
            codes.append(b)
        detail.append(d)
    _mark_contig_breaks(codes, best_detail, strain, spans, k)
    return StrainCells(base=base, codes=codes, detail=detail, in_place_copies=in_place_copies)


def _mark_contig_breaks(codes: list, anchors: list, strain: str, spans: ContigSpans, k: int) -> None:
    """Turn ABSENT/ELSEWHERE into CONTIG_BREAK where the locus runs off the
    assembly (spec section 4, Ruling R7): the column lies outside the column
    range of the strain's in-place columns, and the nearest in-place copy is
    within `k` ranks of either end of its contig."""
    placed = [ci for ci, a in enumerate(anchors) if a is not None]
    if not placed:
        return
    first, last = placed[0], placed[-1]
    for ci, code in enumerate(codes):
        if code not in (ABSENT, ELSEWHERE) or first <= ci <= last:
            continue
        j = first if ci < first else last
        contig, rank = anchors[j][0], anchors[j][1]
        cmin, cmax = spans.get((strain, contig), (None, None))
        if cmin is None:
            continue
        if rank - cmin <= k or cmax - rank <= k:
            codes[ci] = CONTIG_BREAK


# ---- 5. row classes and breakpoints -----------------------------------------
def flank_pair(positions: Positions, strain: str, left: list[str], right: list[str],
               n_locus: int, k: int = DEFAULT_K,
               in_place: set | None = None) -> tuple[str, int, int] | None:
    """(contig, lower rank, upper rank) of the closest left-flank/right-flank
    copy pair on one contig within n_locus + 2k ranks, or None (spec section
    5, "flanks intact"). With `in_place` (a set of (contig, rank) in-place
    copies) both copies must be in place; without it any copies count, which
    is the feasibility script's rule."""
    span = n_locus + 2 * k
    best = None
    for f in left:
        for c1, r1 in positions.get((strain, f), ()):
            if in_place is not None and (c1, r1) not in in_place:
                continue
            for g in right:
                for c2, r2 in positions.get((strain, g), ()):
                    if c2 != c1 or abs(r1 - r2) > span:
                        continue
                    if in_place is not None and (c2, r2) not in in_place:
                        continue
                    cand = (abs(r1 - r2), c1, min(r1, r2), max(r1, r2))
                    if best is None or cand < best:
                        best = cand
    if best is None:
        return None
    return best[1], best[2], best[3]


def row_class(base: list, n_left: int, n_locus: int, intact: bool,
              empty_frac: float = DEFAULT_EMPTY_FRAC) -> str:
    """full / partial / empty / uninformative (spec section 5), from the
    base states (rescue counts as its base state; contig break as absent).
    Falls to "partial" (the `return "partial"` below) whenever the block is
    neither all in-place nor empty-fraction absent, including a block with
    some column "elsewhere" and none "absent" -- matches the feasibility
    prototype's behaviour; spec section 5's own wording for partial is
    "some present, some absent"."""
    if not intact:
        return "uninformative"
    block = base[n_left:n_left + n_locus]
    if block and sum(b == ABSENT for b in block) / len(block) >= empty_frac:
        return "empty"
    if all(b == IN_PLACE for b in block):
        return "full"
    return "partial"


def breakpoint_track(strain_rows: list[tuple[str, str, str]], n_cols: int) -> list[dict]:
    """Per boundary b (between columns b-1 and b): flank-intact strains whose
    state changes between in place and absent, by species; and strains of any
    class whose state changes between in place and contig break (Ruling R8).
    `strain_rows` = [(species, row_class, codes)]. Only non-zero boundaries."""
    out = []
    for b in range(1, n_cols):
        indel: dict[str, int] = collections.Counter()
        brk = 0
        for species, cls, codes in strain_rows:
            pair = {codes[b - 1], codes[b]}
            placed = bool(pair & IN_PLACE_CODES)
            if placed and ABSENT in pair and cls != "uninformative":
                indel[species] += 1
            if placed and CONTIG_BREAK in pair:
                brk += 1
        if indel or brk:
            out.append({"b": b, "indel": dict(sorted(indel.items())), "contig_break": brk})
    return out


# Display order of row classes. "model_difference" occurs only with the DNA
# presence check (section 6b, spec section 4b).
ROW_ORDER = ("full", "partial", "empty", "model_difference", "variable_gap", "uninformative")


def collapse_rows(per_strain: dict[str, tuple[str, str]],
                  details: dict[str, list] | None = None,
                  max_detail_rows: int = MAX_DETAIL_ROWS) -> list[dict]:
    """Collapse strains with the same (row class, codes) into one row.
    Sorted by ROW_ORDER, then codes. `rep` is the per-column detail
    of the row's first strain, kept for the `max_detail_rows` rows with the
    most strains (None beyond that, to bound the payload); see
    encode_detail() for its shape."""
    groups: dict[tuple[str, str], list[str]] = collections.defaultdict(list)
    for strain, key in per_strain.items():
        groups[key].append(strain)
    rows = []
    for (cls, codes), strains in groups.items():
        strains.sort()
        rows.append({"row_class": cls, "codes": codes, "count": len(strains), "strains": strains,
                     "rep": None})
    rows.sort(key=lambda r: (ROW_ORDER.index(r["row_class"]), r["codes"]))
    if details is not None:
        by_count = sorted(range(len(rows)), key=lambda i: (-rows[i]["count"], i))
        for i in by_count[:max_detail_rows]:
            rep = details.get(rows[i]["strains"][0])
            rows[i]["rep"] = encode_detail(rep) if rep else None
    return rows


def encode_detail(detail: list) -> dict:
    """Compact per-column detail for one strain: {"c": [contig names],
    "d": [cell per column]}. A cell is null (no copy), [contig index, rank]
    (a copy with no neighbour) or [contig index, rank, neighbour column,
    signed rank delta] (in place)."""
    contigs: list[str] = []
    index: dict[str, int] = {}
    cells = []
    for d in detail:
        if d is None:
            cells.append(None)
            continue
        contig, rank, nbr, delta = d
        if contig not in index:
            index[contig] = len(contigs)
            contigs.append(contig)
        cells.append([index[contig], rank] if nbr is None else [index[contig], rank, nbr, delta])
    return {"c": contigs, "d": cells}


def informative_score(counts: dict[str, int]) -> int:
    """Spec section 6: min(empty, full) when >= 10 empty-site and >= 2
    full-locus strains, else -1."""
    empty, full = counts.get("empty", 0), counts.get("full", 0)
    if empty >= INFORMATIVE_MIN_EMPTY and full >= INFORMATIVE_MIN_FULL:
        return min(empty, full)
    return -1


def rank_key(result: dict, rank_by: str) -> tuple:
    """Sort key over computed locus results (dicts with informative_score,
    n_carriers, size, locus_id); best first."""
    if rank_by == "size":
        return (-result["size"], -result["n_carriers"], result["locus_id"])
    if rank_by == "strains":
        return (-result["n_carriers"], -result["size"], result["locus_id"])
    return (-result["informative_score"], -result["n_carriers"], -result["size"], result["locus_id"])


# ---- 6. one locus, all strains -------------------------------------------------
DEFAULT_MIN_COLUMN_STRAINS = 2
DEFAULT_MODEL_DIFF_MIN_FRAC = 0.5  # a DNA-present column is a "model difference" only if its gene is in place in >= this fraction of strains


def informative_columns(support: list[int], n_left: int, n_locus: int,
                        min_strains: int) -> list[int]:
    """1 for each locus column in place in at least `min_strains` strains, else 0.

    A column in place in a single strain (usually the exemplar's own gene model that no
    other strain annotates) carries no information about variation among the other strains,
    but it turns each of them into "model difference" or "partial" for that column. Such
    columns are drawn but not used to classify strains. If no locus column reaches the
    threshold, every column is used (a locus with no shared gene still needs a class).
    """
    block = support[n_left:n_left + n_locus]
    mask = [1 if s >= min_strains else 0 for s in block]
    return mask if any(mask) else [1] * len(block)


def _classify_block(seq, n_left: int, n_locus: int, mask: list[int]):
    """(seq with the uninformative locus columns removed, new n_locus)."""
    kept = [seq[n_left + i] for i in range(n_locus) if mask[i]]
    return list(seq[:n_left]) + kept + list(seq[n_left + n_locus:]), len(kept)


def compute_locus(locus: Locus, placement: Placement, tier: str, columns: Columns,
                  strains: list[str], positions: Positions, spans: ContigSpans,
                  matrix, species_of: dict[str, str], k: int = DEFAULT_K,
                  empty_frac: float = DEFAULT_EMPTY_FRAC,
                  min_column_strains: int = DEFAULT_MIN_COLUMN_STRAINS) -> dict:
    """Every strain's cells, row class and flank pair for one locus, plus the
    collapsed rows, counts and breakpoint track. `matrix` has
    .call(family, strain) -> "present" | "genome_only" | "absent" (a
    lib.pangenome_matrix.PresenceMatrix). Keys starting with "_" are for
    the caller (Part B regions) and never go into the page payload."""
    fams = columns.families
    n_left, n_locus = len(columns.left), len(columns.locus)
    per_strain: dict[str, tuple[str, str]] = {}
    details: dict[str, list] = {}
    cells: dict[str, StrainCells] = {}
    pairs: dict[str, tuple[str, int, int]] = {}
    counts = {c: 0 for c in ROW_CLASSES}
    by_species: dict[str, dict[str, int]] = {}
    track_rows = []
    n_carriers = 0
    counts_all = {c: 0 for c in ROW_CLASSES}
    row_class_all: dict[str, str] = {}
    # Pass 1: every strain's cells, so column support (strains with the column in place) is known.
    for s in strains:
        calls = {f: matrix.call(f, s) for f in set(fams)}
        genome_only = frozenset(f for f, c in calls.items() if c == "genome_only")
        unplaced = frozenset(f for f, c in calls.items()
                             if c != "absent" and not positions.get((s, f)))
        sc = strain_cells(s, fams, positions, spans, k, genome_only, unplaced)
        in_place = {(c, r) for _, c, r in sc.in_place_copies}
        pair = flank_pair(positions, s, list(columns.left), list(columns.right), n_locus, k,
                          in_place=in_place)
        cells[s] = sc
        details[s] = sc.detail
        if pair is not None:
            pairs[s] = pair
    support = [sum(1 for sc in cells.values() if sc.base[j] == IN_PLACE) for j in range(len(fams))]
    mask = informative_columns(support, n_left, n_locus, min_column_strains)
    # Pass 2: classify on the informative columns; the all-columns class is kept for comparison.
    for s in strains:
        sc = cells[s]
        intact = s in pairs
        base_f, n_locus_f = _classify_block(sc.base, n_left, n_locus, mask)
        cls = row_class(base_f, n_left, n_locus_f, intact, empty_frac)
        cls_all = row_class(sc.base, n_left, n_locus, intact, empty_frac)
        codes = "".join(sc.codes)
        per_strain[s] = (cls, codes)
        row_class_all[s] = cls_all
        counts[cls] += 1
        counts_all[cls_all] += 1
        sp = species_of.get(s, "")
        by_species.setdefault(sp, {c: 0 for c in ROW_CLASSES})[cls] += 1
        track_rows.append((sp, cls, codes))
        if any(b == IN_PLACE for b in sc.base[n_left:n_left + n_locus]):
            n_carriers += 1
    return {
        "locus_id": locus.locus_id,
        "size": locus.size,
        "n_variants": len(locus.variants),
        "variant_strains": locus.carrier_proxy(len(strains)),
        "exemplar": placement.strain,
        "exemplar_contig": placement.contig,
        "tier": tier,
        "families": fams,
        "n_left": n_left,
        "n_locus": n_locus,
        "n_right": len(columns.right),
        "counts": counts,
        "counts_all_columns": counts_all,
        "locus_mask": mask,
        "column_strains": support,
        "n_excluded_columns": mask.count(0),
        "min_column_strains": min_column_strains,
        "counts_by_species": dict(sorted(by_species.items())),
        "n_carriers": n_carriers,
        "informative_score": informative_score(counts),
        "rows": collapse_rows(per_strain, details),
        "breakpoints": breakpoint_track(track_rows, len(fams)),
        "_cells": cells,
        "_pairs": pairs,
        "_placement": placement,
        "_row_class": {s: v[0] for s, v in per_strain.items()},
        "_row_class_all": row_class_all,
        "_per_strain": per_strain,
    }


def locus_payload(result: dict, key: str, bins: dict[str, str],
                  family_classes: list[str], family_domains: list[str],
                  dominant: str, family_locations: list | None = None,
                  exemplar_span: tuple[int, int] | None = None,
                  ranks: dict | None = None) -> dict:
    """The JSON-safe page entry for one computed locus (drops "_" keys).
    `ranks` (from locus_ranks()) becomes "ranks" (the five scores) and
    "within_species" (spec 4b/6, ranks-brief.md)."""
    out = {k: v for k, v in result.items() if not k.startswith("_")}
    out["key"] = key
    out["family_bins"] = [bins.get(f, "") for f in result["families"]]
    out["family_classes"] = family_classes
    out["family_domains"] = family_domains
    out["dominant_class"] = dominant
    out["exemplar_span"] = ({"start": exemplar_span[0], "end": exemplar_span[1]}
                            if exemplar_span else None)
    if family_locations is not None:
        out["family_locations"] = family_locations
    out["family_labels"], out["label_kinds"] = column_labels(result["families"], family_locations)
    if ranks is not None:
        out["ranks"] = {k: ranks[k] for k in
                        ("whole_annot", "whole_dna", "presence", "species", "within")}
        out["within_species"] = ranks["within_species"]
    return out


# ---- 6b. DNA presence check (spec section 4b) ---------------------------------
# Two more cell codes. The 3-bit payload codes (spec "Data and wiring") now
# use all eight values 0-7.
DNA_PRESENT = "6"  # absent or elsewhere, but the exemplar gene's DNA is at the site
DNA_ABSENT = "7"  # absent, and the site lacks the exemplar gene's DNA
DNA_STATE_LABELS = {DNA_PRESENT: "absent, DNA present", DNA_ABSENT: "absent, DNA absent"}
# Cells the check can change (spec 4b: an "absent" or "elsewhere" locus cell).
DNA_CHECKED_CODES = frozenset({ABSENT, ELSEWHERE, RESCUE_ELSEWHERE})
DEFAULT_DNA_MIN_ID = 90.0
DEFAULT_DNA_MIN_COV = 80.0
DNA_MIN_TARGET_BP = 50  # guard only: a target holds two flank genes, so it is never this short
# Breakpoint track with the check on: only in place <-> DNA absent counts
# (spec 4b), so DNA absent is fed to breakpoint_track() as ABSENT and every
# other not-in-place code as ELSEWHERE, which the track never counts.
_DNA_TRACK = str.maketrans({ABSENT: ELSEWHERE, DNA_PRESENT: ELSEWHERE, DNA_ABSENT: ABSENT})


def dna_checked_strains(result: dict) -> list[str]:
    """Strains the DNA check covers for one computed locus (spec 4b): flank
    intact, with at least one CLASSIFYING locus column not in place. Columns that
    are not used to classify strains (locus_mask 0, see informative_columns) do not
    trigger a check: a strain whose only gap is an exemplar-private column is already a
    full locus. Sorted."""
    nl, nb = result["n_left"], result["n_locus"]
    mask = result.get("locus_mask") or [1] * nb
    return sorted(s for s in result["_pairs"]
                  if any(c not in IN_PLACE_CODES and m
                         for c, m in zip(result["_cells"][s].codes[nl:nl + nb], mask)))


def copy_bp(positions: Positions, gene_locs: dict, strain: str, family: str,
            contig: str, rank: int, rescue_spans: dict | None = None) -> tuple[int, int] | None:
    """(start, end) in bp of the copy of `family` at (contig, rank).

    A family's copies on one contig in rank order are its gene_positions
    copies, plus its TBLASTN rescue hits in `rescue_spans` ({(strain,
    family): [(contig, start, end)]}, Ruling R26), in start order, because
    ranks enumerate genes and rescue hits sorted by (contig, start) (Ruling
    R19). None when the copy has no span (a rescue hit not in
    `rescue_spans`) or the copy counts differ."""
    ranks = sorted(r for c, r in positions.get((strain, family), ()) if c == contig)
    spans = [(s, e) for _pid, c, s, e in gene_locs.get((strain, family), ()) if c == contig]
    spans += [(s, e) for c, s, e in (rescue_spans or {}).get((strain, family), ()) if c == contig]
    spans.sort()
    if rank not in ranks or len(ranks) != len(spans):
        return None
    return spans[ranks.index(rank)]


def rescue_hit_spans(rescue_rows, tblastn_lines) -> dict[tuple[str, str], list[tuple]]:
    """{(strain, family): [(contig, start, end)]} for TBLASTN rescue hits
    (Ruling R26). `rescue_rows` are rescue_positions.tsv rows (strain,
    family, contig, start), where start is min(sstart, send) of the chosen
    HSP; `tblastn_lines` are the strains' own tblastn lines (outfmt `6 std
    qcovs`; subject = `strain|contig`). The end is max(sstart, send) of the
    HSP of that family on that contig that starts at the recorded start (the
    highest bitscore if several). A row with no such HSP has no span."""
    want = {(s, f, c, int(st)) for s, f, c, st in rescue_rows}
    fams = {f for _s, f, _c, _st in want}
    best: dict[tuple, tuple[float, int]] = {}
    for line in tblastn_lines:
        p = line.rstrip("\n").split("\t")
        if len(p) < 12 or p[0] not in fams:
            continue
        strain, _, contig = p[1].partition("|")
        try:
            s0, s1, bits = int(p[8]), int(p[9]), float(p[11])
        except ValueError:
            continue
        key = (strain, p[0], contig, min(s0, s1))
        if key in want and (key not in best or bits > best[key][0]):
            best[key] = (bits, max(s0, s1))
    out: dict[tuple[str, str], list[tuple]] = {}
    for (s, f, c, st), (_bits, end) in sorted(best.items()):
        out.setdefault((s, f), []).append((c, st, end))
    return out


def dna_query(result: dict, exemplar_ranks: list[int], positions: Positions,
              gene_locs: dict, rescue_spans: dict | None = None) -> tuple[str, int, int, list] | None:
    """The exemplar's locus DNA (spec 4b): (contig, start, end, genes) with
    genes = [(column, gene start, gene end)], from the first locus gene's
    start to the last one's end. `exemplar_ranks` are the ranks of the locus
    columns in the exemplar, in column order. A column whose exemplar copy
    is a TBLASTN rescue hit uses the hit's span from `rescue_spans` (spec
    4b, Ruling R26). A column with no span at all is left out and stays
    unchecked (Ruling R21). None if no column has a span."""
    nl = result["n_left"]
    fams = result["families"]
    genes = []
    for i, rank in enumerate(exemplar_ranks):
        bp = copy_bp(positions, gene_locs, result["exemplar"], fams[nl + i],
                     result["exemplar_contig"], rank, rescue_spans)
        if bp is not None:
            genes.append((nl + i, bp[0], bp[1]))
    if not genes:
        return None
    return (result["exemplar_contig"], min(g[1] for g in genes), max(g[2] for g in genes), genes)


def dna_target(result: dict, strain: str, positions: Positions, gene_locs: dict,
               k: int = DEFAULT_K) -> tuple[str, int, int] | None:
    """The strain's target DNA (spec 4b): (contig, start, end), 1-based,
    from the start of its innermost in-place left-flank gene to the end of
    its innermost in-place right-flank gene, on its flank-pair contig. The
    two flank genes are included, because a gene model can extend over the
    locus DNA. Only copies within n_locus + 2k ranks of the pair count (as
    in Ruling R11), and only annotated copies (a bp position is needed).
    Innermost = the flank column nearest the locus; ties between copies of
    that column go to the copy nearest the pair. None when either side has
    no annotated in-place copy (Ruling R20)."""
    nl, nb = result["n_left"], result["n_locus"]
    fams = result["families"]
    contig, p_lo, p_hi = result["_pairs"][strain]
    w = nb + 2 * k
    mid = (p_lo + p_hi) / 2
    left, right = [], []
    for ci, c, r in result["_cells"][strain].in_place_copies:
        if c != contig or not p_lo - w <= r <= p_hi + w or nl <= ci < nl + nb:
            continue
        bp = copy_bp(positions, gene_locs, strain, fams[ci], c, r)
        if bp is not None:
            (left if ci < nl else right).append((ci, abs(r - mid), bp))
    if not left or not right:
        return None
    a = min(left, key=lambda x: (-x[0], x[1]))[2]
    b = min(right, key=lambda x: (x[0], x[1]))[2]
    return contig, min(a[0], b[0]), max(a[1], b[1])


def merge_intervals(intervals) -> list[tuple[int, int]]:
    """Union of closed integer intervals, sorted."""
    out: list[list[int]] = []
    for lo, hi in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(lo, hi) for lo, hi in out]


def gene_coverage(hsps, gene: tuple[int, int], min_id: float = DEFAULT_DNA_MIN_ID) -> float:
    """Fraction of `gene` (query coordinates, closed) covered by the union of
    HSPs [(qstart, qend, pident)] with pident >= min_id (spec 4b: HSPs are
    merged where they overlap, so two HSPs over one gene add up)."""
    g0, g1 = gene
    kept = merge_intervals((q0, q1) for q0, q1, pid in hsps if pid >= min_id)
    covered = sum(max(0, min(hi, g1) - max(lo, g0) + 1) for lo, hi in kept)
    return covered / (g1 - g0 + 1)


def row_class_dna(codes, n_left: int, n_locus: int, intact: bool,
                  empty_frac: float = DEFAULT_EMPTY_FRAC, common=None) -> str:
    """Row class of a DNA-checked strain (spec 4b): empty site = >= empty_frac
    of the locus columns DNA absent; full = all in place; model difference =
    every column in place or DNA present, at least one DNA present; anything
    else is partial (Ruling R21). `common` (one bool per locus column, or None)
    marks the columns whose gene is in place in many strains. When given, a
    strain whose DNA-present columns are all rare ones is "variable_gap"
    (a rare gene lacking in this annotation), not "model_difference"
    (2026-10-08)."""
    if not intact:
        return "uninformative"
    block = codes[n_left:n_left + n_locus]
    if block and sum(c == DNA_ABSENT for c in block) / len(block) >= empty_frac:
        return "empty"
    if all(c in IN_PLACE_CODES for c in block):
        return "full"
    if all(c in IN_PLACE_CODES or c == DNA_PRESENT for c in block):
        if common is not None and not any(c == DNA_PRESENT and common[i] for i, c in enumerate(block)):
            return "variable_gap"
        return "model_difference"
    return "partial"


def apply_dna_calls(result: dict, calls: dict[str, dict[int, str]], species_of: dict[str, str],
                    empty_frac: float = DEFAULT_EMPTY_FRAC,
                    model_diff_min_frac: float = DEFAULT_MODEL_DIFF_MIN_FRAC) -> None:
    """Apply one locus's DNA calls in place (spec 4b).

    `calls` = {strain: {column: "present" | "absent" | "unchecked"}}. In a
    strain with at least one present/absent call, each ABSENT / ELSEWHERE /
    RESCUE_ELSEWHERE locus cell with a call becomes DNA_PRESENT or
    DNA_ABSENT, and the row class comes from row_class_dna(). Other strains
    keep their section 4-5 states. Recomputes counts (with a
    "model_difference" count), counts_by_species, rows, the breakpoint track
    (in place <-> DNA absent only) and the informative score, which counts
    only DNA-confirmed empty-site strains (Ruling R22). Adds `dna` =
    {checked, unchecked, empty_confirmed, empty_to_model_difference}."""
    nl, nb = result["n_left"], result["n_locus"]
    mask = result.get("locus_mask") or [1] * nb
    support = result.get("column_strains")
    common = None
    if support and model_diff_min_frac > 0 and result["_cells"]:
        floor = model_diff_min_frac * len(result["_cells"])
        common = [support[nl + i] >= floor for i in range(nb)]
    common_f = [c for c, m in zip(common, mask) if m] if common else None
    need = set(dna_checked_strains(result))
    per_strain: dict[str, tuple[str, str]] = {}
    details: dict[str, list] = {}
    counts = {c: 0 for c in ROW_ORDER}
    counts_all = {c: 0 for c in ROW_ORDER}
    row_class_all: dict[str, str] = {}
    by_species: dict[str, dict[str, int]] = {}
    track = []
    summary = {"checked": 0, "unchecked": 0, "empty_confirmed": 0, "empty_to_model_difference": 0}
    for s, sc in result["_cells"].items():
        cls = result["_row_class"][s]
        cls_all = result.get("_row_class_all", result["_row_class"]).get(s, cls)
        codes = list(sc.codes)
        col_calls = calls.get(s, {}) if s in need else {}
        if any(v in ("present", "absent") for v in col_calls.values()):
            for ci, v in col_calls.items():
                if nl <= ci < nl + nb and codes[ci] in DNA_CHECKED_CODES and v in ("present", "absent"):
                    codes[ci] = DNA_PRESENT if v == "present" else DNA_ABSENT
            codes_f, nb_f = _classify_block(codes, nl, nb, mask)
            new_cls = row_class_dna(codes_f, nl, nb_f, True, empty_frac, common_f)
            cls_all = row_class_dna(codes, nl, nb, True, empty_frac, common)
            summary["checked"] += 1
            summary["empty_confirmed"] += new_cls == "empty"
            summary["empty_to_model_difference"] += cls == "empty" and new_cls == "model_difference"
            cls = new_cls
        elif s in need:
            summary["unchecked"] += 1
        code_str = "".join(codes)
        per_strain[s] = (cls, code_str)
        details[s] = sc.detail
        counts[cls] += 1
        counts_all[cls_all] += 1
        row_class_all[s] = cls_all
        sp = species_of.get(s, "")
        by_species.setdefault(sp, {c: 0 for c in ROW_ORDER})[cls] += 1
        track.append((sp, cls, code_str.translate(_DNA_TRACK)))
    result["counts"] = counts
    result["counts_all_columns"] = counts_all
    result["_row_class_all"] = row_class_all
    result["counts_by_species"] = dict(sorted(by_species.items()))
    result["rows"] = collapse_rows(per_strain, details)
    result["breakpoints"] = breakpoint_track(track, len(result["families"]))
    result["informative_score"] = informative_score(
        {"empty": summary["empty_confirmed"], "full": counts["full"]})
    result["dna"] = summary
    result["_row_class"] = {s: v[0] for s, v in per_strain.items()}
    result["_per_strain"] = per_strain


# ---- 6c. Multi-rank scores (ranks-brief.md, changed 2026-09-27) --------------
# Five rankings replace the single "informative" ranking (spec section 6,
# 4b "Ranking"). "presence" generalises the old informative_score() from
# empty/full counts to a carriers/losses partition of every non-uninformative
# strain, so an existing empty_frac-partial strain with a missing locus
# column counts as a loss too, not only fully "empty" strains.
DEFAULT_FIXED_DIFF = 0.95
DEFAULT_POLY_MIN_STRAINS = 20
DEFAULT_POLY_MIN_FRAC = 0.05
DEFAULT_POLY_MAX_FRAC = 0.95


def locus_species_stats(per_strain: dict[str, tuple[str, str]], species_of: dict[str, str],
                        n_left: int, n_locus: int, missing_code: str) -> dict:
    """carriers/losses over one locus's non-uninformative strains, overall and
    per species (ranks-brief.md "Definitions"). `missing_code` is DNA_ABSENT
    with the DNA check on, ABSENT otherwise -- the code a partial strain's
    locus block must NOT contain to still count as a carrier."""
    carriers = losses = 0
    by_species: dict[str, dict[str, int]] = {}
    for s, (cls, codes) in per_strain.items():
        if cls == "uninformative":
            continue
        if cls in ("full", "model_difference", "variable_gap"):
            is_carrier = True
        elif cls == "empty":
            is_carrier = False
        else:
            is_carrier = missing_code not in codes[n_left:n_left + n_locus]
        sp = species_of.get(s, "")
        d = by_species.setdefault(sp, {"car": 0, "loss": 0})
        if is_carrier:
            carriers += 1
            d["car"] += 1
        else:
            losses += 1
            d["loss"] += 1
    return {"carriers": carriers, "losses": losses, "by_species": by_species}


def score_whole_annot(empty_n: int, full: int) -> int:
    """Rank A (spec 4b/6): min(empty_n, full) if empty_n >= 10 and full >= 2,
    else -1. Identical to informative_score() given the same counts."""
    if empty_n >= INFORMATIVE_MIN_EMPTY and full >= INFORMATIVE_MIN_FULL:
        return min(empty_n, full)
    return -1


def score_whole_dna(empty_n: int, full: int, model_difference: int) -> int:  # pass model_difference + variable_gap
    """Rank B: like A, but the "full-locus-like" side also counts model
    difference strains (annotated carriers plus a gene-model difference)."""
    return score_whole_annot(empty_n, full + model_difference)


def score_presence(carriers: int, losses: int) -> int:
    """Rank C (the page default): min(losses, carriers) if losses >= 10 and
    carriers >= 2, else -1."""
    if losses >= INFORMATIVE_MIN_EMPTY and carriers >= INFORMATIVE_MIN_FULL:
        return min(losses, carriers)
    return -1


def score_species(by_species: dict[str, dict[str, int]], n_species_total: int,
                  fixed_diff: float = DEFAULT_FIXED_DIFF, min_n: int = 10) -> float | int:
    """Rank D: max_s f_s - min_s f_s over species with n_s >= min_n, when at
    least 2 such species exist and the difference is >= fixed_diff; else -1.
    Always -1 when fewer than 2 species exist in the samplesheet at all."""
    if n_species_total < 2:
        return -1
    fracs = []
    for d in by_species.values():
        n = d["car"] + d["loss"]
        if n >= min_n:
            fracs.append(d["loss"] / n)
    if len(fracs) < 2:
        return -1
    diff = max(fracs) - min(fracs)
    return round(diff, 3) if diff >= fixed_diff else -1


def score_within(by_species: dict[str, dict[str, int]],
                 poly_min_strains: int = DEFAULT_POLY_MIN_STRAINS,
                 poly_min_frac: float = DEFAULT_POLY_MIN_FRAC,
                 poly_max_frac: float = DEFAULT_POLY_MAX_FRAC) -> tuple[int, str | None]:
    """Rank E: (score, species) where score = the largest min(car_s, loss_s)
    over species with n_s >= poly_min_strains and a loss fraction inside
    [poly_min_frac, poly_max_frac]; (-1, None) if none qualify. Ties go to
    the alphabetically first species name."""
    best = None
    for sp, d in sorted(by_species.items()):
        n = d["car"] + d["loss"]
        if n < poly_min_strains:
            continue
        f = d["loss"] / n
        if not (poly_min_frac <= f <= poly_max_frac):
            continue
        score = min(d["car"], d["loss"])
        if best is None or score > best[1]:
            best = (sp, score)
    return (best[1], best[0]) if best else (-1, None)


def locus_ranks(result: dict, species_of: dict[str, str], dna_on: bool, n_species_total: int,
                fixed_diff: float = DEFAULT_FIXED_DIFF,
                poly_min_strains: int = DEFAULT_POLY_MIN_STRAINS,
                poly_min_frac: float = DEFAULT_POLY_MIN_FRAC,
                poly_max_frac: float = DEFAULT_POLY_MAX_FRAC) -> dict:
    """All five rank scores for one computed locus (ranks-brief.md), plus the
    carriers/losses totals used for the shared tie-break (higher strain count,
    then locus_id) and within_species (rank E's species, or None)."""
    nl, nb = result["n_left"], result["n_locus"]
    counts = result["counts"]
    empty_n = result["dna"]["empty_confirmed"] if dna_on and "dna" in result else counts["empty"]
    full = counts["full"]
    model_difference = counts.get("model_difference", 0) + counts.get("variable_gap", 0)
    missing_code = DNA_ABSENT if dna_on else ABSENT
    stats = locus_species_stats(result["_per_strain"], species_of, nl, nb, missing_code)
    within_score, within_species = score_within(stats["by_species"], poly_min_strains,
                                                poly_min_frac, poly_max_frac)
    return {
        "whole_annot": score_whole_annot(empty_n, full),
        "whole_dna": score_whole_dna(empty_n, full, model_difference),
        "presence": score_presence(stats["carriers"], stats["losses"]),
        "species": score_species(stats["by_species"], n_species_total, fixed_diff),
        "within": within_score,
        "within_species": within_species,
        "carriers": stats["carriers"],
        "losses": stats["losses"],
    }


def rank_sort_key(ranks: dict, key: str, locus_id: str) -> tuple:
    """Sort key for one rank score, best first: score descending, then
    strain count (carriers + losses) descending, then locus_id (ranks-brief.md
    "Ties inside a rank")."""
    return (-ranks[key], -(ranks["carriers"] + ranks["losses"]), locus_id)


# ---- 8. clinker panel: strains and regions (spec section 8) --------------------
DEFAULT_CLINKER_MAX_STRAINS = 12
CLINKER_CLASSES = ("full", "partial", "empty", "model_difference", "variable_gap")


def strain_region(result: dict, strain: str, spans: ContigSpans, flank: int,
                  k: int = DEFAULT_K) -> dict | None:
    """The rank range to draw for one strain (spec section 8, Ruling R11).

    On the contig of the strain's flank pair: from the lowest to the highest
    in-place column copy within n_locus + 2k ranks of the pair, then `flank`
    genes further on each side, clipped at the contig ends. The exemplar
    without a flank pair uses its Placement block instead. None for any
    other strain without a flank pair (it cannot be anchored)."""
    cells = result["_cells"][strain]
    pair = result["_pairs"].get(strain)
    fams = result["families"]
    if pair is not None:
        contig, p_lo, p_hi = pair
        w = result["n_locus"] + 2 * k
        anchors = sorted((r, fams[ci]) for ci, c, r in cells.in_place_copies
                         if c == contig and p_lo - w <= r <= p_hi + w)
        lo, hi = anchors[0][0], anchors[-1][0]
    elif strain == result["exemplar"]:
        place = result["_placement"]
        contig, lo, hi = place.contig, place.lo, place.hi
        anchors = sorted((r, fams[ci]) for ci, c, r in cells.in_place_copies
                         if c == contig and lo <= r <= hi)
    else:
        return None
    cmin, cmax = spans.get((strain, contig), (lo, hi))
    return {"contig": contig, "rank_lo": max(cmin, lo - flank), "rank_hi": min(cmax, hi + flank),
            "anchors": anchors}


def select_clinker_strains(result: dict, n50: dict[str, int], species_of: dict[str, str],
                           spans: ContigSpans, flank: int, k: int = DEFAULT_K,
                           max_strains: int = DEFAULT_CLINKER_MAX_STRAINS) -> list[dict]:
    """Up to `max_strains` strains for the clinker figure, in spec order:
    the exemplar; per row class (full, partial, empty, model difference)
    and species the flank-intact strain with the best quality_key(); then
    more full or partial strains by quality_key() until the cap.
    Uninformative strains
    are never chosen (the exemplar is kept even when its flanks are not
    intact). Each pick carries its reason and its strain_region()."""
    pairs = result["_pairs"]
    row_cls = result["_row_class"]

    def qkey(s):
        contig = pairs[s][0] if s in pairs else result["exemplar_contig"]
        cmin, cmax = spans.get((s, contig), (0, -1))
        return quality_key(s, n50, cmax - cmin + 1)

    picks: list[dict] = []
    chosen: set[str] = set()

    def add(strain, reason):
        if len(picks) >= max_strains or strain in chosen:
            return
        region = strain_region(result, strain, spans, flank, k)
        if region is None:
            return
        chosen.add(strain)
        picks.append({"strain": strain, "reason": reason, "row_class": row_cls[strain],
                      "species": species_of.get(strain, ""), **region})

    add(result["exemplar"], "exemplar")
    species = sorted({species_of.get(s, "") for s in pairs})
    for cls in CLINKER_CLASSES:
        for sp in species:
            pool = [s for s in pairs if row_cls[s] == cls and species_of.get(s, "") == sp
                    and s not in chosen]
            if pool:
                add(min(pool, key=qkey), "best_" + cls)
    fill = sorted((s for s in pairs if row_cls[s] in ("full", "partial") and s not in chosen),
                  key=qkey)
    for s in fill:
        if len(picks) >= max_strains:
            break
        add(s, "fill")
    return picks
