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
    base states (rescue counts as its base state; contig break as absent)."""
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
ROW_ORDER = ("full", "partial", "empty", "model_difference", "uninformative")


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
def compute_locus(locus: Locus, placement: Placement, tier: str, columns: Columns,
                  strains: list[str], positions: Positions, spans: ContigSpans,
                  matrix, species_of: dict[str, str], k: int = DEFAULT_K,
                  empty_frac: float = DEFAULT_EMPTY_FRAC) -> dict:
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
    for s in strains:
        calls = {f: matrix.call(f, s) for f in set(fams)}
        genome_only = frozenset(f for f, c in calls.items() if c == "genome_only")
        unplaced = frozenset(f for f, c in calls.items()
                             if c != "absent" and not positions.get((s, f)))
        sc = strain_cells(s, fams, positions, spans, k, genome_only, unplaced)
        in_place = {(c, r) for _, c, r in sc.in_place_copies}
        pair = flank_pair(positions, s, list(columns.left), list(columns.right), n_locus, k,
                          in_place=in_place)
        cls = row_class(sc.base, n_left, n_locus, pair is not None, empty_frac)
        codes = "".join(sc.codes)
        per_strain[s] = (cls, codes)
        details[s] = sc.detail
        cells[s] = sc
        if pair is not None:
            pairs[s] = pair
        counts[cls] += 1
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
        "counts_by_species": dict(sorted(by_species.items())),
        "n_carriers": n_carriers,
        "informative_score": informative_score(counts),
        "rows": collapse_rows(per_strain, details),
        "breakpoints": breakpoint_track(track_rows, len(fams)),
        "_cells": cells,
        "_pairs": pairs,
        "_row_class": {s: v[0] for s, v in per_strain.items()},
    }


def locus_payload(result: dict, key: str, bins: dict[str, str],
                  family_classes: list[str], family_domains: list[str],
                  dominant: str, family_locations: list | None = None,
                  exemplar_span: tuple[int, int] | None = None) -> dict:
    """The JSON-safe page entry for one computed locus (drops "_" keys)."""
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
    return out
