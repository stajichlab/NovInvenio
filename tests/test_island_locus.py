"""Unit tests for lib/island_locus.py on synthetic strains (spec validation
plan item 1). Positions are (contig, rank); a contig's ranks run from its
min to its max rank (see the module docstring). Each section imports the
names it adds, so the file grows task by task."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import Locus, candidate_loci, group_loci  # noqa: E402


def isl(locus_id, members, n_strains=2, example="S1", contig="c1"):
    return {"locus_id": locus_id, "member_families": ",".join(members),
            "n_strains": str(n_strains), "example_strain": example, "locus_contig": contig}


# ---- group_loci -------------------------------------------------------------

def test_island_with_half_its_families_in_a_larger_island_joins_its_locus():
    loci = group_loci([isl("big", ["a", "b", "c", "d"]), isl("small", ["a", "x"])])
    assert [loc.locus_id for loc in loci] == ["big"]
    assert [v["locus_id"] for v in loci[0].variants] == ["big", "small"]


def test_island_below_containment_starts_its_own_locus():
    loci = group_loci([isl("big", ["a", "b", "c", "d"]), isl("small", ["a", "x", "y"])])
    assert sorted(loc.locus_id for loc in loci) == ["big", "small"]


def test_equal_size_islands_never_merge():
    loci = group_loci([isl("one", ["a", "b"]), isl("two", ["a", "b"])])
    assert len(loci) == 2


def test_grouping_is_transitive_to_the_root():
    rows = [isl("root", ["a", "b", "c", "d", "e", "f"]), isl("mid", ["a", "b", "c", "z"]),
            isl("leaf", ["z", "q"])]
    loci = group_loci(rows)
    assert [loc.locus_id for loc in loci] == ["root"]
    assert [v["locus_id"] for v in loci[0].variants] == ["root", "mid", "leaf"]


def test_tie_between_two_larger_islands_goes_to_the_earlier_one():
    rows = [isl("first", ["a", "b", "c"]), isl("second", ["x", "y", "z"]), isl("pair", ["a", "x"])]
    loci = group_loci(rows)
    by_id = {loc.locus_id: [v["locus_id"] for v in loc.variants] for loc in loci}
    assert by_id == {"first": ["first", "pair"], "second": ["second"]}


def test_unplaced_islands_are_not_grouped():
    assert group_loci([isl("-", ["a", "b"])]) == []


def test_carrier_proxy_is_capped_at_the_strain_count():
    loc = Locus(root={"members": ["a"], "n_strains": "8"},
                variants=[{"n_strains": "8"}, {"n_strains": "5"}])
    assert loc.carrier_proxy(10) == 10


def test_candidates_filter_on_min_strains_and_rank_by_proxy():
    loci = group_loci([isl("L1", ["a", "b"], 1), isl("L2", ["c", "d"], 9), isl("L3", ["e", "f", "g"], 3)])
    assert [loc.locus_id for loc in candidate_loci(loci, 20, "informative", 5, 2)] == ["L2", "L3"]
    assert [loc.locus_id for loc in candidate_loci(loci, 20, "size", 5, 2)] == ["L3", "L2"]


from island_locus import Placement, carrier_placements, choose_exemplar, quality_key  # noqa: E402


# ---- exemplar -----------------------------------------------------------------

SPANS = {("S1", "c1"): (0, 19), ("S2", "c1"): (0, 19), ("S3", "c2"): (100, 104)}


def test_carrier_needs_every_member_within_the_span_window():
    pos = {("S1", "a"): [("c1", 8)], ("S1", "b"): [("c1", 9)],
           ("S2", "a"): [("c1", 1)], ("S2", "b"): [("c1", 19)]}
    places = carrier_placements(["a", "b"], pos, ["S1", "S2"], SPANS, k=10)
    assert [(p.strain, p.lo, p.hi, p.left_avail, p.right_avail) for p in places] == [("S1", 8, 9, 8, 10)]


def test_strain_missing_a_member_is_not_a_carrier():
    pos = {("S1", "a"): [("c1", 8)]}
    assert carrier_placements(["a", "b"], pos, ["S1"], SPANS) == []


def test_tightest_window_wins_for_a_multicopy_member():
    pos = {("S1", "a"): [("c1", 2), ("c1", 12)], ("S1", "b"): [("c1", 13)]}
    (p,) = carrier_placements(["a", "b"], pos, ["S1"], SPANS)
    assert (p.lo, p.hi) == (12, 13)


def place(strain, left, right, genes=20):
    return Placement(strain, "c1", left, left + 1, left, right, genes)


def test_exemplar_prefers_full_flanks_then_n50_then_name():
    picks = [place("B", 6, 6), place("A", 6, 6), place("C", 2, 9)]
    assert choose_exemplar(picks, {"A": 10, "B": 50}) == (picks[0], "full")
    assert choose_exemplar(picks, {"A": 50, "B": 50})[0].strain == "A"


def test_exemplar_falls_back_to_short_flanks():
    picks = [place("A", 4, 3), place("B", 9, 1)]
    assert choose_exemplar(picks, {}, flank=5, flank_min=3) == (picks[0], "short_flanks")


def test_exemplar_at_contig_end_takes_the_most_flank_genes():
    picks = [place("A", 0, 2), place("B", 1, 2)]
    assert choose_exemplar(picks, {}) == (picks[1], "contig_end")


def test_strains_without_n50_rank_after_known_n50_then_by_contig_genes():
    assert quality_key("X", {"Y": 1}, 500) > quality_key("Y", {"Y": 1}, 10)
    assert quality_key("X", {}, 500) < quality_key("W", {}, 100)


def test_no_carrier_gives_no_exemplar():
    assert choose_exemplar([], {}) is None


from island_locus import locus_columns  # noqa: E402


# ---- columns ------------------------------------------------------------------

ORDER = [(r, f"g{r}") for r in range(20)]


def test_columns_take_the_block_and_flank_genes_each_side():
    cols = locus_columns(ORDER, ["g9", "g10"], 5, block=(9, 10))
    assert cols.left == ("g4", "g5", "g6", "g7", "g8")
    assert cols.locus == ("g9", "g10")
    assert cols.right == ("g11", "g12", "g13", "g14", "g15")
    assert (cols.left_avail, cols.right_avail) == (9, 9)


def test_columns_clip_at_the_contig_end():
    cols = locus_columns(ORDER, ["g1"], 5, block=(1, 1))
    assert cols.left == ("g0",) and len(cols.right) == 5


def test_block_holds_every_gene_between_its_ends_even_non_members():
    cols = locus_columns(ORDER, ["g9", "g11"], 1, block=(9, 11))
    assert cols.locus == ("g9", "g10", "g11")


def test_without_block_every_member_copy_on_the_contig_sets_the_block():
    cols = locus_columns(ORDER, ["g3", "g15"], 1)
    assert cols.locus[0] == "g3" and cols.locus[-1] == "g15"


def test_no_member_on_the_contig_gives_none():
    assert locus_columns(ORDER, ["zz"], 5) is None


from island_locus import (  # noqa: E402
    ABSENT, CONTIG_BREAK, ELSEWHERE, IN_PLACE, RESCUE_ELSEWHERE, RESCUE_IN_PLACE, strain_cells,
)


# ---- cell states ------------------------------------------------------------------

COLS = ["L1", "A", "B", "R1"]


def test_in_place_needs_another_column_within_k_on_the_same_contig():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "A"): [("c1", 6)], ("S1", "B"): [("c2", 6)],
           ("S1", "R1"): [("c1", 30)]}
    sc = strain_cells("S1", COLS, pos, {}, k=10)
    assert sc.codes == [IN_PLACE, IN_PLACE, ELSEWHERE, ELSEWHERE]
    assert sc.detail[0] == ("c1", 5, 1, 1)
    assert sc.detail[2] == ("c2", 6, None, None)


def test_absent_when_no_copy():
    sc = strain_cells("S1", COLS, {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 7)]}, {})
    assert sc.codes == [IN_PLACE, ABSENT, ABSENT, IN_PLACE]


def test_a_genome_only_family_is_drawn_as_rescue():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "A"): [("c1", 6)], ("S1", "B"): [("c9", 1)]}
    sc = strain_cells("S1", COLS, pos, {}, genome_only=frozenset({"A", "B"}))
    assert sc.codes[:3] == [IN_PLACE, RESCUE_IN_PLACE, RESCUE_ELSEWHERE]
    assert sc.base[:3] == [IN_PLACE, IN_PLACE, ELSEWHERE]


def test_present_without_a_position_is_elsewhere_but_base_absent():
    sc = strain_cells("S1", COLS, {}, {}, present_unplaced=frozenset({"A"}))
    assert sc.codes[1] == ELSEWHERE and sc.base[1] == ABSENT


def test_contig_break_beyond_the_last_anchor_near_a_contig_end():
    pos = {("S1", "L1"): [("c1", 97)], ("S1", "A"): [("c1", 98)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)}, k=10)
    assert sc.codes == [IN_PLACE, IN_PLACE, CONTIG_BREAK, CONTIG_BREAK]
    assert sc.base[2:] == [ABSENT, ABSENT]


def test_absence_between_anchors_is_never_a_contig_break():
    pos = {("S1", "L1"): [("c1", 97)], ("S1", "R1"): [("c1", 98)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)})
    assert sc.codes == [IN_PLACE, ABSENT, ABSENT, IN_PLACE]


def test_absence_far_from_a_contig_end_is_absent():
    pos = {("S1", "L1"): [("c1", 50)], ("S1", "A"): [("c1", 51)]}
    sc = strain_cells("S1", COLS, pos, {("S1", "c1"): (0, 99)})
    assert sc.codes[2:] == [ABSENT, ABSENT]


def test_a_tandem_paralog_is_not_its_own_neighbour():
    # Review Focus 1: family P fills two columns (two copies in the exemplar).
    # A strain with one lone copy of P must not be "in place" because the
    # same copy sits in both columns.
    cols = ["L1", "P", "P", "R1"]
    sc = strain_cells("S1", cols, {("S1", "P"): [("c1", 40)]}, {})
    assert sc.codes == [ABSENT, ELSEWHERE, ELSEWHERE, ABSENT]


from island_locus import (  # noqa: E402
    breakpoint_track, collapse_rows, flank_pair, informative_score, rank_key, row_class,
)


# ---- flanks, row classes, breakpoints ------------------------------------------------

def test_flank_pair_on_one_contig_within_locus_plus_2k():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 27)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, k=10) == ("c1", 5, 27)
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 1, k=10) is None


def test_flank_pair_on_different_contigs_is_not_intact():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c2", 6)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2) is None


def test_flank_pair_can_require_in_place_copies():
    pos = {("S1", "L1"): [("c1", 5)], ("S1", "R1"): [("c1", 6)]}
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, in_place={("c1", 5)}) is None
    assert flank_pair(pos, "S1", ["L1"], ["R1"], 2, in_place={("c1", 5), ("c1", 6)}) == ("c1", 5, 6)


def test_row_classes():
    full = [IN_PLACE] * 4
    assert row_class(full, 1, 2, True) == "full"
    assert row_class([IN_PLACE, IN_PLACE, ABSENT, IN_PLACE], 1, 2, True) == "partial"
    assert row_class([IN_PLACE, ABSENT, ABSENT, IN_PLACE], 1, 2, True) == "empty"
    assert row_class(full, 1, 2, False) == "uninformative"


def test_empty_site_threshold_is_80_percent_of_locus_columns():
    four_of_five = [IN_PLACE] + [ABSENT] * 4 + [IN_PLACE] + [IN_PLACE]
    three_of_five = [IN_PLACE] + [ABSENT] * 3 + [IN_PLACE] * 3
    assert row_class(four_of_five, 1, 5, True, 0.8) == "empty"
    assert row_class(three_of_five, 1, 5, True, 0.8) == "partial"


def test_elsewhere_block_is_partial_not_full():
    assert row_class([IN_PLACE, ELSEWHERE, IN_PLACE], 1, 1, True) == "partial"


def test_breakpoint_track_counts_flank_intact_indels_by_species_and_contig_breaks():
    rows = [("sp1", "empty", "1001"), ("sp1", "full", "1111"), ("sp2", "partial", "1101"),
            ("sp2", "uninformative", "1001"), ("sp1", "uninformative", "1155")]
    track = breakpoint_track(rows, 4)
    assert track == [
        {"b": 1, "indel": {"sp1": 1}, "contig_break": 0},
        {"b": 2, "indel": {"sp2": 1}, "contig_break": 1},
        {"b": 3, "indel": {"sp1": 1, "sp2": 1}, "contig_break": 0},
    ]


def test_collapse_groups_same_class_and_codes_and_orders_by_class():
    rows = collapse_rows({"a": ("empty", "100"), "b": ("full", "111"), "c": ("full", "111"),
                          "d": ("uninformative", "111")},
                         details={"b": [("c1", 1, 1, 1), None, ("c1", 3, None, None)]})
    assert [(r["row_class"], r["count"], r["strains"]) for r in rows] == [
        ("full", 2, ["b", "c"]), ("empty", 1, ["a"]), ("uninformative", 1, ["d"])]
    assert rows[0]["rep"] == {"c": ["c1"], "d": [[0, 1, 1, 1], None, [0, 3]]}


def test_detail_is_kept_only_for_the_largest_rows():
    per = {f"s{i}": ("full", format(i, "03b")) for i in range(4)}
    per["s0b"] = ("full", "000")
    rows = collapse_rows(per, details={s: [None, None, None] for s in per}, max_detail_rows=1)
    assert [r["rep"] is not None for r in rows] == [True, False, False, False]


def test_informative_score_needs_10_empty_and_2_full():
    assert informative_score({"empty": 10, "full": 2}) == 2
    assert informative_score({"empty": 9, "full": 50}) == -1
    assert informative_score({"empty": 30, "full": 1}) == -1


def test_rank_key_orders():
    a = {"informative_score": 5, "n_carriers": 10, "size": 2, "locus_id": "a"}
    b = {"informative_score": -1, "n_carriers": 90, "size": 9, "locus_id": "b"}
    assert sorted([b, a], key=lambda r: rank_key(r, "informative"))[0] is a
    assert sorted([a, b], key=lambda r: rank_key(r, "strains"))[0] is b
    assert sorted([a, b], key=lambda r: rank_key(r, "size"))[0] is b


from island_locus import Columns, compute_locus  # noqa: E402


# ---- compute_locus end to end --------------------------------------------------------

class FakeMatrix:
    def __init__(self, calls):
        self.calls = calls

    def call(self, fam, strain):
        return self.calls.get((fam, strain), "absent")


def test_compute_locus_classifies_four_synthetic_strains():
    cols = Columns(left=("L1",), locus=("A", "B"), right=("R1",), left_avail=5, right_avail=5)
    pos = {}
    calls = {}

    def put(strain, contig, ranks):
        for fam, r in zip(["L1", "A", "B", "R1"], ranks):
            if r is not None:
                pos[(strain, fam)] = [(contig, r)]
                calls[(fam, strain)] = "present"

    put("FULL", "c1", [10, 11, 12, 13])
    put("EMPTY", "c1", [10, None, None, 11])
    put("PART", "c1", [10, 11, None, 13])
    put("FRAG", "c2", [0, 1, None, None])
    spans = {("FULL", "c1"): (0, 50), ("EMPTY", "c1"): (0, 50), ("PART", "c1"): (0, 50),
             ("FRAG", "c2"): (0, 2)}
    loc = Locus(root={"members": ["A", "B"], "locus_id": "FULL:c1:1-9"},
                variants=[{"n_strains": "2"}])
    res = compute_locus(loc, Placement("FULL", "c1", 11, 12, 11, 38, 51), "full", cols,
                        ["EMPTY", "FRAG", "FULL", "PART"], pos, spans, FakeMatrix(calls),
                        {"FULL": "sp1", "EMPTY": "sp2", "PART": "sp1"})
    assert res["counts"] == {"full": 1, "partial": 1, "empty": 1, "uninformative": 1}
    assert res["_row_class"] == {"EMPTY": "empty", "FRAG": "uninformative", "FULL": "full",
                                 "PART": "partial"}
    frag = [r for r in res["rows"] if r["strains"] == ["FRAG"]][0]
    assert frag["codes"] == "1155"
    assert res["n_carriers"] == 3
    assert res["counts_by_species"]["sp2"]["empty"] == 1
    assert res["breakpoints"] == [
        {"b": 1, "indel": {"sp2": 1}, "contig_break": 0},
        {"b": 2, "indel": {"sp1": 1}, "contig_break": 1},
        {"b": 3, "indel": {"sp1": 1, "sp2": 1}, "contig_break": 0},
    ]
    assert res["informative_score"] == -1
