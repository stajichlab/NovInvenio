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
