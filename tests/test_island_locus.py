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
