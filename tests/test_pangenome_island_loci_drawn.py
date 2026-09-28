"""choose_drawn() in bin/pangenome_island_loci.py: the union drawn-set rule
(ranks-brief.md "Drawn set") -- top --per_rank loci under whole_annot,
whole_dna, species, within (score >= 0 only), then fill with loci in
presence (C) order, including -1, until --top_loci; deduplicated; the
result is always in C order."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from pangenome_island_loci import choose_drawn  # noqa: E402


def r(locus_id):
    return {"locus_id": locus_id}


def rk(presence=-1, whole_annot=-1, whole_dna=-1, species=-1, within=-1, carriers=0, losses=0):
    return {"presence": presence, "whole_annot": whole_annot, "whole_dna": whole_dna,
            "species": species, "within": within, "within_species": None,
            "carriers": carriers, "losses": losses}


def test_union_of_per_rank_winners_plus_c_fill():
    # L1 tops presence and whole_annot; L2 only tops "within"; L3 only tops
    # "species"; L4 is a C-order filler with no other rank informative.
    results = [r("L1"), r("L2"), r("L3"), r("L4")]
    ranks = {
        "L1": rk(presence=10, whole_annot=5, carriers=12, losses=10),
        "L2": rk(presence=1, within=8, carriers=5, losses=5),
        "L3": rk(presence=2, species=0.97, carriers=4, losses=4),
        "L4": rk(presence=0, carriers=1, losses=1),
    }
    drawn = choose_drawn(results, ranks, top_loci=4, per_rank=1)
    assert [x["locus_id"] for x in drawn] == ["L1", "L3", "L2", "L4"]  # presence (C) order


def test_per_rank_caps_how_many_join_from_each_non_presence_rank():
    results = [r(f"L{i}") for i in range(1, 6)]
    ranks = {
        "L1": rk(presence=1, whole_annot=9, carriers=1, losses=1),
        "L2": rk(presence=2, whole_annot=8, carriers=1, losses=1),
        "L3": rk(presence=3, whole_annot=7, carriers=1, losses=1),  # 3rd-best whole_annot, excluded
        "L4": rk(presence=4, carriers=1, losses=1),
        "L5": rk(presence=5, carriers=1, losses=1),
    }
    drawn = choose_drawn(results, ranks, top_loci=2, per_rank=2)
    ids = [x["locus_id"] for x in drawn]
    assert set(ids) == {"L1", "L2"}  # top 2 whole_annot winners already fill top_loci=2


def test_negative_scores_never_join_a_per_rank_top_but_can_fill_via_c():
    results = [r("L1"), r("L2")]
    ranks = {"L1": rk(presence=-1, carriers=1, losses=1), "L2": rk(presence=-1, carriers=2, losses=2)}
    drawn = choose_drawn(results, ranks, top_loci=5, per_rank=20)
    # both -1 under every rank; fill takes all of them in C order (ties: strain count desc)
    assert [x["locus_id"] for x in drawn] == ["L2", "L1"]


def test_dedup_a_locus_topping_two_ranks_is_drawn_once():
    results = [r("L1"), r("L2")]
    ranks = {
        "L1": rk(presence=1, whole_annot=9, within=9, carriers=1, losses=1),
        "L2": rk(presence=2, carriers=1, losses=1),
    }
    drawn = choose_drawn(results, ranks, top_loci=2, per_rank=1)
    assert [x["locus_id"] for x in drawn] == ["L2", "L1"]


def test_result_size_can_exceed_top_loci_when_per_rank_unions_already_do():
    results = [r(f"L{i}") for i in range(1, 5)]
    ranks = {
        "L1": rk(presence=1, whole_annot=9, carriers=1, losses=1),
        "L2": rk(presence=2, whole_dna=9, carriers=1, losses=1),
        "L3": rk(presence=3, species=0.99, carriers=1, losses=1),
        "L4": rk(presence=4, within=9, carriers=1, losses=1),
    }
    drawn = choose_drawn(results, ranks, top_loci=1, per_rank=1)
    assert len(drawn) == 4  # each locus alone tops one of the four other ranks
