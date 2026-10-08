"""Clinker strain selection and per-strain regions (spec section 8)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import (  # noqa: E402
    Columns, Locus, Placement, apply_dna_calls, compute_locus, select_clinker_strains,
    strain_region,
)

FAMS = ["L1", "L2", "A", "B", "R1", "R2"]


class FakeMatrix:
    def __init__(self, calls):
        self.calls = calls

    def call(self, fam, strain):
        return self.calls.get((fam, strain), "absent")


def build(strains):
    """strains: {name: (contig, [rank or None per column], species)}."""
    pos, calls, spans, species = {}, {}, {}, {}
    for s, (contig, ranks, sp) in strains.items():
        species[s] = sp
        for fam, r in zip(FAMS, ranks):
            if r is not None:
                pos[(s, fam)] = [(contig, r)]
                calls[(fam, s)] = "present"
        spans[(s, contig)] = (0, 100)
    cols = Columns(left=("L1", "L2"), locus=("A", "B"), right=("R1", "R2"),
                   left_avail=20, right_avail=20)
    loc = Locus(root={"members": ["A", "B"], "locus_id": "X"}, variants=[{"n_strains": "3"}])
    ex = sorted(strains)[0]
    res = compute_locus(loc, Placement(ex, strains[ex][0], 22, 23, 20, 77, 101), "full", cols,
                        sorted(strains), pos, spans, FakeMatrix(calls), species)
    return res, spans, species


FULL = [20, 21, 22, 23, 24, 25]
EMPTY = [20, 21, None, None, 22, 23]
PART = [20, 21, 22, None, 24, 25]


def test_region_spans_in_place_copies_plus_flank_genes():
    res, spans, _ = build({"A1": ("c1", FULL, "sp")})
    assert strain_region(res, "A1", spans, flank=5) == {
        "contig": "c1", "rank_lo": 15, "rank_hi": 30,
        "anchors": [(20, "L1"), (21, "L2"), (22, "A"), (23, "B"), (24, "R1"), (25, "R2")]}


def test_empty_site_region_is_the_two_flank_blocks():
    res, spans, _ = build({"A1": ("c1", FULL, "sp"), "E1": ("c1", EMPTY, "sp")})
    reg = strain_region(res, "E1", spans, flank=5)
    assert (reg["rank_lo"], reg["rank_hi"]) == (15, 28)
    assert [f for _, f in reg["anchors"]] == ["L1", "L2", "R1", "R2"]


def test_region_is_clipped_at_the_contig_end():
    res, spans, _ = build({"A1": ("c1", [1, 2, 3, 4, 5, 6], "sp")})
    assert strain_region(res, "A1", spans, flank=5)["rank_lo"] == 0


def test_unanchored_strain_has_no_region():
    res, spans, _ = build({"A1": ("c1", FULL, "sp"), "U1": ("c2", [None, None, 5, None, None, None], "sp")})
    assert strain_region(res, "U1", spans, flank=5) is None


def test_selection_order_exemplar_then_class_species_then_fill():
    strains = {"A1": ("c1", FULL, "sp1"), "B1": ("c1", FULL, "sp2"), "C1": ("c1", PART, "sp1"),
               "D1": ("c1", EMPTY, "sp2"), "E1": ("c1", FULL, "sp1"), "F1": ("c1", PART, "sp1"),
               "U1": ("c2", [None, None, 5, None, None, None], "sp1")}
    res, spans, species = build(strains)
    picks = select_clinker_strains(res, {"E1": 900, "A1": 100, "C1": 50, "F1": 60}, species,
                                   spans, flank=5)
    assert [(p["strain"], p["reason"]) for p in picks] == [
        ("A1", "exemplar"), ("E1", "best_full"), ("B1", "best_full"), ("F1", "best_partial"),
        ("D1", "best_empty"), ("C1", "fill")]
    assert "U1" not in [p["strain"] for p in picks]


def test_selection_is_capped():
    strains = {f"S{i:02d}": ("c1", FULL, "sp") for i in range(20)}
    res, spans, species = build(strains)
    assert len(select_clinker_strains(res, {}, species, spans, flank=5, max_strains=12)) == 12


def test_each_pick_carries_its_region_and_class():
    res, spans, species = build({"A1": ("c1", FULL, "sp")})
    (pick,) = select_clinker_strains(res, {}, species, spans, flank=5)
    assert pick["row_class"] == "full" and pick["contig"] == "c1" and pick["rank_hi"] == 30


def test_a_model_difference_strain_is_chosen_after_the_empty_sites():
    # Spec section 8 item 2 with section 4b: model difference is a row class.
    res, spans, species = build({"A1": ("c1", FULL, "sp"), "D1": ("c1", EMPTY, "sp"),
                                 "D2": ("c1", EMPTY, "sp")})
    apply_dna_calls(res, {"D1": {2: "present", 3: "present"}}, species, model_diff_min_frac=0)
    picks = select_clinker_strains(res, {}, species, spans, flank=5)
    assert [(p["strain"], p["reason"], p["row_class"]) for p in picks] == [
        ("A1", "exemplar", "full"), ("D2", "best_empty", "empty"),
        ("D1", "best_model_difference", "model_difference")]
