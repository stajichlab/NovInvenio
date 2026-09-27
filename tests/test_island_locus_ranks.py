"""Multi-rank scores in lib/island_locus.py (ranks-brief.md, spec sections
4b/6, changed 2026-09-27): five rankings (A whole_annot, B whole_dna,
C presence [default, alias "informative"], D species, E within) replace the
single informative-polymorphism ranking."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import (  # noqa: E402
    ABSENT, DNA_ABSENT, locus_ranks, locus_species_stats, rank_sort_key, score_presence,
    score_species, score_whole_annot, score_whole_dna, score_within,
)


# ---- locus_species_stats: carriers/losses partition --------------------------

def test_full_and_model_difference_are_carriers_empty_is_a_loss():
    per_strain = {"A": ("full", "11"), "B": ("model_difference", "16"), "C": ("empty", "77")}
    stats = locus_species_stats(per_strain, {}, 0, 2, ABSENT)
    assert stats["carriers"] == 2 and stats["losses"] == 1


def test_partial_with_a_missing_locus_column_is_a_loss():
    # n_left=0, n_locus=2: locus block is codes[0:2].
    per_strain = {"A": ("partial", "10"), "B": ("partial", "12")}
    stats = locus_species_stats(per_strain, {}, 0, 2, ABSENT)
    assert stats["carriers"] == 1 and stats["losses"] == 1  # "12" has no "0"; "10" does


def test_dna_on_uses_dna_absent_as_the_missing_code_not_absent():
    # With the DNA check on, a literal "0" in the locus block is unchecked,
    # not "missing" -- only DNA_ABSENT ("7") counts as missing.
    per_strain = {"A": ("partial", "10"), "B": ("partial", "17")}
    stats = locus_species_stats(per_strain, {}, 0, 2, DNA_ABSENT)
    assert stats["carriers"] == 1 and stats["losses"] == 1  # "10" carries; "17" is a loss


def test_uninformative_strains_are_excluded():
    per_strain = {"A": ("uninformative", "55"), "B": ("full", "11")}
    stats = locus_species_stats(per_strain, {}, 0, 2, ABSENT)
    assert stats["carriers"] == 1 and stats["losses"] == 0
    assert "" not in stats["by_species"] or stats["by_species"][""]["car"] + \
        stats["by_species"][""]["loss"] == 1


def test_per_species_counts():
    per_strain = {"A": ("full", "11"), "B": ("empty", "00"), "C": ("full", "11")}
    species_of = {"A": "sp1", "B": "sp1", "C": "sp2"}
    stats = locus_species_stats(per_strain, species_of, 0, 2, ABSENT)
    assert stats["by_species"] == {"sp1": {"car": 1, "loss": 1}, "sp2": {"car": 1, "loss": 0}}


# ---- score functions ----------------------------------------------------------

def test_score_whole_annot_threshold():
    assert score_whole_annot(10, 2) == 2
    assert score_whole_annot(9, 50) == -1
    assert score_whole_annot(30, 1) == -1


def test_score_whole_dna_counts_model_difference_with_full():
    assert score_whole_dna(10, 1, 1) == 2  # full=1 + model_difference=1 = 2
    assert score_whole_dna(10, 0, 1) == -1  # full + model_difference = 1 < 2
    assert score_whole_dna(9, 5, 5) == -1


def test_score_presence_threshold():
    assert score_presence(losses=10, carriers=2) == 2
    assert score_presence(losses=9, carriers=50) == -1
    assert score_presence(losses=30, carriers=1) == -1


def test_score_species_needs_two_species_at_min_n_and_the_fixed_diff():
    by_species = {"sp1": {"car": 0, "loss": 20}, "sp2": {"car": 20, "loss": 0}}  # n=20 each
    assert score_species(by_species, n_species_total=2) == 1.0  # 1.0 - 0.0
    by_species2 = {"sp1": {"car": 1, "loss": 19}, "sp2": {"car": 19, "loss": 1}}  # diff 0.9 < 0.95
    assert score_species(by_species2, n_species_total=2) == -1


def test_score_species_is_always_minus_one_with_fewer_than_two_species_total():
    by_species = {"sp1": {"car": 0, "loss": 20}, "sp2": {"car": 20, "loss": 0}}
    assert score_species(by_species, n_species_total=1) == -1


def test_score_species_needs_two_qualifying_species_even_with_more_in_samplesheet():
    by_species = {"sp1": {"car": 0, "loss": 20}, "sp2": {"car": 1, "loss": 1}}  # sp2 n=2 < 10
    assert score_species(by_species, n_species_total=2) == -1


def test_score_within_picks_the_species_with_the_larger_minor_count():
    by_species = {
        "sp1": {"car": 15, "loss": 5},   # n=20, f=0.25, min=5
        "sp2": {"car": 30, "loss": 10},  # n=40, f=0.25, min=10
        "sp3": {"car": 5, "loss": 95},   # n=100, f=0.95 (boundary, included), min=5
        "sp4": {"car": 25, "loss": 5},   # n=30, f<0.05 boundary excluded (0.1667 not <0.05, included) min=5
    }
    score, sp = score_within(by_species)
    assert (score, sp) == (10, "sp2")


def test_score_within_excludes_species_below_poly_min_strains_or_outside_the_frac_window():
    by_species = {"small": {"car": 14, "loss": 1}, "extreme": {"car": 20, "loss": 0}}
    assert score_within(by_species) == (-1, None)  # small: n=15 < 20; extreme: frac 0 outside window


def test_score_within_ties_go_to_the_alphabetically_first_species():
    by_species = {"zzz": {"car": 10, "loss": 10}, "aaa": {"car": 10, "loss": 10}}
    assert score_within(by_species) == (10, "aaa")


# ---- locus_ranks + rank_sort_key -----------------------------------------------

def test_locus_ranks_one_species_samplesheet_gives_species_minus_one():
    result = {
        "n_left": 0, "n_locus": 2,
        "counts": {"full": 2, "partial": 0, "empty": 10, "uninformative": 0},
        "_per_strain": dict(
            [(f"F{i}", ("full", "11")) for i in range(2)] +
            [(f"E{i}", ("empty", "00")) for i in range(10)]),
    }
    species_of = {**{f"F{i}": "sp1" for i in range(2)}, **{f"E{i}": "sp1" for i in range(10)}}
    ranks = locus_ranks(result, species_of, dna_on=False, n_species_total=1)
    assert ranks["species"] == -1
    assert ranks["presence"] == 2  # min(losses=10, carriers=2)
    assert ranks["whole_annot"] == 2
    assert ranks["carriers"] == 2 and ranks["losses"] == 10


def test_rank_sort_key_orders_score_desc_then_strains_desc_then_locus_id():
    a = {"presence": 5, "carriers": 10, "losses": 0}
    b = {"presence": 5, "carriers": 90, "losses": 0}
    c = {"presence": -1, "carriers": 1, "losses": 0}
    keyed = sorted([("b", b), ("a", a), ("c", c)],
                  key=lambda kv: rank_sort_key(kv[1], "presence", kv[0]))
    assert [k for k, _ in keyed] == ["b", "a", "c"]
