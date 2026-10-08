"""Preferred reference strains (decision 5, 2026-10-08): the display strain of an island or locus is
the first listed strain that carries it; otherwise the existing rule applies."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import Placement, choose_exemplar  # noqa: E402
from pangenome_report_tables import add_island_locus  # noqa: E402


def _island(example="OTHER"):
    return {"n_strains": "3", "example_strain": example, "island_size": "2",
            "member_families": "OTHER|a,OTHER|b", "n_supporting_pairs": "1", "classifications": "x"}


def _member_to_rep():
    # member -> family representative; both strains' genes belong to the two families
    m = {}
    for s in ("OTHER", "REF", "REF2"):
        m[f"{s}|{s}_a"] = "OTHER|a"
        m[f"{s}|{s}_b"] = "OTHER|b"
    return m


def _positions(strains, missing=()):
    gp = {}
    for s in strains:
        for g, (st, en) in (("a", (100, 200)), ("b", (300, 400))):
            if (s, g) in missing:
                continue
            gp[(s, f"{s}_{g}")] = {"contig": f"c_{s}", "start": st, "end": en}
    return gp


def test_no_reference_list_keeps_the_example_strain():
    rows = add_island_locus([_island()], _member_to_rep(), _positions(["OTHER", "REF"]))
    assert rows[0]["example_strain"] == "OTHER"
    assert rows[0]["locus_id"].startswith("OTHER:")


def test_first_listed_strain_that_carries_the_island_becomes_the_display_strain():
    rows = add_island_locus([_island()], _member_to_rep(), _positions(["OTHER", "REF"]),
                            reference_strains=["REF"])
    assert rows[0]["example_strain"] == "REF"
    assert rows[0]["locus_id"] == "REF:c_REF:100-400"
    assert rows[0]["first_strain"] == "OTHER"                 # the strain BUILD_ISLANDS found first


def test_a_reference_that_lacks_a_family_is_skipped_for_the_next_one():
    gp = _positions(["OTHER", "REF", "REF2"], missing={("REF", "b")})
    rows = add_island_locus([_island()], _member_to_rep(), gp, reference_strains=["REF", "REF2"])
    assert rows[0]["example_strain"] == "REF2"


def test_if_no_listed_strain_carries_it_the_existing_strain_is_kept():
    gp = _positions(["OTHER", "REF"], missing={("REF", "a")})
    rows = add_island_locus([_island()], _member_to_rep(), gp, reference_strains=["REF"])
    assert rows[0]["example_strain"] == "OTHER"


def _place(strain, left, right):
    return Placement(strain, "c1", left, left + 1, left, right, 20)


def test_locus_exemplar_prefers_the_listed_strain_within_a_tier():
    picks = [_place("BEST_N50", 6, 6), _place("REF", 6, 6)]
    n50 = {"BEST_N50": 9_000_000, "REF": 1_000}
    assert choose_exemplar(picks, n50)[0].strain == "BEST_N50"
    assert choose_exemplar(picks, n50, preferred=["REF"])[0].strain == "REF"


def test_locus_exemplar_tier_still_comes_before_the_preference():
    picks = [_place("REF", 1, 1), _place("OTHER", 6, 6)]       # only OTHER has full flanks
    assert choose_exemplar(picks, {}, preferred=["REF"]) == (picks[1], "full")


def test_the_parameter_is_wired_through_nextflow():
    root = Path(__file__).parent.parent
    cfg = (root / "nextflow.config").read_text()
    assert "pangenome_island_reference_strains" in cfg
    for mod in ("report.nf", "island_loci.nf", "island_dna_check.nf"):
        nf = (root / "modules" / "pangenome" / mod).read_text()
        assert "--island_reference_strains" in nf, mod
