import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_pair_classification import load_captain_families


def test_load_captain_families_default_id_sep(tmp_path):
    tblout = tmp_path / "captain.tblout"
    tblout.write_text(
        "s1|proteinA - - - - - - - - - - - - - - - -\n"
        "s2|proteinB - - - - - - - - - - - - - - - -\n"
    )
    member_to_rep = {"s1|proteinA": "famA", "s2|proteinB": "famB"}
    result = load_captain_families(str(tblout), member_to_rep)
    assert result == {"s1": {"famA"}, "s2": {"famB"}}


def test_load_captain_families_threads_custom_id_sep(tmp_path):
    # F6: pre-existing bug this consolidation fixes -- a study configured
    # with a non-default --pangenome_id_sep previously got zero captain
    # hits, silently, because load_captain_families hardcoded id_sep="|".
    tblout = tmp_path / "captain.tblout"
    tblout.write_text("s1__proteinA - - - - - - - - - - - - - - - -\n")
    member_to_rep = {"s1__proteinA": "famA"}
    assert load_captain_families(str(tblout), member_to_rep) == {}  # default "|" finds nothing
    assert load_captain_families(str(tblout), member_to_rep, id_sep="__") == {"s1": {"famA"}}


# --- Decision 1 (2026-10-08): representatives only for the linkage evidence -------------------
from pangenome_pair_classification import classify_pair, restrict_to_representatives  # noqa: E402


def _positions(strains, pair_close=True):
    """{strain: {family: [(contig, rank)]}}: famA and famB adjacent (linked) in each strain."""
    return {s: {"famA": [("c1", 1)], "famB": [("c1", 2 if pair_close else 500)]} for s in strains}


def test_restrict_to_representatives_keeps_only_listed_strains():
    gp = _positions(["R1", "R2", "N1", "N2", "N3"])
    out = restrict_to_representatives(gp, {"R1", "R2"})
    assert set(out) == {"R1", "R2"}


def test_clones_no_longer_supply_the_minimum_co_carrying_count():
    # 2 representatives + 3 near-identical strains all carry both families: 5 carriers, which
    # meets min_co_carrying = 5 when every strain counts but not when only the 2 representatives do.
    gp = _positions(["R1", "R2", "N1", "N2", "N3"])
    args = dict(family_a="famA", family_b="famB", permutation_p=0.001, clade_composition={"x": 2, "y": 2},
                captain_families={})
    all_strains = classify_pair(gene_position=gp, **args)
    reps_only = classify_pair(gene_position=restrict_to_representatives(gp, {"R1", "R2"}), **args)
    assert all_strains[0] != "insufficient_data"
    assert reps_only == ("insufficient_data", 0.0)


def test_linkage_fraction_is_computed_over_representatives():
    # Representatives are not linked; the three clones are linked. All-strain fraction is 3/5,
    # representatives-only fraction is 0.
    gp = {"R1": {"famA": [("c1", 1)], "famB": [("c1", 500)]}}
    gp |= {f"R{i}": {"famA": [("c1", 1)], "famB": [("c1", 500)]} for i in range(2, 6)}
    gp |= {f"N{i}": {"famA": [("c1", 1)], "famB": [("c1", 2)]} for i in range(1, 4)}
    reps = {f"R{i}" for i in range(1, 6)}
    _, frac_all = classify_pair("famA", "famB", 0.001, {"x": 1, "y": 1}, gp, {})
    _, frac_reps = classify_pair("famA", "famB", 0.001, {"x": 1, "y": 1}, restrict_to_representatives(gp, reps), {})
    assert abs(frac_all - 3 / 8) < 1e-9 and frac_reps == 0.0


def test_process_passes_the_inventory_and_the_evidence_setting():
    nf = (Path(__file__).parent.parent / "modules" / "pangenome" / "pair_classification.nf").read_text()
    assert "--evidence_strains ${params.pangenome_pair_class_evidence_strains}" in nf
    assert "--inventory ${strain_inventory}" in nf
    assert "stageAs: 'strain_inventory.tsv'" in nf


def test_nextflow_config_default_is_representatives():
    cfg = (Path(__file__).parent.parent / "nextflow.config").read_text()
    assert "pangenome_pair_class_evidence_strains   = 'representatives'" in cfg
