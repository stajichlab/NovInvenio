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
