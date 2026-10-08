"""Uniform island column labels: the display strain's own gene IDs, not mmseqs representatives."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_labels import LABEL_REF, LABEL_REP, column_labels  # noqa: E402
from island_locus import locus_payload  # noqa: E402
from island_synteny import build_payload  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402
from pangenome_matrix import PresenceMatrix  # noqa: E402

# Representatives come from three different strains, as mmseqs picks them.
REPS = ["StrainZ|rep1", "StrainQ|rep2", "StrainK|rep3"]


def test_label_is_the_display_strains_gene_for_each_family():
    locs = [{"protein": "AFUA_1", "contig": "c1", "start": 1, "end": 9},
            {"protein": "AFUA_2", "contig": "c1", "start": 20, "end": 29}]
    labels, kinds = column_labels(REPS[:2], locs)
    assert labels == ["AFUA_1", "AFUA_2"]
    assert kinds == [LABEL_REF, LABEL_REF]


def test_family_without_a_gene_in_the_display_strain_keeps_its_representative_id():
    locs = [{"protein": "AFUA_1"}, None, {"rescued": True, "contig": "c1", "start": 5}]
    labels, kinds = column_labels(REPS, locs)
    assert labels == ["AFUA_1", "StrainQ|rep2", "StrainK|rep3"]
    assert kinds == [LABEL_REF, LABEL_REP, LABEL_REP]


def test_no_locations_at_all_gives_representative_labels():
    labels, kinds = column_labels(REPS, None)
    assert labels == REPS and kinds == [LABEL_REP] * 3


def test_island_payload_labels_follow_the_example_strain_and_keep_family_ids():
    fams = ["famA", "famB"]
    m = PresenceMatrix(families=fams, strains=["S1", "S2"])
    for f in fams:
        m.set_call(f, "S1", "present")
    row = {"locus_id": "S1:c1:100-400", "locus_contig": "c1", "locus_start": "100",
           "locus_end": "400", "example_strain": "S1", "island_size": "2",
           "n_strains": "2", "member_families": "famA,famB", "pfam_domains": "-"}
    positions = {("S1", "famA"): [("c1", 1)], ("S1", "famB"): [("c1", 2)]}
    locs = {("S1", "famA"): [("S1_gene_a", "c1", 100, 200)]}      # famB: no gene in S1
    isl = build_payload([row], m, positions, "demo", gene_locations=locs)["islands"][0]
    assert isl["families"] == fams                                 # stable key unchanged
    assert isl["family_labels"] == ["S1_gene_a", "famB"]
    assert isl["label_kinds"] == [LABEL_REF, LABEL_REP]


def test_island_payload_without_gene_locations_has_representative_labels():
    m = PresenceMatrix(families=["famA", "famB"], strains=["S1", "S2"])
    for f in ("famA", "famB"):
        m.set_call(f, "S1", "present")
    row = {"locus_id": "S1:c1:1-9", "locus_contig": "c1", "locus_start": "1", "locus_end": "9",
           "example_strain": "S1", "island_size": "2", "n_strains": "2",
           "member_families": "famA,famB", "pfam_domains": "-"}
    isl = build_payload([row], m, {("S1", "famA"): [("c1", 1)], ("S1", "famB"): [("c1", 2)]},
                        "demo")["islands"][0]
    assert isl["family_labels"] == ["famA", "famB"]
    assert isl["label_kinds"] == [LABEL_REP, LABEL_REP]


def test_locus_payload_labels_follow_the_exemplars_genes():
    result = {"families": ["famA", "famB"]}
    locs = [{"protein": "EX_1"}, None]
    out = locus_payload(result, "k", {}, ["unannotated"] * 2, ["", ""], "unannotated",
                        family_locations=locs)
    assert out["family_labels"] == ["EX_1", "famB"]
    assert out["label_kinds"] == [LABEL_REF, LABEL_REP]
    assert out["families"] == ["famA", "famB"]


def test_page_draws_labels_not_representative_ids_and_explains_italics():
    t = ISLAND_SYNTENY_TEMPLATE
    assert "function colLabel(" in t and "LABEL_FONT_REP" in t
    assert "ellipsize(hctx, colLabel(isl, i), LABEL_MAX_W)" in t        # synteny head
    assert "ellipsize(lhctx, colLabel(locus, i), LABEL_MAX_W)" in t     # locus head
    assert "ellipsize(hctx, fam, LABEL_MAX_W)" not in t                 # old rep-ID draw gone
    assert "an italic label is a family" in t                           # legend sentence, both views
