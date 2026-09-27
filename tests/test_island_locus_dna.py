"""DNA presence check in lib/island_locus.py (spec section 4b), on the
four synthetic strains of tests/test_island_locus.py's compute_locus test.
A strain's gene at rank r sits at bp r*1000+1 .. r*1000+800."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_locus import (  # noqa: E402
    Columns, Locus, Placement, apply_dna_calls, compute_locus, copy_bp, dna_checked_strains,
    dna_query, dna_target, gene_coverage, merge_intervals, rescue_hit_spans, row_class_dna,
)

FAMS = ["L1", "A", "B", "R1"]
SPECIES = {"FULL": "sp1", "EMPTY": "sp2", "PART": "sp1"}


class FakeMatrix:
    def __init__(self, calls):
        self.calls = calls

    def call(self, fam, strain):
        return self.calls.get((fam, strain), "absent")


def fixture():
    pos, calls, locs = {}, {}, {}

    def put(strain, contig, ranks):
        for fam, r in zip(FAMS, ranks):
            if r is not None:
                pos[(strain, fam)] = [(contig, r)]
                calls[(fam, strain)] = "present"
                locs[(strain, fam)] = [(f"{strain}_{fam}", contig, r * 1000 + 1, r * 1000 + 800)]

    put("FULL", "c1", [10, 11, 12, 13])
    put("EMPTY", "c1", [10, None, None, 11])
    put("PART", "c1", [10, 11, None, 13])
    put("FRAG", "c2", [0, 1, None, None])
    spans = {("FULL", "c1"): (0, 50), ("EMPTY", "c1"): (0, 50), ("PART", "c1"): (0, 50),
             ("FRAG", "c2"): (0, 2)}
    cols = Columns(left=("L1",), locus=("A", "B"), right=("R1",), left_avail=5, right_avail=5)
    loc = Locus(root={"members": ["A", "B"], "locus_id": "FULL:c1:1-9"},
                variants=[{"n_strains": "2"}])
    res = compute_locus(loc, Placement("FULL", "c1", 11, 12, 11, 38, 51), "full", cols,
                        ["EMPTY", "FRAG", "FULL", "PART"], pos, spans, FakeMatrix(calls), SPECIES)
    return res, pos, locs


def test_checked_strains_are_flank_intact_and_not_full():
    res, _, _ = fixture()
    assert dna_checked_strains(res) == ["EMPTY", "PART"]


def test_copy_bp_matches_copies_in_rank_and_start_order():
    pos = {("S", "F"): [("c1", 7), ("c1", 3), ("c2", 1)]}
    locs = {("S", "F"): [("p2", "c1", 900, 990), ("p1", "c1", 100, 190), ("p3", "c2", 5, 50)]}
    assert copy_bp(pos, locs, "S", "F", "c1", 3) == (100, 190)
    assert copy_bp(pos, locs, "S", "F", "c1", 7) == (900, 990)


def test_copy_bp_is_none_for_a_rescue_copy_or_a_count_mismatch():
    pos = {("S", "F"): [("c1", 3)], ("S", "G"): [("c1", 4), ("c1", 5)]}
    locs = {("S", "G"): [("p1", "c1", 10, 20)]}
    assert copy_bp(pos, locs, "S", "F", "c1", 3) is None
    assert copy_bp(pos, locs, "S", "G", "c1", 4) is None


def test_query_spans_the_exemplar_locus_genes():
    res, pos, locs = fixture()
    assert dna_query(res, [11, 12], pos, locs) == (
        "c1", 11001, 12800, [(1, 11001, 11800), (2, 12001, 12800)])


def test_query_leaves_out_a_column_without_a_gene_model_or_a_rescue_span():
    res, pos, locs = fixture()
    del locs[("FULL", "B")]
    assert dna_query(res, [11, 12], pos, locs)[3] == [(1, 11001, 11800)]
    assert dna_query(res, [11, 12], pos, locs, {("FULL", "B"): [("c9", 5, 90)]})[3] == [
        (1, 11001, 11800)]


def test_query_uses_the_rescue_hit_span_of_a_rescue_only_exemplar_column():
    # Locus 1M0:scaffold_390:15977-16664 (plan Task 19): the exemplar's
    # second locus gene is a TBLASTN rescue hit. Spec 4b: its span comes from rescue_positions.tsv
    # (Ruling R26), so the column gets a DNA call like any other.
    res, pos, locs = fixture()
    del locs[("FULL", "B")]
    spans = {("FULL", "B"): [("c1", 12101, 12900)]}
    assert copy_bp(pos, locs, "FULL", "B", "c1", 12, spans) == (12101, 12900)
    assert dna_query(res, [11, 12], pos, locs, spans) == (
        "c1", 11001, 12900, [(1, 11001, 11800), (2, 12101, 12900)])


def test_rescue_hit_spans_take_the_end_of_the_hsp_at_the_recorded_start():
    # rescue_positions.tsv keeps only min(sstart, send) of the best qualifying
    # HSP (bin/pangenome_extract_rescue_positions.py); the end is that HSP's
    # max(sstart, send) in the strain's own tblastn output (outfmt 6 std qcovs).
    rows = [("EX", "famR", "c1", 15568), ("EX", "famQ", "c1", 16013), ("EX", "famZ", "c1", 1)]
    lines = [
        # a higher-bitscore HSP of famR that does not start at 15568 is ignored
        "famR\tEX|c1\t85.4\t192\t28\t0\t28\t219\t16588\t16013\t1e-157\t300\t100\n",
        "famR\tEX|c1\t100.0\t132\t0\t0\t214\t345\t15963\t15568\t1e-157\t272\t100\n",
        "famR\tEX|c1\t99.0\t100\t0\t0\t214\t300\t15568\t15700\t1e-50\t150\t100\n",
        "famQ\tEX|c1\t94.9\t176\t9\t0\t20\t195\t16540\t16013\t1e-130\t310\t91\n",
        "famQ\tOTHER|c1\t99.0\t176\t9\t0\t20\t195\t16013\t16900\t1e-130\t400\t91\n",
        "# a comment\n", "short\tline\n",
    ]
    assert rescue_hit_spans(rows, lines) == {("EX", "famR"): [("c1", 15568, 15963)],
                                             ("EX", "famQ"): [("c1", 16013, 16540)]}


def test_target_spans_the_innermost_flank_genes_themselves():
    res, pos, locs = fixture()
    assert dna_target(res, "EMPTY", pos, locs) == ("c1", 10001, 11800)
    assert dna_target(res, "PART", pos, locs) == ("c1", 10001, 13800)


def test_target_keeps_locus_dna_that_lies_inside_a_flank_gene_model():
    # Guerrero_1 at scaffold_390 (plan Task 19): the strain's flank gene model
    # extends over the site where the exemplar has a locus gene. The target
    # must include that model, or the locus DNA is missed.
    res, pos, locs = fixture()
    locs[("EMPTY", "R1")] = [("EMPTY_R1", "c1", 10850, 12900)]
    contig, start, end = dna_target(res, "EMPTY", pos, locs)
    assert (contig, start, end) == ("c1", 10001, 12900)
    assert start <= 11001 and 12800 <= end  # the exemplar's locus genes' span fits inside


def test_no_target_without_an_annotated_in_place_flank_gene():
    res, pos, locs = fixture()
    del locs[("EMPTY", "R1")]
    assert dna_target(res, "EMPTY", pos, locs) is None


def test_merge_intervals_joins_overlaps_and_neighbours():
    assert merge_intervals([(5, 9), (1, 3), (4, 4), (20, 12)]) == [(1, 9), (12, 20)]


def test_coverage_merges_hsps_and_ignores_low_identity():
    hsps = [(1, 60, 99.0), (50, 100, 95.0), (101, 200, 80.0)]
    assert gene_coverage(hsps, (1, 100), 90) == 1.0
    assert gene_coverage(hsps, (1, 200), 90) == 0.5
    assert gene_coverage([], (1, 10), 90) == 0.0


def test_row_class_dna():
    assert row_class_dna("1661", 1, 2, True) == "model_difference"
    assert row_class_dna("1771", 1, 2, True) == "empty"
    assert row_class_dna("1171", 1, 2, True) == "partial"
    assert row_class_dna("1671", 1, 2, True) == "partial"
    assert row_class_dna("1601", 1, 2, True) == "partial"
    assert row_class_dna("1111", 1, 2, True) == "full"
    assert row_class_dna("1661", 1, 2, False) == "uninformative"


def test_dna_present_turns_an_empty_site_into_a_model_difference():
    res, _, _ = fixture()
    apply_dna_calls(res, {"EMPTY": {1: "present", 2: "present"}, "PART": {2: "absent"}}, SPECIES)
    assert res["counts"] == {"full": 1, "partial": 1, "empty": 0, "model_difference": 1,
                             "uninformative": 1}
    assert res["_row_class"]["EMPTY"] == "model_difference"
    assert [(r["row_class"], r["codes"]) for r in res["rows"]] == [
        ("full", "1111"), ("partial", "1171"), ("model_difference", "1661"),
        ("uninformative", "1155")]
    assert res["dna"] == {"checked": 2, "unchecked": 0, "empty_confirmed": 0,
                          "empty_to_model_difference": 1}
    assert res["counts_by_species"]["sp2"]["model_difference"] == 1


def test_breakpoints_count_only_in_place_to_dna_absent():
    res, _, _ = fixture()
    apply_dna_calls(res, {"EMPTY": {1: "present", 2: "present"}, "PART": {2: "absent"}}, SPECIES)
    assert res["breakpoints"] == [
        {"b": 2, "indel": {"sp1": 1}, "contig_break": 1},
        {"b": 3, "indel": {"sp1": 1}, "contig_break": 0},
    ]


def test_dna_absent_confirms_the_empty_site():
    res, _, _ = fixture()
    apply_dna_calls(res, {"EMPTY": {1: "absent", 2: "absent"}}, SPECIES)
    assert res["_row_class"]["EMPTY"] == "empty"
    assert [r["codes"] for r in res["rows"] if r["row_class"] == "empty"] == ["1771"]
    assert res["dna"]["empty_confirmed"] == 1 and res["dna"]["unchecked"] == 1


def test_unchecked_strains_keep_their_states_and_do_not_count_as_confirmed():
    res, _, _ = fixture()
    apply_dna_calls(res, {"EMPTY": {1: "unchecked", 2: "unchecked"}}, SPECIES)
    assert res["counts"]["empty"] == 1 and res["counts"]["model_difference"] == 0
    assert res["dna"] == {"checked": 0, "unchecked": 2, "empty_confirmed": 0,
                          "empty_to_model_difference": 0}


def test_informative_score_uses_dna_confirmed_empty_sites_only():
    res, _, _ = fixture()
    for i in range(12):
        res["_cells"][f"E{i}"] = res["_cells"]["EMPTY"]
        res["_pairs"][f"E{i}"] = res["_pairs"]["EMPTY"]
        res["_row_class"][f"E{i}"] = "empty"
    for i in range(3):
        res["_cells"][f"F{i}"] = res["_cells"]["FULL"]
        res["_row_class"][f"F{i}"] = "full"
    absent = {f"E{i}": {1: "absent", 2: "absent"} for i in range(12)}
    apply_dna_calls(res, absent, SPECIES)
    assert res["informative_score"] == 4
    apply_dna_calls(res, {}, SPECIES)
    assert res["informative_score"] == -1
