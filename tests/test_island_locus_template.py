"""Locus view JS in island_synteny.html: pure helpers run under node, plus
structural checks. Same brace-matching extraction as
tests/test_island_synteny.py."""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from island_locus_template import LOCUS_VIEW_CSS, LOCUS_VIEW_HTML, LOCUS_VIEW_JS  # noqa: E402
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402


def _extract_function(src: str, name: str) -> str:
    start = src.index("function " + name + "(")
    depth = 0
    for i in range(src.index("{", start), len(src)):
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                return src[start:i + 1]
    raise AssertionError(f"unterminated function {name}()")


def run_node(names, cases_js):
    node = shutil.which("node")
    if not node:
        pytest.skip("node not available")
    funcs = "\n".join(_extract_function(ISLAND_SYNTENY_TEMPLATE, n) for n in names)
    proc = subprocess.run([node, "--input-type=commonjs", "-e", funcs + "\n" + cases_js],
                          capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout.strip().splitlines()[-1])


def test_page_embeds_the_three_fragments():
    assert LOCUS_VIEW_CSS in ISLAND_SYNTENY_TEMPLATE
    assert LOCUS_VIEW_HTML in ISLAND_SYNTENY_TEMPLATE
    assert LOCUS_VIEW_JS in ISLAND_SYNTENY_TEMPLATE
    assert ISLAND_SYNTENY_TEMPLATE.index(LOCUS_VIEW_JS) < ISLAND_SYNTENY_TEMPLATE.index("// ---- init")


def test_island_view_is_wrapped_so_it_can_be_hidden():
    page = ISLAND_SYNTENY_TEMPLATE
    assert page.index('<div id="island-view">') < page.index('id="isv-body"')


def test_locus_fragments_use_no_innerhtml_and_no_hex_colours():
    assert "innerHTML" not in LOCUS_VIEW_JS
    for frag in (LOCUS_VIEW_CSS, LOCUS_VIEW_HTML, LOCUS_VIEW_JS):
        assert not re.findall(r"(?<![\w-])#[0-9a-fA-F]{3,8}(?![\w-])", frag)


def test_state_styles_encode_evidence_by_token():
    out = run_node(["locusStateStyle"], """
      console.log(JSON.stringify(["0","1","2","3","4","5"].map(locusStateStyle)));""")
    assert [o["token"] for o in out] == ["--grid", "--series-1", "--series-1", "--series-2",
                                         "--series-2", "--warn"]
    assert [o["hatch"] for o in out] == [False, False, False, True, True, True]
    assert out[2]["alpha"] < 1 and out[1]["alpha"] == 1


def test_rows_group_by_class_then_species_then_count():
    rows = [
        {"row_class": "empty", "codes": "10", "count": 5, "strains": ["e1"]},
        {"row_class": "full", "codes": "11", "count": 1, "strains": ["b1"]},
        {"row_class": "full", "codes": "12", "count": 9, "strains": ["a1"]},
        {"row_class": "uninformative", "codes": "00", "count": 50, "strains": ["u1"]},
        {"row_class": "full", "codes": "13", "count": 3, "strains": ["x1"]},
    ]
    species = {"a1": "sp2", "b1": "sp1", "e1": "sp1", "u1": "sp1"}
    out = run_node(["speciesCounts", "haplotypeSpecies", "locusSortedRows"], f"""
      var rows = {json.dumps(rows)};
      console.log(JSON.stringify(locusSortedRows(rows, {json.dumps(species)}).map(function (r) {{ return r.strains[0]; }})));""")
    assert out == ["b1", "a1", "x1", "e1", "u1"]


def test_sidebar_orders():
    loci = [{"informative_score": 5, "n_carriers": 9, "size": 2, "locus_id": "b"},
            {"informative_score": 40, "n_carriers": 3, "size": 7, "locus_id": "a"},
            {"informative_score": -1, "n_carriers": 99, "size": 1, "locus_id": "c"}]
    out = run_node(["locusSidebarOrder"], f"""
      var L = {json.dumps(loci)};
      console.log(JSON.stringify(["informative","strains","size","name"].map(function (k) {{
        return locusSidebarOrder(L, [0, 1, 2], k); }})));""")
    assert out == [[1, 0, 2], [2, 0, 1], [1, 0, 2], [1, 0, 2]]


def test_breakpoint_bars_are_one_per_species_then_contig_break():
    out = run_node(["breakpointBars"], """
      console.log(JSON.stringify(breakpointBars({b: 2, indel: {sp2: 4}, contig_break: 3}, ["sp1", "sp2"])));""")
    assert out == [{"kind": "indel", "species": "sp1", "n": 0},
                   {"kind": "indel", "species": "sp2", "n": 4},
                   {"kind": "contig_break", "species": "", "n": 3}]


def test_cell_reason_names_the_neighbour_and_rank():
    locus = {"families": ["F1", "A", "B"]}
    row = {"codes": "112", "count": 3, "strains": ["S1", "S2", "S3"],
           "rep": {"c": ["scaffold_9"], "d": [[0, 40, 1, 1], [0, 41, 0, -1], [0, 90]]}}
    out = run_node(["locusStateLabel", "cellReasonLines"], f"""
      var locus = {json.dumps(locus)}, row = {json.dumps(row)};
      console.log(JSON.stringify([cellReasonLines(locus, row, 1, 10), cellReasonLines(locus, row, 2, 10)]));""")
    assert out[0] == ["State: in place", "In S1 (first of 3 strains): scaffold_9, gene rank 41",
                      "Neighbour F1 at rank 40 (-1)"]
    assert out[1][2] == "No other column within 10 genes on this contig"


def test_cell_reason_for_contig_break_and_missing_detail():
    out = run_node(["locusStateLabel", "cellReasonLines"], """
      var locus = {families: ["F1"]};
      var r1 = {codes: "5", count: 1, strains: ["S1"], rep: {c: [], d: [null]}};
      var r2 = {codes: "0", count: 1, strains: ["S1"], rep: null};
      console.log(JSON.stringify([cellReasonLines(locus, r1, 0, 10), cellReasonLines(locus, r2, 0, 10)]));""")
    assert "not evidence" in out[0][1]
    assert out[1][1] == "Position detail not stored for this row"


def test_tier_text():
    out = run_node(["locusTierText"], """
      var p = {flank_min: 3};
      console.log(JSON.stringify([locusTierText("full", p), locusTierText("short_flanks", p), locusTierText("contig_end", p)]));""")
    assert out == ["", "short flanks (3 genes per side)", "exemplar at contig end"]


# ---- DNA presence check (spec section 4b) ----

def test_dna_states_are_hatched_grey_and_plain_absent():
    out = run_node(["locusStateStyle", "locusStateLabel"], """
      console.log(JSON.stringify([locusStateStyle("6"), locusStateStyle("7"),
        locusStateLabel("6"), locusStateLabel("7")]));""")
    assert out[0] == {"token": "--text-secondary", "alpha": 0.45, "hatch": True}
    assert out[1] == {"token": "--grid", "alpha": 1, "hatch": False}
    assert out[2].startswith("absent, DNA present") and out[3] == "absent, DNA absent"


def test_model_difference_rows_sort_after_empty_sites():
    rows = [{"row_class": "uninformative", "codes": "0", "count": 1, "strains": ["u"]},
            {"row_class": "model_difference", "codes": "6", "count": 1, "strains": ["m"]},
            {"row_class": "empty", "codes": "7", "count": 1, "strains": ["e"]}]
    out = run_node(["speciesCounts", "haplotypeSpecies", "locusSortedRows", "locusClassLabel"], f"""
      var r = locusSortedRows({json.dumps(rows)}, {{}});
      console.log(JSON.stringify([r.map(function (x) {{ return x.strains[0]; }}),
        locusClassLabel("model_difference")]));""")
    assert out == [["e", "m", "u"], "model difference"]


def test_sidebar_shows_the_model_difference_count_only_with_the_check():
    base = {"size": 2, "n_variants": 1,
            "counts": {"full": 3, "partial": 1, "empty": 0, "uninformative": 2}}
    dna = dict(base, counts=dict(base["counts"], model_difference=9))
    out = run_node(["locusSidebarStats"], f"""
      console.log(JSON.stringify([locusSidebarStats({json.dumps(base)}),
        locusSidebarStats({json.dumps(dna)})]));""")
    assert "model difference" not in out[0]
    assert "model difference 9" in out[1]


def test_note_says_whether_empty_site_is_dna_confirmed():
    out = run_node(["locusDnaNote"], """
      console.log(JSON.stringify([locusDnaNote({dna_check: true, dna_min_id: 90, dna_min_cov: 80}),
        locusDnaNote({dna_check: false}), locusDnaNote({})]));""")
    assert "DNA-confirmed" in out[0] and ">= 90% identity" in out[0]
    assert "not DNA-confirmed" in out[1] and "not DNA-confirmed" in out[2]


def test_cell_reason_explains_the_dna_state():
    out = run_node(["locusStateLabel", "cellReasonLines"], """
      var locus = {families: ["F1", "A"]};
      var row = {codes: "16", count: 1, strains: ["S1"], rep: {c: [], d: [null, null]}};
      var row7 = {codes: "17", count: 1, strains: ["S1"], rep: {c: [], d: [null, null]}};
      console.log(JSON.stringify([cellReasonLines(locus, row, 1, 10), cellReasonLines(locus, row7, 1, 10)]));""")
    assert "not a deletion" in out[0][-1]
    assert out[1][-1] == "The site lacks the exemplar gene's DNA in S1 (blastn)"


# F1: counts.empty includes unchecked strains once the DNA check has run
# (Ruling R22 only credits DNA-confirmed empty sites), so the sidebar and the
# main note must switch to locus.dna.empty_confirmed and surface
# locus.dna.unchecked separately instead of quietly folding it into "empty".
def test_sidebar_uses_dna_confirmed_empty_count_and_flags_unchecked():
    locus = {"size": 2, "n_variants": 1,
             "counts": {"full": 3, "partial": 1, "empty": 5, "model_difference": 9,
                        "uninformative": 2},
             "dna": {"checked": 8, "unchecked": 6, "empty_confirmed": 2,
                     "empty_to_model_difference": 9}}
    out = run_node(["locusSidebarStats"], f"""
      var locus = {json.dumps(locus)};
      console.log(JSON.stringify([locusSidebarStats(locus, {{dna_check: true}}),
        locusSidebarStats(locus, {{dna_check: false}}), locusSidebarStats(locus)]));""")
    assert "empty 2" in out[0] and "not checked 6" in out[0]
    assert "empty 5" in out[1] and "not checked" not in out[1]
    assert "empty 5" in out[2] and "not checked" not in out[2]


def test_note_reports_unchecked_strains_separately_when_dna_check_ran():
    out = run_node(["locusDnaNote"], """
      console.log(JSON.stringify([
        locusDnaNote({dna_check: true, dna_min_id: 90, dna_min_cov: 80},
                     {empty_confirmed: 7, unchecked: 3}),
        locusDnaNote({dna_check: true, dna_min_id: 90, dna_min_cov: 80},
                     {empty_confirmed: 7, unchecked: 0}),
        locusDnaNote({dna_check: true, dna_min_id: 90, dna_min_cov: 80})]));""")
    assert "DNA-confirmed for 7 of 10 checked strains" in out[0] and "3 had no DNA call" in out[0]
    assert "DNA-confirmed for" not in out[1] and "DNA-confirmed:" in out[1]
    assert "DNA-confirmed for" not in out[2] and "DNA-confirmed:" in out[2]


# F4: the legend was missing code 4 (rescue, elsewhere) and, with the DNA
# check on, code 0's default swatch (--grid, no hatch) is visually identical
# to code 7's -- only the label told them apart, and "absent" alone reads as
# if it were itself a DNA-confirmed call.
def test_legend_adds_rescue_elsewhere_code():
    out = run_node(["legendCodes"], """
      console.log(JSON.stringify([legendCodes(false), legendCodes(true)]));""")
    assert out[0] == ["1", "2", "3", "4", "0", "5"]
    assert out[1] == ["1", "2", "3", "4", "0", "5", "6", "7"]


def test_legend_label_distinguishes_not_checked_from_dna_absent():
    out = run_node(["locusStateLabel", "legendLabel"], """
      console.log(JSON.stringify([legendLabel("0", false), legendLabel("0", true),
        legendLabel("7", true)]));""")
    assert out == ["absent", "absent (not checked)", "absent, DNA absent"]
