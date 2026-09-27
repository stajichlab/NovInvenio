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
    out = run_node(["locusIsRankKey", "locusSidebarOrder"], f"""
      var L = {json.dumps(loci)};
      console.log(JSON.stringify(["informative","strains","size","name"].map(function (k) {{
        return locusSidebarOrder(L, [0, 1, 2], k); }})));""")
    assert out == [[1, 0, 2], [2, 0, 1], [1, 0, 2], [1, 0, 2]]


# ---- Multi-rank sort (ranks-brief.md, changed 2026-09-27) --------------------

def test_sort_select_has_the_five_ranks_in_order_plus_strains_size_name():
    m = re.findall(r'<option value="([\w-]+)"', LOCUS_VIEW_HTML)
    assert m == ["presence", "within", "species", "whole_dna", "whole_annot",
                "strains", "size", "name"]
    assert 'id="lv-sort-species"' in LOCUS_VIEW_HTML


RANK_FUNCS = ["locusIsRankKey", "locusRankLabel", "locusRankValue"]


def test_rank_sort_orders_by_locus_ranks_key_with_minus_one_last():
    loci = [{"ranks": {"whole_annot": 5}}, {"ranks": {"whole_annot": 40}},
            {"ranks": {"whole_annot": -1}}]
    out = run_node(RANK_FUNCS + ["locusSidebarOrder"], f"""
      var L = {json.dumps(loci)};
      console.log(JSON.stringify(locusSidebarOrder(L, [0, 1, 2], "whole_annot")));""")
    assert out == [1, 0, 2]


def test_rank_sort_ties_keep_the_payload_order():
    loci = [{"ranks": {"presence": 5}}, {"ranks": {"presence": 5}}, {"ranks": {"presence": -1}}]
    out = run_node(RANK_FUNCS + ["locusSidebarOrder"], f"""
      var L = {json.dumps(loci)};
      console.log(JSON.stringify(locusSidebarOrder(L, [0, 1, 2], "presence")));""")
    assert out == [0, 1, 2]


def test_rank_note_says_not_informative_for_minus_one():
    out = run_node(RANK_FUNCS + ["locusRankNote"], """
      console.log(JSON.stringify(locusRankNote({ranks: {whole_dna: -1}}, "whole_dna")));""")
    assert out == "not informative for this sort"


def test_rank_note_shows_species_difference_and_within_species_name():
    out = run_node(RANK_FUNCS + ["locusRankNote"], """
      var species = locusRankNote({ranks: {species: 0.973}}, "species");
      var within = locusRankNote({ranks: {within: 12}, within_species: "Sp one"}, "within");
      console.log(JSON.stringify([species, within]));""")
    assert out[0] == "species-specific: difference 0.973"
    assert out[1] == "within-species: 12 (Sp one)"


def test_rank_note_is_blank_for_the_legacy_sort_keys():
    out = run_node(RANK_FUNCS + ["locusRankNote"], """
      console.log(JSON.stringify([locusRankNote({}, "strains"), locusRankNote({}, "size"),
        locusRankNote({}, "name")]));""")
    assert out == ["", "", ""]


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
    # V3: code 7 (absent, DNA absent -- a confirmed deletion) must not share
    # code 0's swatch (absent / not checked); it gets its own dark solid fill.
    assert out[1] == {"token": "--text-primary", "alpha": 0.7, "hatch": False}
    assert out[2].startswith("absent, DNA present") and out[3] == "absent, DNA absent"


def test_code_7_is_visually_distinct_from_code_0():
    out = run_node(["locusStateStyle"], """
      console.log(JSON.stringify([locusStateStyle("0"), locusStateStyle("7")]));""")
    assert out[0] != out[1]


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
    # V3: the note explains the hatched-grey swatch (code 6) but, before this
    # fix, said nothing about the dark solid swatch (code 7, confirmed
    # deletion) it sits right next to in the legend.
    assert "dark = DNA absent (confirmed)" in out[0]


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


# V1: hatch() drew diagonals from x-h to x+w+h with no clip, so white lines
# spilled onto neighbouring cells (a solid-blue cell shows a white diagonal
# belonging to the cell next to it).
def test_hatch_clips_to_its_own_cell():
    out = run_node(["hatch"], """
      var calls = [];
      var ctx = {
        save: function () { calls.push("save"); },
        restore: function () { calls.push("restore"); },
        beginPath: function () { calls.push("beginPath"); },
        rect: function (x, y, w, h) { calls.push("rect:" + [x, y, w, h].join(",")); },
        clip: function () { calls.push("clip"); },
        moveTo: function () { calls.push("moveTo"); },
        lineTo: function () { calls.push("lineTo"); },
        stroke: function () { calls.push("stroke"); },
      };
      hatch(ctx, 10, 20, 30, 8, "red");
      console.log(JSON.stringify(calls));""")
    assert out[:4] == ["save", "beginPath", "rect:10,20,30,8", "clip"]
    assert out.index("clip") < out.index("moveTo") < out.index("stroke")
    assert out[-1] == "restore"


# V2: rows are grouped by row_class with only a thin separator line, giving
# no visible label for what a band of rows is. locusGridSlots() inserts a
# header slot before each band's first row, carrying the class label and the
# band's total strain count (sum of row.count), so the grid-drawing and
# hover-mapping code can walk one array of fixed-height slots.
def test_grid_slots_insert_one_header_per_band_with_its_strain_count():
    rows = [
        {"row_class": "empty", "codes": "0", "count": 3, "strains": ["e1"]},
        {"row_class": "empty", "codes": "0", "count": 57, "strains": ["e2"]},
        {"row_class": "model_difference", "codes": "6", "count": 2, "strains": ["m1"]},
    ]
    out = run_node(["locusGridSlots"], f"""
      console.log(JSON.stringify(locusGridSlots({json.dumps(rows)})));""")
    assert [s["header"] for s in out] == [True, False, False, True, False]
    assert out[0]["cls"] == "empty" and out[0]["count"] == 60
    assert out[3]["cls"] == "model_difference" and out[3]["count"] == 2
    assert out[1]["row"]["strains"] == ["e1"]
    assert out[4]["row"]["strains"] == ["m1"]


def test_row_at_slot_maps_through_header_offsets_to_the_right_data_row():
    rows = [
        {"row_class": "full", "codes": "1", "count": 2, "strains": ["S1", "S2"]},
        {"row_class": "empty", "codes": "0", "count": 1, "strains": ["S3"]},
    ]
    out = run_node(["locusGridSlots", "locusRowAtSlot"], f"""
      var slots = locusGridSlots({json.dumps(rows)});
      console.log(JSON.stringify(slots.map(function (_, i) {{
        var r = locusRowAtSlot(slots, i); return r ? r.strains : null; }})));""")
    assert out == [None, ["S1", "S2"], None, ["S3"]]


# ---- clinker panel (Part B) ----

def test_clinker_panel_states():
    out = run_node(["clinkerPanelState"], """
      var on = {enabled: true, keys: ["L001"]};
      console.log(JSON.stringify([
        clinkerPanelState({key: "L001"}, on), clinkerPanelState({key: "L002"}, on),
        clinkerPanelState({key: "L001"}, {enabled: false, keys: ["L001"]}),
        clinkerPanelState({key: "../x"}, {enabled: true, keys: ["../x"]})]));""")
    assert out == [{"mode": "ok", "src": "clinker/L001.html"}, {"mode": "missing", "src": ""},
                   {"mode": "off", "src": ""}, {"mode": "missing", "src": ""}]


# N3: NII's per-run sync publishes only the top N loci' clinker pages and
# sets window.CLINKER_PUBLISHED to record which keys made the cut. Undefined
# (no `window` global at all under plain node, or the property unset in
# jsdom) must behave exactly like today (every key in `clinker.keys` gets
# its iframe); a key present in CLINKER_PUBLISHED still does; one absent
# from it switches to the new "not_published" mode instead of "ok".
def test_clinker_panel_state_respects_clinker_published():
    out = run_node(["clinkerPanelState"], """
      global.window = {};
      var on = {enabled: true, keys: ["L001", "L002"]};
      var results = [];
      results.push(clinkerPanelState({key: "L001"}, on));
      window.CLINKER_PUBLISHED = ["L001"];
      results.push(clinkerPanelState({key: "L001"}, on));
      results.push(clinkerPanelState({key: "L002"}, on));
      console.log(JSON.stringify(results));""")
    assert out == [
        {"mode": "ok", "src": "clinker/L001.html"},
        {"mode": "ok", "src": "clinker/L001.html"},
        {"mode": "not_published", "src": ""},
    ]


def test_clinker_reason_and_region_text():
    out = run_node(["locusClassLabel", "clinkerReason", "clinkerRegionText"], """
      var a = {reason: "exemplar", row_class: "full", species: "", contig: "c1", rank_lo: 3, rank_hi: 9};
      var b = {reason: "best_empty", row_class: "empty", species: "Sp two", contig: "c1",
               bp_start: 100, bp_end: 900, n_genes: 12};
      var c = {reason: "fill", row_class: "partial", species: "Sp one", contig: "c2", rank_lo: 1, rank_hi: 2};
      console.log(JSON.stringify([clinkerReason(a), clinkerReason(b), clinkerReason(c),
                                  clinkerRegionText(a), clinkerRegionText(b)]));""")
    assert out == ["locus exemplar", "best-assembled empty site strain of Sp two",
                   "more partial strains, best assembly first", "c1 gene ranks 3-9",
                   "c1:100-900, 12 genes"]


def test_clinker_panel_markup_is_below_the_grid():
    page = ISLAND_SYNTENY_TEMPLATE
    assert page.index('id="lv-grid"') < page.index('id="lv-clinker"')
    assert "Synteny (clinker)" in page
    assert "No synteny figure for this locus." in page


# ---- clinker panel gap note + strain popup (Part B, C1/C3) ----

def test_clinker_list_line_text_appends_the_gap_note_when_split():
    out = run_node(["locusClassLabel", "clinkerReason", "clinkerRegionText", "clinkerGapNote",
                    "clinkerListLineText"], """
      var whole = {strain: "S1", reason: "exemplar", row_class: "full", species: "",
                   contig: "c1", bp_start: 1, bp_end: 9000, n_genes: 5, n_blocks: 1, gap_bp: 0};
      var split = {strain: "S2", reason: "fill", row_class: "full", species: "",
                   contig: "c1", bp_start: 1, bp_end: 338667, n_genes: 5, n_blocks: 4,
                   gap_bp: 241208};
      console.log(JSON.stringify([clinkerListLineText(whole), clinkerListLineText(split)]));""")
    assert out[0] == "S1: locus exemplar; c1:1-9000, 5 genes"
    assert out[1] == ("S2: more full locus strains, best assembly first; c1:1-338667, 5 genes; "
                      "4 blocks, 241.2 kb without genes not drawn")


def test_clinker_popup_lines_use_species_over_pick_species():
    out = run_node(["locusClassLabel", "clinkerReason", "clinkerRegionText", "clinkerGapNote",
                    "clinkerPopupLines"], """
      var pick = {strain: "S1", reason: "exemplar", row_class: "full", species: "stale",
                  contig: "c1", bp_start: 100, bp_end: 900, n_genes: 12, n_blocks: 3,
                  gap_bp: 5000, drawn_bp: 700};
      console.log(JSON.stringify([clinkerPopupLines(pick, "Coccidioides immitis"),
                                  clinkerPopupLines(pick, "")]));""")
    with_species, without_species = out
    assert with_species == ["S1", "Coccidioides immitis", "locus exemplar", "c1:100-900, 12 genes",
                            "3 blocks", "3 blocks, 5.0 kb without genes not drawn",
                            "drawn: 700 bp"]
    assert without_species[:2] == ["S1", "stale"]


def test_clinker_popup_lines_unknown_species_and_unsplit_region():
    out = run_node(["locusClassLabel", "clinkerReason", "clinkerRegionText", "clinkerGapNote",
                    "clinkerPopupLines"], """
      var pick = {strain: "S1", reason: "exemplar", row_class: "full", species: "",
                  contig: "c1", rank_lo: 0, rank_hi: 9, n_blocks: 1};
      console.log(JSON.stringify(clinkerPopupLines(pick, "")));""")
    assert out == ["S1", "Unknown species", "locus exemplar", "c1 gene ranks 0-9", "1 block"]
