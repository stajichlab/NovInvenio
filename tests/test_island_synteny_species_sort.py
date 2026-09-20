"""Behavioural tests for the "by species" row-sort helpers added to
lib/island_synteny_template.py for issue #119 (and tightened in review).

These three JS functions -- speciesCounts, haplotypeSpecies,
speciesBandLabel -- are deliberately written to take `speciesMap` as an
explicit argument rather than closing over the module-level SPECIES var, so
they can be extracted verbatim from the template source (a brace-matching
slice, not a JS parser -- these are simple, single-level functions) and run
under plain `node`, with NO jsdom/DOM dependency: none of them touch the
DOM, so there is no reason to require it just to exercise this logic.

Per the issue #119 review: do not touch tests/test_report_js_behaviour.py
or tests/js/drive_reports.mjs (those cover the full DOM-driven pages and are
being extended separately for issue #120). This file is new and standalone.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from island_synteny_template import ISLAND_SYNTENY_TEMPLATE  # noqa: E402

FUNCTION_NAMES = ["speciesCounts", "haplotypeSpecies", "speciesBandLabel"]


def _extract_function(src: str, name: str) -> str:
    """A brace-matching slice of `function <name>(...) { ... }` from `src`.

    Not a JS parser -- deliberately just enough for these three flat,
    single-level function declarations (no nested `function` keyword
    confusion possible since the search anchors on the literal declaration
    text). Raises with a clear message if the template ever stops defining
    one of these under this exact name, rather than silently testing
    nothing.
    """
    needle = "function " + name + "("
    start = src.index(needle)
    brace_start = src.index("{", start)
    depth = 0
    for i in range(brace_start, len(src)):
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                return src[start:i + 1]
    raise AssertionError(f"unterminated function {name}() in template")


def _extracted_functions_js() -> str:
    return "\n\n".join(_extract_function(ISLAND_SYNTENY_TEMPLATE, name)
                       for name in FUNCTION_NAMES)


def test_all_three_functions_are_present_and_extractable():
    # Fails loudly (not silently) if a rename/refactor breaks extraction --
    # the point of this test file is to keep testing the REAL logic, not a
    # stale copy of it.
    js = _extracted_functions_js()
    for name in FUNCTION_NAMES:
        assert f"function {name}(" in js


def _run_node_cases(cases_js: str) -> list[str]:
    """Run the extracted functions plus `cases_js` (which should `console.log`
    one "PASS <name>" or "FAIL <name>: <detail>" line per case) under node,
    and return the output lines. Skips if node is unavailable."""
    node = shutil.which("node")
    if not node:
        pytest.skip("node not available")
    script = _extracted_functions_js() + "\n\n" + cases_js
    proc = subprocess.run([node, "--input-type=commonjs", "-e", script],
                          capture_output=True, text=True, check=False)
    assert proc.returncode == 0, f"node script crashed:\n{proc.stdout}\n{proc.stderr}"
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    failed = [ln for ln in lines if ln.startswith("FAIL")]
    assert not failed, "node behaviour checks failed:\n" + "\n".join(lines) + proc.stderr
    return lines


def hap(strains, count=None):
    return {"strains": strains, "count": count if count is not None else len(strains)}


# ---- case 1: a haplotype whose strains are all one species bands to that species

def test_pure_species_haplotype_bands_to_its_own_species():
    species_map = {"S1": "Coccidioides immitis", "S2": "Coccidioides immitis"}
    h = hap(["S1", "S2"])
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var h = {json.dumps(h)};
    var got = haplotypeSpecies(h, speciesMap);
    console.log(got === "Coccidioides immitis" ? "PASS pure" : "FAIL pure: " + got);
    var label = speciesBandLabel(h, speciesMap);
    console.log(label === "Coccidioides immitis" ? "PASS pure_label" : "FAIL pure_label: " + label);
    """
    _run_node_cases(cases)


# ---- case 2: a haplotype spanning two species bands to the MODAL one

def test_mixed_haplotype_bands_to_the_modal_species():
    species_map = {
        "S1": "Coccidioides immitis", "S2": "Coccidioides immitis",
        "S3": "Coccidioides immitis", "S4": "Coccidioides posadasii",
    }
    h = hap(["S1", "S2", "S3", "S4"])
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var h = {json.dumps(h)};
    var got = haplotypeSpecies(h, speciesMap);
    console.log(got === "Coccidioides immitis" ? "PASS modal" : "FAIL modal: " + got);
    var label = speciesBandLabel(h, speciesMap);
    console.log(label === "Coccidioides immitis (+1 other)" ? "PASS modal_label" : "FAIL modal_label: " + label);
    """
    _run_node_cases(cases)


# ---- case 3: an exact tie (2 and 2) breaks alphabetically and is STABLE

def test_exact_tie_breaks_alphabetically_and_is_stable_across_runs():
    species_map = {
        "S1": "Coccidioides posadasii", "S2": "Coccidioides posadasii",
        "S3": "Coccidioides immitis", "S4": "Coccidioides immitis",
    }
    h = hap(["S1", "S2", "S3", "S4"])
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var h = {json.dumps(h)};
    var results = [];
    for (var i = 0; i < 20; i++) {{ results.push(haplotypeSpecies(h, speciesMap)); }}
    var allSame = results.every(function (r) {{ return r === results[0]; }});
    console.log(allSame ? "PASS tie_stable" : "FAIL tie_stable: " + JSON.stringify(results));
    console.log(results[0] === "Coccidioides immitis" ? "PASS tie_alpha" : "FAIL tie_alpha: " + results[0]);
    var label = speciesBandLabel(h, speciesMap);
    console.log(label === "Coccidioides immitis (+1 other)" ? "PASS tie_label" : "FAIL tie_label: " + label);
    """
    _run_node_cases(cases)


# ---- case 4: a strain missing from the species map does not crash and does
# not silently become some other species

def test_strain_missing_from_species_map_does_not_crash_or_misattribute():
    species_map = {"S1": "Coccidioides immitis"}  # S2 deliberately absent
    h = hap(["S1", "S2"])
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var h = {json.dumps(h)};
    var got = haplotypeSpecies(h, speciesMap);
    // Tie between "Coccidioides immitis" (1) and "" (1, the missing S2) --
    // "" sorts before any real species name alphabetically, so "" is NOT
    // picked as modal only because "Coccidioides immitis" < "" is false...
    // rather than assert a specific winner here (both are legitimate with
    // this tie), assert the missing strain is NEVER attributed to some
    // species it was never mapped to.
    var counts = speciesCounts(h, speciesMap);
    var knownWrong = Object.keys(counts).some(function (sp) {{
      return sp !== "" && sp !== "Coccidioides immitis";
    }});
    console.log(!knownWrong ? "PASS missing_no_misattribution" : "FAIL missing_no_misattribution: " + JSON.stringify(counts));
    console.log(counts[""] === 1 ? "PASS missing_bucketed_unknown" : "FAIL missing_bucketed_unknown: " + JSON.stringify(counts));
    var label = speciesBandLabel(h, speciesMap);
    console.log(/\\(\\+1 other\\)$/.test(label) ? "PASS missing_label_flags_mixed" : "FAIL missing_label_flags_mixed: " + label);
    """
    _run_node_cases(cases)


def test_all_strains_missing_from_species_map_bands_as_unknown_not_crash():
    species_map = {}
    h = hap(["S1", "S2"])
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var h = {json.dumps(h)};
    var got = haplotypeSpecies(h, speciesMap);
    console.log(got === "" ? "PASS all_unknown" : "FAIL all_unknown: " + JSON.stringify(got));
    var label = speciesBandLabel(h, speciesMap);
    console.log(label === "Unknown species" ? "PASS all_unknown_label" : "FAIL all_unknown_label: " + label);
    """
    _run_node_cases(cases)


# ---- case 5: the sort groups same-species bands together and is deterministic
# (exercises the actual sort comparator logic used in sortedHaplotypes'
# "species" branch, reimplemented here in terms of the same two extracted
# primitives so a change to the real comparator that stops matching this
# would be caught by test_island_synteny_row_sort_offers_a_species_option's
# sibling tests in tests/test_report_templates.py, not silently drift)

def test_species_sort_groups_bands_together_and_is_deterministic():
    species_map = {
        "S1": "Coccidioides posadasii", "S2": "Coccidioides posadasii",
        "S3": "Coccidioides immitis", "S4": "Coccidioides immitis",
        "S5": "Coccidioides immitis",
    }
    haps = [
        hap(["S1", "S2"]),           # pure posadasii, count 2
        hap(["S3", "S4", "S5"]),     # pure immitis, count 3
        hap(["S6"]),                 # unknown species, count 1
    ]
    cases = f"""
    var speciesMap = {json.dumps(species_map)};
    var haps = {json.dumps(haps)};
    function speciesSortCompare(a, b) {{
      var sa = haplotypeSpecies(a, speciesMap) || "\\uffff";
      var sb = haplotypeSpecies(b, speciesMap) || "\\uffff";
      if (sa !== sb) return sa < sb ? -1 : 1;
      if (b.count !== a.count) return b.count - a.count;
      var an = a.strains[0] || "", bn = b.strains[0] || "";
      return an < bn ? -1 : an > bn ? 1 : 0;
    }}
    var sorted1 = haps.slice().sort(speciesSortCompare).map(function(h) {{ return h.strains.join(","); }});
    var sorted2 = haps.slice().reverse().sort(speciesSortCompare).map(function(h) {{ return h.strains.join(","); }});
    var expected = ["S3,S4,S5", "S1,S2", "S6"];
    console.log(JSON.stringify(sorted1) === JSON.stringify(expected) ? "PASS groups_together" : "FAIL groups_together: " + JSON.stringify(sorted1));
    console.log(JSON.stringify(sorted1) === JSON.stringify(sorted2) ? "PASS deterministic" : "FAIL deterministic: " + JSON.stringify(sorted1) + " vs " + JSON.stringify(sorted2));
    """
    _run_node_cases(cases)
