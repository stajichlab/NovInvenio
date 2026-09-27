"""lib/clinker_html.py: strip sequences from clinker's embedded data."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import inject_ui_fixes, slim_clinker_html  # noqa: E402


def gene(uid):
    return {"uid": uid, "label": uid, "names": {"locus_tag": uid}, "start": 0, "end": 9,
            "strand": 1, "sequence": "ATG" * 50, "translation": "M" * 50}


DATA = {"clusters": [{"uid": "c1", "name": "S1", "loci": [{"uid": "l1", "genes": [gene("g1")]}]}],
        "links": [{"uid": "k1", "query": gene("g1"), "target": gene("g2"), "identity": 0.9}],
        "groups": [{"uid": "gr1", "label": "famA", "genes": ["g1", "g2"]}]}
PAGE = ("<html><head></head><body><div id=\"plot\"></div><script>const data=" +
        json.dumps(DATA) + ";function serialise(svg) { return 1; }\nplot(data)</script></body></html>")


def embedded(html):
    start = html.index("const data=") + len("const data=")
    return json.JSONDecoder().raw_decode(html, start)[0]


def test_sequences_are_removed_from_clusters_and_links():
    data = embedded(slim_clinker_html(PAGE))
    genes = [data["clusters"][0]["loci"][0]["genes"][0], data["links"][0]["query"],
             data["links"][0]["target"]]
    for g in genes:
        assert "sequence" not in g and "translation" not in g
        assert g["label"] and g["end"] == 9


def test_everything_else_is_kept():
    slim = slim_clinker_html(PAGE)
    # C2/C3: slim_clinker_html also applies inject_ui_fixes()'s <style> (see
    # its own tests below); PAGE's <head></head> has no sidebar markup to
    # collapse.
    assert slim.startswith(
        "<html><head><style>body{overflow:auto;margin-left:16px}</style></head>"
        "<body><div id=\"plot\"></div><script>const data=")
    assert slim.endswith(";function serialise(svg) { return 1; }\nplot(data)</script></body></html>")
    data = embedded(slim)
    assert data["groups"] == DATA["groups"]
    assert data["links"][0]["identity"] == 0.9


def test_a_closing_script_tag_in_a_label_is_escaped():
    d = json.loads(json.dumps(DATA))
    d["groups"][0]["label"] = "</script><b>x"
    page = PAGE.replace(json.dumps(DATA), json.dumps(d))
    blob = slim_clinker_html(page).split("const data=", 1)[1].split(";function", 1)[0]
    assert "</script>" not in blob


def test_non_clinker_page_raises():
    with pytest.raises(ValueError):
        slim_clinker_html("<html></html>")


# ---- C2/C3: sidebar starts collapsed; scrolling + left margin (Part B) ----
# Markup verified 2026-09-27 against a real clinker 0.0.32 page (spike run
# under $SCRATCH/clinker-venv): the sidebar is #div-floater, its toggle is
# <button class="collapsible active" onclick="toggleActive()">clinker</button>
# followed by <div id="div-summary"> (clinker's own unmodified toggleActive()
# reads menu.style.display, so setting it inline here is what makes the
# very next click reopen the panel correctly). The body sets
# `overflow: hidden`, which clips a figure taller than the iframe.
REAL_CLINKER_PAGE = (
    "<html><head><style>body {\n  margin: 0;\n  padding: 0;\n  overflow: hidden;\n"
    "  font-family: sans-serif;\n}\n</style></head><body><main><div id='plot'></div>"
    "<div id=\"div-floater\">\n"
    "      <button class=\"collapsible active\" type=\"button\" onclick=\"toggleActive()\">"
    "clinker</button>\n      <div id=\"div-summary\">\n"
    "        <p id=\"p-result-summary\"></p>\n      </div>\n    </div>\n"
    "    <script>\n      const coll = document.querySelector(\".collapsible\");\n"
    "      const menu = document.querySelector(\"#div-summary\");\n"
    "      function toggleActive(event) {\n        coll.classList.toggle(\"active\");\n"
    "        menu.style.display = menu.style.display === \"none\" ? \"block\" : \"none\";\n"
    "      }\n    </script>\n    <script>const data=" + json.dumps(DATA) +
    ";</script></main></body></html>"
)


def test_sidebar_starts_collapsed_but_its_button_still_works():
    out = inject_ui_fixes(REAL_CLINKER_PAGE)
    assert 'class="collapsible active"' not in out
    assert out.count('class="collapsible"') == 1
    assert out.count('<div id="div-summary" style="display:none">') == 1
    # clinker's own toggleActive() is untouched: menu.style.display (the
    # inline style we just set) is what its first click reads.
    assert 'menu.style.display === "none" ? "block" : "none"' in out


def test_body_overflow_is_no_longer_hidden_and_page_still_has_a_left_margin():
    out = inject_ui_fixes(REAL_CLINKER_PAGE)
    assert out.count("<style>body{overflow:auto;margin-left:16px}</style>") == 1
    assert out.index("<style>body{overflow:auto;margin-left:16px}</style>") < out.index("</head>")


def test_inject_ui_fixes_is_idempotent_on_a_page_missing_the_markers():
    # A no-op (not a crash) if a future clinker version's markup differs.
    assert inject_ui_fixes("<html><head></head><body>x</body></html>") == (
        "<html><head><style>body{overflow:auto;margin-left:16px}</style></head><body>x</body></html>")


def test_slim_clinker_html_also_applies_the_ui_fixes():
    out = slim_clinker_html(REAL_CLINKER_PAGE)
    assert 'class="collapsible active"' not in out
    assert '<div id="div-summary" style="display:none">' in out
    assert "<style>body{overflow:auto;margin-left:16px}</style>" in out
    # and the data blob was still slimmed.
    assert "sequence" not in embedded(out)["clusters"][0]["loci"][0]["genes"][0]
