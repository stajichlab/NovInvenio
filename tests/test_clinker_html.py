"""lib/clinker_html.py: strip sequences from clinker's embedded data."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from clinker_html import slim_clinker_html  # noqa: E402


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
    assert slim.startswith("<html><head></head><body><div id=\"plot\"></div><script>const data=")
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
