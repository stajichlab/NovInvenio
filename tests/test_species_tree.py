import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "bin" / "species_tree.py"


def _load():
    spec = importlib.util.spec_from_file_location("species_tree", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


st = _load()


def _groups(out, ins):
    g = {t: "OUT" for t in out}
    g.update({t: "IN" for t in ins})
    return g


def _side_tips(newick):
    """Tip sets under the two children of the root of a rooted Newick."""
    root = st.parse_newick(newick)
    assert len(root.children) == 2
    return [set(st.tip_labels(c)) for c in root.children]


MONO = "((o1:1,o2:1):2,((a:1,b:1):1,(c:1,d:1):1):1);"


def test_monophyletic_out_roots_on_out_clade():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    # unrooted presentation: root is a trifurcation
    nwk = "((o1:1,o2:1):2,(a:1,b:1):1,(c:1,d:1):1);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert info["method"] == "outgroup"
    assert info["n_out"] == 2
    assert info["n_out_outside"] == 0
    assert {"o1", "o2"} in _side_tips(rooted)


def test_non_monophyletic_out_reports_outside_count():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    nwk = "((o1:1,a:1):1,(b:1,c:1):1,(o2:1,d:1):1);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert info["method"] == "outgroup_best_split"
    assert info["n_out"] == 2
    assert info["n_out_outside"] >= 1
    assert set(st.tip_labels(st.parse_newick(rooted))) == {"o1", "o2", "a", "b", "c", "d"}


def test_no_out_tips_midpoint():
    g = _groups([], ["a", "b", "c", "d"])
    nwk = "((a:1,b:1):1,c:1,d:6);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert info["method"] == "midpoint"
    assert info["n_out"] == 0
    # longest path is d..a/b = 8; midpoint 4 from d lies on the d edge
    root = st.parse_newick(rooted)
    assert {"d"} in [set(st.tip_labels(c)) for c in root.children]


def test_rooting_preserves_tip_distances():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    nwk = "((o1:1,o2:1):2,(a:1,b:1):1,(c:1,d:3):1);"
    rooted, _ = st.root_on_outgroup(nwk, g)
    assert st.tip_distance(st.parse_newick(nwk), "o1", "d") == \
        st.tip_distance(st.parse_newick(rooted), "o1", "d")


def test_labels_with_special_characters_roundtrip():
    labels = ["A.fum_1", "B|2", "C 3", "D:4"]
    g = _groups(["D:4"], ["A.fum_1", "B|2", "C 3"])
    nwk = "('A.fum_1':1,'B|2':1,('C 3':1,'D:4':1):1);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert set(st.tip_labels(st.parse_newick(rooted))) == set(labels)
    svg = st.render_svg(rooted, g)
    for lab in labels:
        assert f">{lab}</text>" in svg


def test_render_svg_has_one_label_per_tip():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    svg = st.render_svg(MONO, g)
    assert svg.lstrip().startswith("<svg")
    assert "http" not in svg.replace("http://www.w3.org/2000/svg", "")
    for t in g:
        assert len(re.findall(rf">{t}</text>", svg)) == 1
    assert svg.count('class="tip"') == 6
    assert "OUT" in svg and "IN" in svg  # legend


def test_cli_writes_outputs(tmp_path):
    cfg = tmp_path / "config.csv"
    cfg.write_text(
        "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
        "IN,x,1,p,d,a,T\nIN,x,2,p,d,b,T\nIN,x,3,p,d,c,T\n"
        "DISC_OUT,y,1,p,d,o1,T\nDISC_OUT,y,2,p,d,o2,T\n"
    )
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((o1:1,o2:1):2,(a:1,b:1):1,c:1);\n")
    prefix = tmp_path / "species_tree"
    subprocess.run([sys.executable, str(SCRIPT), "--newick", str(nwk),
                    "--config", str(cfg), "--out_prefix", str(prefix)], check=True)
    assert (tmp_path / "species_tree.rooted.nwk").read_text().strip().endswith(";")
    assert (tmp_path / "species_tree.svg").read_text().lstrip().startswith("<svg")
    info = json.loads((tmp_path / "species_tree.rooting.json").read_text())
    assert info["method"] == "outgroup" and info["n_out"] == 2
