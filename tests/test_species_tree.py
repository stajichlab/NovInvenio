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
        "OUT,y,1,p,d,o1,T\nOUT,y,2,p,d,o2,T\n"
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


def test_near_ingroup_not_treated_as_out():
    g = _groups(["o1", "o2"], ["a", "b"])
    g["n1"] = "NEAR_INGROUP"
    g["d1"] = "DISCOVERY_OUT"
    nwk = "((o1:1,o2:1):2,((a:1,b:1):1,n1:1):1,d1:1);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert info["n_out"] == 2
    assert info["method"] == "outgroup"
    assert info["n_out_outside"] == 0
    assert {"o1", "o2"} in _side_tips(rooted)


def test_cli_fails_when_no_tip_label_matches_config(tmp_path):
    cfg = tmp_path / "config.csv"
    cfg.write_text(
        "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
        "IN,x,1,p,d,a,T\nIN,x,2,p,d,b,T\nOUT,y,1,o1,d,o1,T\n"
    )
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((o_1:1,a_:1):2,b_:1,c_:1);\n")
    r = subprocess.run([sys.executable, str(SCRIPT), "--newick", str(nwk),
                        "--config", str(cfg), "--out_prefix", str(tmp_path / "x")],
                       capture_output=True, text=True)
    assert r.returncode == 3
    assert "tree tip labels do not match config Short names" in r.stderr


def test_impure_split_reports_non_out_inside():
    # x sits between o1 and o2; the best OUT split puts x on the OUT side
    g = _groups(["o1", "o2"], ["x", "a", "b", "c"])
    nwk = "((o1:1,x:1):1,o2:1,(a:1,b:1,c:1):9);"
    rooted, info = st.root_on_outgroup(nwk, g)
    assert info["method"] == "outgroup_best_split"
    assert info["n_out_outside"] == 0
    assert info["n_non_out_inside"] == 1


def test_pure_split_has_zero_non_out_inside():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    _, info = st.root_on_outgroup("((o1:1,o2:1):2,(a:1,b:1):1,(c:1,d:1):1);", g)
    assert info["n_non_out_inside"] == 0


# ---- tip display labels ---------------------------------------------------

def _width(svg):
    return int(re.search(r'<svg[^>]* width="(\d+)"', svg).group(1))


def test_labels_mapping_shows_display_text_not_id():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    svg = st.render_svg(MONO, g, labels={"a": "Aspergillus fumigatus", "o1": "Out one"})
    assert ">Aspergillus fumigatus</text>" in svg
    assert ">Out one</text>" in svg
    assert ">a</text>" not in svg and ">o1</text>" not in svg


def test_labels_missing_tip_falls_back_to_id():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    svg = st.render_svg(MONO, g, labels={"a": "Name A"})
    assert ">b</text>" in svg and ">o2</text>" in svg
    assert svg.count('class="tip"') == 6


def test_labels_are_xml_escaped():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    svg = st.render_svg(MONO, g, labels={"a": "A & B <x> \"q\" / z"})
    assert "A &amp; B &lt;x&gt;" in svg
    assert "A & B" not in svg
    st_xml = __import__("xml.dom.minidom", fromlist=["x"]).parseString(svg)
    assert st_xml.documentElement.tagName == "svg"


def test_width_uses_displayed_text():
    g = _groups(["o1", "o2"], ["a", "b", "c", "d"])
    plain = st.render_svg(MONO, g)
    long_ = st.render_svg(MONO, g, labels={"a": "Some very long species name strain ABC-123456"})
    assert _width(long_) > _width(plain)
    # a short display text must not leave the SVG wider than the ID-sized one
    short = st.render_svg(MONO, g, labels={t: "x" for t in g})
    assert _width(short) <= _width(plain)


def _cli_cfg(tmp_path, rows):
    cfg = tmp_path / "config.csv"
    cfg.write_text("GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n" + "".join(
        f"{grp},{sp},{stn},p,d,{sh},T\n" for grp, sp, stn, sh in rows))
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((o1:1,o2:1):2,(a:1,b:1):1,c:1);\n")
    return cfg, nwk


def _run_cli(tmp_path, cfg, nwk, prefix, *extra):
    subprocess.run([sys.executable, str(SCRIPT), "--newick", str(nwk), "--config", str(cfg),
                    "--out_prefix", str(tmp_path / prefix), *extra], check=True)
    return (tmp_path / f"{prefix}.svg").read_text()


_ROWS = [("IN", "Sp x", "S1", "a"), ("IN", "Sp x", "", "b"), ("IN", "Sp z", "S3", "c"),
         ("OUT", "", "S4", "o1"), ("OUT", "Sp y", "S5", "o2")]


def test_cli_tip_labels_species_disambiguates(tmp_path):
    cfg, nwk = _cli_cfg(tmp_path, _ROWS)
    svg = _run_cli(tmp_path, cfg, nwk, "sp", "--tip_labels", "species")
    assert ">Sp x S1</text>" in svg       # duplicate species, strain present
    assert ">Sp x (b)</text>" in svg      # duplicate species, strain empty
    assert ">Sp z</text>" in svg          # unique species
    assert ">o1</text>" in svg            # empty Species -> Short
    assert ">Sp y</text>" in svg
    assert ">a</text>" not in svg


def test_cli_default_labels_are_short_ids(tmp_path):
    cfg, nwk = _cli_cfg(tmp_path, _ROWS)
    svg = _run_cli(tmp_path, cfg, nwk, "dflt")
    for t in ("a", "b", "c", "o1", "o2"):
        assert f">{t}</text>" in svg
    assert "Sp x" not in svg
    explicit = _run_cli(tmp_path, cfg, nwk, "short", "--tip_labels", "short")
    assert explicit == svg


def test_cli_rooting_identical_between_modes(tmp_path):
    cfg, nwk = _cli_cfg(tmp_path, _ROWS)
    _run_cli(tmp_path, cfg, nwk, "m_short", "--tip_labels", "short")
    _run_cli(tmp_path, cfg, nwk, "m_sp", "--tip_labels", "species")
    for ext in ("rooted.nwk", "rooting.json"):
        assert (tmp_path / f"m_short.{ext}").read_text() == (tmp_path / f"m_sp.{ext}").read_text()


def test_cli_rejects_unknown_tip_labels_value(tmp_path):
    cfg, nwk = _cli_cfg(tmp_path, _ROWS)
    r = subprocess.run([sys.executable, str(SCRIPT), "--newick", str(nwk), "--config", str(cfg),
                        "--tip_labels", "bogus"], capture_output=True, text=True)
    assert r.returncode != 0
