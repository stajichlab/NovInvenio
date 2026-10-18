"""make_index_report.py: species tree section of report.html."""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "bin" / "make_index_report.py"

CONFIG = (
    "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
    "IN,Sp a,s1,a.faa,a.fna,a,fam\n"
    "IN,Sp b,s2,b.faa,b.fna,b,fam\n"
    "OUT,Sp o,s3,o.faa,o.fna,o,fam\n"
)
SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><g/></svg>\n'


def _run(tmp_path, svg=None, status=None):
    cfg = tmp_path / "config.csv"
    cfg.write_text(CONFIG)
    out = tmp_path / "report.html"
    cmd = [sys.executable, str(SCRIPT), "--config", str(cfg), "--output", str(out)]
    if svg is not None:
        p = tmp_path / "species_tree.svg"
        p.write_text(svg)
        cmd += ["--species_tree_svg", str(p)]
    if status is not None:
        p = tmp_path / "species_tree.status.txt"
        p.write_text(status)
        cmd += ["--species_tree_status", str(p)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out.read_text()


def _ok(method="outgroup", n_out=1, outside=0, markers=10):
    return (f"ok\nmarkers: {markers}\n"
            f"rooting: method={method} n_out={n_out} n_out_outside={outside}\n")


def test_ok_status_embeds_svg_and_method_line(tmp_path):
    html = _run(tmp_path, SVG, _ok())
    assert 'id="species-tree"' in html
    assert "<svg" in html
    assert "FastTree, 10 markers (top treeness/RCV), fungi_odb12" in html


def test_marker_count_comes_from_status(tmp_path):
    html = _run(tmp_path, SVG, _ok(markers=7))
    assert "FastTree, 7 markers" in html
    assert "FastTree, 10 markers" not in html


def test_failed_status_shows_reason_and_no_svg(tmp_path):
    html = _run(tmp_path, SVG, "failed: no shared markers\n")
    assert 'id="species-tree"' in html
    assert "tree not available: no shared markers" in html
    assert "<svg" not in html


def test_skipped_status_shows_reason(tmp_path):
    html = _run(tmp_path, "", "skipped: fewer than 4 species\n")
    assert "tree not available: fewer than 4 species" in html


def test_no_arguments_omits_section(tmp_path):
    html = _run(tmp_path)
    assert "species-tree" not in html


def test_only_status_given_still_renders_section(tmp_path):
    html = _run(tmp_path, None, "skipped: disabled\n")
    assert 'id="species-tree"' in html
    assert "tree not available: disabled" in html


def test_ok_status_with_empty_svg_degrades(tmp_path):
    html = _run(tmp_path, "", _ok())
    assert "tree not available" in html
    assert "<svg" not in html


def test_ok_status_with_missing_svg_degrades(tmp_path):
    html = _run(tmp_path, None, _ok())
    assert "tree not available" in html


def test_status_text_is_escaped(tmp_path):
    html = _run(tmp_path, None, "failed: <script>alert(1)</script>\n")
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html


def test_rooting_outgroup(tmp_path):
    html = _run(tmp_path, SVG, _ok("outgroup", 1, 0))
    assert "rooted on OUT group (monophyletic)" in html


def test_rooting_best_split(tmp_path):
    html = _run(tmp_path, SVG, _ok("outgroup_best_split", 5, 2))
    assert "rooted on best OUT/non-OUT split: 2 of 5 OUT tips outside the OUT side" in html


def test_rooting_best_split_zero_outside(tmp_path):
    html = _run(tmp_path, SVG, _ok("outgroup_best_split", 5, 0))
    assert "rooted on best OUT/non-OUT split: 0 of 5 OUT tips outside the OUT side" in html


def test_rooting_midpoint(tmp_path):
    html = _run(tmp_path, SVG, _ok("midpoint", 0, 0))
    assert "midpoint rooted (no OUT tips)" in html


def test_stub_status_for_incomplete_job_renders_reason(tmp_path):
    html = _run(tmp_path, "\n", "failed: species tree job did not complete\n")
    assert "tree not available: species tree job did not complete" in html
    assert "<svg" not in html
