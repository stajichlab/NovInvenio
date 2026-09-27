"""bin/pangenome_island_clinker.py with a stand-in clinker executable."""
import stat
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
from pangenome_island_clinker import main  # noqa: E402

# Writes a minimal clinker-shaped page to the -p argument; fails for L002.
FAKE = """#!/usr/bin/env python3
import json, sys
args = sys.argv[1:]
out = args[args.index("-p") + 1]
if "L002" in args[0]:
    sys.exit(3)
gene = {"uid": "g", "label": "g", "sequence": "ATG" * 30, "translation": "M" * 30}
data = {"clusters": [{"loci": [{"genes": [gene]}]}], "links": [], "groups": []}
page = ("<html><head></head><body><div id=\\"div-floater\\">"
       "<button class=\\"collapsible active\\" onclick=\\"toggleActive()\\">clinker</button>"
       "<div id=\\"div-summary\\"></div></div>"
       "<script>const data=" + json.dumps(data) + ";plot(data)</script></body></html>")
open(out, "w").write(page)
"""


def setup(tmp_path):
    fake = tmp_path / "fake_clinker"
    fake.write_text(FAKE)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    dirs = []
    for key in ("L001", "L002", "L003"):
        d = tmp_path / "gbk" / key
        d.mkdir(parents=True)
        (d / "S1.gbk").write_text("LOCUS\n")
        (d / "groups.csv").write_text("g,famA\n")
        dirs.append(str(d))
    (tmp_path / "gbk" / "L003" / "S1.gbk").unlink()
    return fake, dirs


def test_each_locus_gets_a_slim_page_and_failures_are_skipped(tmp_path, capsys):
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    assert main(["--locus_dirs", *dirs, "--out_dir", str(out), "--clinker", str(fake)]) == 0
    assert sorted(p.name for p in out.iterdir()) == ["L001.html"]
    assert "ATGATG" not in (out / "L001.html").read_text()
    err = capsys.readouterr().err
    assert "L002: clinker exited 3" in err and "L003: no .gbk files" in err


def test_keep_sequences_writes_the_page_unchanged(tmp_path):
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    main(["--locus_dirs", dirs[0], "--out_dir", str(out), "--clinker", str(fake), "--keep_sequences"])
    assert "ATGATG" in (out / "L001.html").read_text()


def test_keep_sequences_still_applies_the_ui_fixes(tmp_path):
    # C2/C3: --keep_sequences skips slim_clinker_html() (and its sequences
    # stay in), but the sidebar-collapse/scrolling/left-margin fixes must
    # still apply -- they come from lib/clinker_html.inject_ui_fixes(),
    # called directly since slimming is skipped.
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    main(["--locus_dirs", dirs[0], "--out_dir", str(out), "--clinker", str(fake), "--keep_sequences"])
    page = (out / "L001.html").read_text()
    assert 'class="collapsible active"' not in page
    assert '<div id="div-summary" style="display:none">' in page
    assert "<style>body{overflow:auto;margin-left:16px}</style>" in page


def test_missing_clinker_executable_skips_every_locus(tmp_path, capsys):
    _fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    assert main(["--locus_dirs", dirs[0], "--out_dir", str(out),
                 "--clinker", str(tmp_path / "no_such_clinker")]) == 0
    assert list(out.iterdir()) == []
    assert "cannot run" in capsys.readouterr().err
