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


def test_missing_clinker_executable_fails_loudly(tmp_path, capsys):
    # I2a (final review): a missing clinker executable is not a per-locus
    # failure (R15's skip-and-continue) -- it means every locus in this
    # ISLAND_CLINKER task would silently produce no page, which used to
    # look like "0 of N loci drawn" instead of the real cause (the
    # published container image has no clinker). Fail loudly instead.
    _fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    missing = tmp_path / "no_such_clinker"
    assert main(["--locus_dirs", dirs[0], "--out_dir", str(out), "--clinker", str(missing)]) == 2
    assert not out.exists()
    err = capsys.readouterr().err
    assert f"ERROR: clinker executable '{missing}' not found" in err
    assert "container image 0.5.0 has no clinker" in err
    assert "--pangenome_clinker false" in err


def test_clinker_on_path_by_name_is_accepted(tmp_path, monkeypatch):
    # shutil.which() resolution: a bare name found on PATH (not just an
    # existing file path) must not be treated as missing.
    fake, dirs = setup(tmp_path)
    monkeypatch.setenv("PATH", str(fake.parent) + ":" + __import__("os").environ["PATH"])
    out = tmp_path / "out"
    assert main(["--locus_dirs", dirs[0], "--out_dir", str(out), "--clinker", fake.name]) == 0
    assert sorted(p.name for p in out.iterdir()) == ["L001.html"]


def test_real_locus_failure_is_still_skipped_not_fatal(tmp_path, capsys):
    # R15 (per-locus skip) is unchanged for a clinker that IS present but
    # fails on one locus (L002) or has no .gbk files (L003): exit 0, only
    # those loci are missing from out_dir.
    fake, dirs = setup(tmp_path)
    out = tmp_path / "out"
    assert main(["--locus_dirs", *dirs, "--out_dir", str(out), "--clinker", str(fake)]) == 0
    assert sorted(p.name for p in out.iterdir()) == ["L001.html"]


# Writes a page with no "const data=" block at all (a clinker version whose
# markup changed, or a truncated/corrupted -p output).
MALFORMED = """#!/usr/bin/env python3
import sys
out = sys.argv[sys.argv.index("-p") + 1]
open(out, "w").write("<html><head></head><body>no data here</body></html>")
"""


def test_a_malformed_clinker_page_is_skipped_not_fatal(tmp_path, capsys):
    # M1 (final review): slim_clinker_html() raises ValueError on a page
    # with no "const data=" block. That must not abort the whole run (R15)
    # -- it is a per-locus failure like any other clinker exit-code
    # failure, just caught one step later.
    fake = tmp_path / "malformed_clinker"
    fake.write_text(MALFORMED)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    d = tmp_path / "gbk" / "L009"
    d.mkdir(parents=True)
    (d / "S1.gbk").write_text("LOCUS\n")
    (d / "groups.csv").write_text("g,famA\n")
    out = tmp_path / "out"
    assert main(["--locus_dirs", str(d), "--out_dir", str(out), "--clinker", str(fake)]) == 0
    assert list(out.iterdir()) == []
    assert "L009" in capsys.readouterr().err


# Records the exact argv it was called with, as one line per arg, to
# out_dir/<key>.cmdline -- lets a test inspect file order and flags without
# a mock.
RECORDING = """#!/usr/bin/env python3
import json, sys
args = sys.argv[1:]
out = args[args.index("-p") + 1]
with open(out + ".cmdline", "w") as fh:
    fh.write("\\n".join(args))
gene = {"uid": "g", "label": "g"}
data = {"clusters": [{"loci": [{"genes": [gene]}]}], "links": [], "groups": []}
open(out, "w").write("<html><head></head><body><script>const data=" + json.dumps(data) +
                     ";plot(data)</script></body></html>")
"""


def test_gbk_files_are_passed_in_pick_order_with_use_file_order(tmp_path):
    # M7 (final review): bin/pangenome_island_gbk_slice.py writes
    # clinker_order.txt in clinker pick order (exemplar first); the .gbk
    # files must be passed to clinker in that order, not alphabetical, and
    # clinker's -ufo/--use_file_order flag must be present so it actually
    # respects that order.
    fake = tmp_path / "recording_clinker"
    fake.write_text(RECORDING)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    d = tmp_path / "gbk" / "L001"
    d.mkdir(parents=True)
    (d / "S2.gbk").write_text("LOCUS\n")
    (d / "S1.gbk").write_text("LOCUS\n")
    (d / "groups.csv").write_text("g,famA\n")
    (d / "clinker_order.txt").write_text("S2.gbk\nS1.gbk\n")
    out = tmp_path / "out"
    assert main(["--locus_dirs", str(d), "--out_dir", str(out), "--clinker", str(fake)]) == 0
    argv = (out / "L001.raw.html.cmdline").read_text().split("\n")
    gbk_args = [a for a in argv if a.endswith(".gbk")]
    assert [Path(a).name for a in gbk_args] == ["S2.gbk", "S1.gbk"]
    assert "-ufo" in argv


def test_missing_order_file_falls_back_to_alphabetical(tmp_path):
    fake = tmp_path / "recording_clinker"
    fake.write_text(RECORDING)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    d = tmp_path / "gbk" / "L001"
    d.mkdir(parents=True)
    (d / "S2.gbk").write_text("LOCUS\n")
    (d / "S1.gbk").write_text("LOCUS\n")
    (d / "groups.csv").write_text("g,famA\n")
    out = tmp_path / "out"
    assert main(["--locus_dirs", str(d), "--out_dir", str(out), "--clinker", str(fake)]) == 0
    argv = (out / "L001.raw.html.cmdline").read_text().split("\n")
    gbk_args = [a for a in argv if a.endswith(".gbk")]
    assert [Path(a).name for a in gbk_args] == ["S1.gbk", "S2.gbk"]
    assert "-ufo" in argv
