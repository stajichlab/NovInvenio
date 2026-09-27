"""Slimmed clinker pages render the same figure as the full page.

Spec section 8: the 0.81 MB slimmed page was not yet checked to render and
behave the same. This runs clinker 0.0.32 on three small GenBank files,
slims the page with lib/clinker_html.py, loads both pages in a
headless Chromium (`--dump-dom`, which returns the DOM after scripts ran)
and compares the drawn gene, cluster and link-path counts.

Skips unless `clinker` is on PATH and a headless Chromium is found:
NOVINVENIO_HEADLESS_CHROME, else the newest
~/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell.
"""
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import FeatureLocation, SeqFeature  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402
from clinker_html import slim_clinker_html  # noqa: E402


def find_chrome() -> str | None:
    env = os.environ.get("NOVINVENIO_HEADLESS_CHROME")
    if env and Path(env).is_file():
        return env
    hits = sorted(Path.home().glob(".cache/ms-playwright/chromium_headless_shell-*/"
                                   "chrome-headless-shell-linux64/chrome-headless-shell"))
    return str(hits[-1]) if hits else None


CHROME = find_chrome()
pytestmark = pytest.mark.skipif(not (shutil.which("clinker") and CHROME),
                                reason="needs clinker on PATH and a headless Chromium")


def write_gbks(d: Path) -> None:
    rng = random.Random(7)
    prots = ["M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(120)) for _ in range(4)]
    for strain, keep in (("S1", [0, 1, 2, 3]), ("S2", [0, 1, 3]), ("S3", [0, 3])):
        feats, pos = [], 0
        for i in keep:
            loc = FeatureLocation(pos, pos + 3 * len(prots[i]) + 3, strand=1)
            tag = f"{strain}_g{i}"
            feats.append(SeqFeature(loc, type="gene", qualifiers={"locus_tag": [tag]}))
            feats.append(SeqFeature(loc, type="CDS", qualifiers={
                "locus_tag": [tag], "translation": [prots[i]], "note": [f"family=f{i}"]}))
            pos += 3 * len(prots[i]) + 200
        rec = SeqRecord(Seq("A" * pos), id=strain, name=strain, features=feats)
        rec.annotations["molecule_type"] = "DNA"
        SeqIO.write(rec, str(d / f"{strain}.gbk"), "genbank")
    with open(d / "groups.csv", "w") as fh:
        for strain, keep in (("S1", [0, 1, 2, 3]), ("S2", [0, 1, 3]), ("S3", [0, 3])):
            for i in keep:
                fh.write(f"{strain}_g{i},f{i}\n")


def dom_counts(page: Path) -> dict:
    dom = subprocess.run([CHROME, "--headless", "--no-sandbox", "--disable-gpu",
                          "--virtual-time-budget=10000", "--window-size=1600,1000",
                          "--dump-dom", str(page)],
                         capture_output=True, text=True, timeout=120).stdout
    return {k: dom.count(k) for k in ('class="gene"', 'class="cluster"', "<path")}


def test_slimmed_page_draws_the_same_figure(tmp_path):
    write_gbks(tmp_path)
    full = tmp_path / "full.html"
    subprocess.run(["clinker", *sorted(str(p) for p in tmp_path.glob("*.gbk")),
                    "-gf", str(tmp_path / "groups.csv"), "-j", "1", "-p", str(full)],
                   check=True, capture_output=True)
    slim = tmp_path / "slim.html"
    slim.write_text(slim_clinker_html(full.read_text()))
    assert slim.stat().st_size < full.stat().st_size
    before, after = dom_counts(full), dom_counts(slim)
    assert before['class="gene"'] == 9 and before['class="cluster"'] == 3
    assert after == before
