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
import io
import json
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


# ---- C1: a multi-record GenBank file is one cluster of several loci -----
# split_blocks()/build_record(block_label=...) (lib/genbank_slice.py) write
# every gene-free-gap block of one strain's region as its own record in the
# same <strain>.gbk file, named "b1", "b2", ... (review fix round 1, item 1:
# never "<strain>_bN" -- that truncates identically for a long strain name).
# This checks clinker 0.0.32 actually reads that as intended -- one cluster,
# named after the FILE (so still identifiable as the strain) -- with one
# locus per block, each locus named after its record -- rather than, say,
# three unrelated clusters.
@pytest.mark.skipif(not shutil.which("clinker"), reason="needs clinker on PATH")
def test_multi_record_gbk_is_one_cluster_with_several_loci(tmp_path):
    from genbank_slice import build_record, cds_exons, rank_entries  # local import: bin/lib on sys.path

    gff = ("##gff-version 3\n"
          "c1\tsrc\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=p1\n"
          "c1\tsrc\tCDS\t100000\t100300\t.\t+\t0\tID=cds2;Parent=p2\n")
    gene_rows = [("S1", "p1", "c1", 1, 300), ("S1", "p2", "c1", 100000, 100300)]
    entries = rank_entries(iter(gene_rows), iter(()), {("S1", "c1")})[("S1", "c1")]
    exons = cds_exons(io.StringIO(gff), {"c1"}, {"p1", "p2"})
    recs = []
    labels = {}
    for i, entry in enumerate(entries):
        rec, summ = build_record("S1", "c1", "A" * 100400, [entry], exons, {}, {"p1": "famA",
                                 "p2": "famA"}, lambda pid: pid, block_label=f"b{i + 1}")
        recs.append(rec)
        labels.update(summ["labels"])
    SeqIO.write(recs, str(tmp_path / "S1.gbk"), "genbank")
    with open(tmp_path / "groups.csv", "w") as fh:
        for label, fam in sorted(labels.items()):
            fh.write(f"{label},{fam}\n")
    out = tmp_path / "out.html"
    subprocess.run(["clinker", str(tmp_path / "S1.gbk"), "-gf", str(tmp_path / "groups.csv"),
                    "-j", "1", "-p", str(out)], check=True, capture_output=True)
    html = out.read_text()
    start = html.index("const data=") + len("const data=")
    data = json.JSONDecoder().raw_decode(html, start)[0]
    assert len(data["clusters"]) == 1
    cluster = data["clusters"][0]
    assert cluster["name"] == "S1"
    assert len(cluster["loci"]) == 2
    assert sorted(locus["name"] for locus in cluster["loci"]) == ["b1", "b2"]
    assert [len(locus["genes"]) for locus in
           sorted(cluster["loci"], key=lambda x: x["name"])] == [1, 1]
