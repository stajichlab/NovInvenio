"""Tests for bin/pangenome_starship_overlap.py on a tiny synthetic pangenome run.

The first five tests are NII tests/test_starship_island_overlap.py (2026-10-17), run
against the ported script unchanged except for the module name and the run fixture
helper; the rest cover what the port added (single-file inputs, 0-byte stubs, the
'strain' column, the gene-level element form, --id_sep, input_checks.tsv)."""
import csv
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import pangenome_starship_overlap as sio  # noqa: E402


def _w(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def make_run(tmp: Path) -> Path:
    rd = tmp / "pangenome"
    # families F1..F7 (named by representative member)
    fams = {"F1": ("core", 10), "F2": ("shell", 6), "F3": ("shell", 5), "F4": ("cloud", 5),
            "F5": ("cloud", 2), "F6": ("core", 10), "F7": ("shell", 8), "F8": ("shell", 9)}
    rows = ["family\tfrequency\tstrain_count\tbin"] + [f"{f}\t0.5\t{n}\t{b}" for f, (b, n) in fams.items()]
    _w(rd / "frequency_table.tsv", "\n".join(rows) + "\n")
    _w(rd / "significant_islands.tsv",
       "n_strains\texample_strain\tisland_size\tmember_families\tn_supporting_pairs\tclassifications\n"
       "2\tA\t3\tF2,F3,F4\t1\tunexplained_physical\n"
       "1\tB\t2\tF5,F7\t1\tunexplained_physical\n")
    # members: A|g1..g6 ; B|x1, B|x2
    m = {"g1": "F1", "g2": "F2", "g3": "F3", "g4": "F4", "g5": "F6", "g6": "F5"}
    cl = [f"{f}\tA|{g}" for g, f in m.items()] + ["F5\tB|x1", "F7\tB|x2", "F8\tB|x3"]
    _w(rd / "cluster" / "tier1_cluster.tsv", "\n".join(cl) + "\n")
    coords = {"g1": (100, 200), "g2": (300, 400), "g3": (500, 600), "g4": (700, 800),
              "g5": (900, 1000), "g6": (1100, 1200)}
    gp = ["Short\tprotein_id\tcontig\tstart\tend"] + [f"A\t{g}\tc1\t{s}\t{e}" for g, (s, e) in coords.items()]
    gp += ["B\tx1\tc9\t100\t200", "B\tx2\tc9\t300\t400", "B\tx3\tc9\t500\t600"]
    _w(rd / "gene_positions.tsv", "\n".join(gp) + "\n")
    fp = ["Short\tfamily\tcontig\trank"] + [f"A\t{m[g]}\tc1\t{i}" for i, g in enumerate(coords)]
    fp += ["B\tF5\tc9\t0", "B\tF7\tc9\t1", "B\tF8\tc9\t2"]
    _w(rd / "family_positions.tsv", "\n".join(fp) + "\n")
    _w(rd / "rescue_positions.tsv", "Short\tfamily\tcontig\tstart\n")
    return rd


ELS = ("element_id\tlabel\tconfidence\tShort\tcontig\tstart\tend\n"
       "E1\tShip h1\thigh\tA\tc1\t480\t820\n"
       "E2\tShip h2\tunknown\tA\tc1\t1050\t1250\n")
TBL = ("# comment\nB|x1 - DUF3435 PF11917.14 1e-50 100 0 1e-50 100 0 1 1 0 0 1 1 1 1 x\n"
       "B|x3 - DUF3435 PF11917.14 1e-2 10 0 1e-2 10 0 1 1 0 0 1 1 1 1 weak\n")


def _summary(out: Path) -> dict:
    return {r["metric"]: r["value"] for r in csv.DictReader(open(out / "summary.tsv"), delimiter="\t")}


def test_overlap_end_to_end(tmp_path):
    rd = make_run(tmp_path)
    ann = tmp_path / "els.tsv"
    _w(ann, ELS)
    tbl = tmp_path / "cap.tblout"
    _w(tbl, TBL)
    out = tmp_path / "out"
    sio.main(["--run-dir", str(rd), "--starship-annotation", str(ann), "--captain-tblout", str(tbl),
              "--out-dir", str(out), "--n-perm", "20", "--min-carriers", "5"])
    summ = _summary(out)
    assert summ["n_islands"] == "2"
    assert summ["n_island_families"] == "5"            # F2 F3 F4 F5 F7
    assert summ["n_starship_all_families"] == "3"      # F3 F4 (E1) + F5 (E2, g6)
    assert summ["n_starship_high_families"] == "2"     # F3 F4
    assert summ["n_starship_all_families_in_island"] == "3"
    assert summ["n_captain_families"] == "1"           # F5 only; x3 hit fails E-value
    assert summ["n_islands_with_captain_family"] == "1"
    assert summ["n_island_instances_in_annotated_strains"] == "1"   # A: F2,F3,F4
    assert summ["n_islands_overlapping_element_in_annotated_strain"] == "1"
    ships = {r["element_id"]: r for r in csv.DictReader(open(out / "per_starship.tsv"), delimiter="\t")}
    assert ships["E1"]["n_genes"] == "2"
    assert ships["E1"]["n_families_in_island"] == "2"
    assert ships["E1"]["n_island_instances_overlapping"] == "1"
    assert ships["E1"]["n_genes_inside_island_instances"] == "2"
    assert ships["E2"]["n_island_instances_overlapping"] == "0"   # g6 run has size 1
    assert ships["E2"]["n_islands_hit_family_level"] == "1"       # island F5,F7 via family F5
    isl = list(csv.DictReader(gzip.open(out / "per_island.tsv.gz", "rt"), delimiter="\t"))
    assert isl[0]["starship_elements_coordinate"] == "E1"
    assert isl[0]["instances_in_annotated_strains"] == "A:c1:300-800"
    assert isl[1]["n_captain_families"] == "1"
    enr = list(csv.DictReader(open(out / "enrichment.tsv"), delimiter="\t"))
    cooc = next(r for r in enr if r["annotation"] == "starship_high" and r["universe"] == "cooc_eligible")
    # cooc universe: shell/cloud with strain_count >= 5: F2 F3 F4 F7 F8 -> N=5; islands F2 F3 F4 F7 -> n=4
    assert (cooc["N_universe"], cooc["K_annotated"], cooc["n_island"], cooc["k_both"]) == ("5", "2", "4", "2")
    assert (out / "PROVENANCE.yaml").read_text().count("sha256") >= 5
    perm = next(csv.DictReader(open(out / "coord_permutation.tsv"), delimiter="\t"))
    assert perm["observed"] == "3"   # g3 g4 (E1) + g6 (E2, F5 is an island family)


def test_extra_gene_positions_add_starship_family(tmp_path):
    rd = make_run(tmp_path)
    # an unplaced annotated-strain member A|g7 of family F8, inside E1
    with open(rd / "cluster" / "tier1_cluster.tsv", "a") as fh:
        fh.write("F8\tA|g7\n")
    ann = tmp_path / "els.tsv"
    _w(ann, "element_id\tlabel\tconfidence\tShort\tcontig\tstart\tend\nE1\tShip h1\thigh\tA\tc1\t480\t820\n")
    extra = tmp_path / "extra.tsv"
    _w(extra, "Short\tprotein_id\tcontig\tstart\tend\nA\tg7\tc1\t650\t660\nA\tg3\tc1\t1\t2\n")
    out = tmp_path / "out"
    sio.main(["--run-dir", str(rd), "--starship-annotation", str(ann), "--extra-gene-positions", str(extra),
              "--out-dir", str(out), "--n-perm", "5"])
    summ = _summary(out)
    assert summ["n_extra_gene_positions_used"] == "1"     # g3 is already placed and is ignored
    assert summ["n_starship_all_families"] == "3"         # F3 F4 + F8
    assert summ["n_island_instances_in_annotated_strains"] == "1"   # unchanged by extra genes


def test_island_instances_requires_exact_set():
    freq = {"F1": {"bin": "core"}, "F2": {"bin": "shell"}, "F3": {"bin": "cloud"}, "F4": {"bin": "shell"}}
    order = {"A": [("c", "F1", 1, 2), ("c", "F2", 3, 4), ("c", "F3", 5, 6), ("c", "F4", 7, 8)]}
    assert sio.island_instances(order, freq, {frozenset({"F2", "F3"}): 0}) == []
    got = sio.island_instances(order, freq, {frozenset({"F2", "F3", "F4"}): 7})
    assert got == [{"Short": "A", "contig": "c", "start": 3, "end": 8, "island_idx": 7, "n_genes": 3}]


def test_contig_break_splits_runs():
    freq = {"F2": {"bin": "shell"}, "F3": {"bin": "cloud"}}
    order = {"A": [("c1", "F2", 1, 2), ("c2", "F3", 5, 6)]}
    assert sio.island_instances(order, freq, {frozenset({"F2", "F3"}): 0}) == []


def test_hyper_row_counts():
    r = sio.hyper_row("x", "u", {"a", "b", "c", "d"}, {"a", "b"}, {"a", "c"})
    assert (r["N_universe"], r["K_annotated"], r["n_island"], r["k_both"]) == (4, 2, 2, 1)
    assert r["expected_k"] == 1.0


# ---------------------------------------------------------------- port additions
def _file_args(rd: Path) -> list[str]:
    return ["--frequency_table", str(rd / "frequency_table.tsv"),
            "--significant_islands", str(rd / "significant_islands.tsv"),
            "--cluster_tsv", str(rd / "cluster" / "tier1_cluster.tsv"),
            "--gene_positions", str(rd / "gene_positions.tsv"),
            "--family_positions", str(rd / "family_positions.tsv"),
            "--rescue_positions", str(rd / "rescue_positions.tsv")]


def test_file_flags_and_strain_column_match_run_dir(tmp_path):
    """Single-file flags (the Nextflow process form) + a 'strain' column give the same
    tables as --run_dir + a 'Short' column."""
    rd = make_run(tmp_path)
    a1, a2 = tmp_path / "a1.tsv", tmp_path / "a2.tsv"
    _w(a1, ELS)
    _w(a2, ELS.replace("\tShort\t", "\tstrain\t"))
    tbl = tmp_path / "cap.tblout"
    _w(tbl, TBL)
    o1, o2 = tmp_path / "o1", tmp_path / "o2"
    sio.main(["--run_dir", str(rd), "--starship_annotation", str(a1), "--captain_tblout", str(tbl),
              "--out_dir", str(o1), "--n_perm", "20"])
    sio.main(_file_args(rd) + ["--starship_annotation", str(a2), "--captain_tblout", str(tbl),
                               "--out_dir", str(o2), "--n_perm", "20"])
    for f in ("summary.tsv", "enrichment.tsv", "class_breakdown.tsv", "per_starship.tsv", "coord_permutation.tsv"):
        assert (o1 / f).read_text() == (o2 / f).read_text(), f
    for f in ("per_island.tsv.gz", "annotated_families.tsv.gz"):
        assert gzip.open(o1 / f, "rt").read() == gzip.open(o2 / f, "rt").read(), f


def test_zero_byte_stubs_mean_not_given(tmp_path):
    rd = make_run(tmp_path)
    (rd / "rescue_positions.tsv").write_text("")   # EMPTY_EVALUES_STUB when rescue is off
    ann = tmp_path / "els.tsv"
    _w(ann, ELS)
    stub = tmp_path / "empty_evalues.tsv"
    stub.write_text("")
    out = tmp_path / "out"
    sio.main(_file_args(rd) + ["--starship_annotation", str(ann), "--captain_tblout", str(stub),
                               "--out_dir", str(out), "--n_perm", "5"])
    summ = _summary(out)
    assert "n_captain_families" not in summ
    assert summ["n_starship_all_families"] == "3"
    checks = {r["check"]: r["value"] for r in csv.DictReader(open(out / "input_checks.tsv"), delimiter="\t")}
    assert checks["gain_loss"] == "not_computed"
    cb = list(csv.DictReader(open(out / "class_breakdown.tsv"), delimiter="\t"))
    assert cb and {r["direction"] for r in cb} == {"NA"}


def test_captain_stub_only_is_an_error(tmp_path):
    rd = make_run(tmp_path)
    stub = tmp_path / "empty_evalues.tsv"
    stub.write_text("")
    with pytest.raises(SystemExit):
        sio.main(_file_args(rd) + ["--captain_tblout", str(stub), "--out_dir", str(tmp_path / "o")])


def test_gain_loss_fills_direction(tmp_path):
    rd = make_run(tmp_path)
    ann = tmp_path / "els.tsv"
    _w(ann, ELS)
    gl = tmp_path / "gl.tsv"
    _w(gl, "family\tdirection\nF3\tgain\nF4\tloss\n")
    out = tmp_path / "out"
    sio.main(_file_args(rd) + ["--starship_annotation", str(ann), "--gain_loss", str(gl),
                               "--out_dir", str(out), "--n_perm", "5"])
    cb = list(csv.DictReader(open(out / "class_breakdown.tsv"), delimiter="\t"))
    dirs = {(r["annotation"], r["bin"], r["direction"]) for r in cb}
    assert ("starship_all", "shell", "gain") in dirs and ("starship_all", "cloud", "loss") in dirs


def test_contig_mismatch_is_counted(tmp_path):
    rd = make_run(tmp_path)
    ann = tmp_path / "els.tsv"
    _w(ann, ELS + "E3\tx\thigh\tA\tchrX\t1\t10\n")
    out = tmp_path / "out"
    sio.main(_file_args(rd) + ["--starship_annotation", str(ann), "--out_dir", str(out), "--n_perm", "5"])
    checks = {r["check"]: r for r in csv.DictReader(open(out / "input_checks.tsv"), delimiter="\t")}
    assert checks["elements_contig_not_in_gene_positions"]["value"] == "1"
    assert "E3" in checks["elements_contig_not_in_gene_positions"]["note"]


def test_gene_level_elements(tmp_path):
    """Gene-level table: membership is the listed genes; the span comes from their positions."""
    rd = make_run(tmp_path)
    ann = tmp_path / "els_genes.tsv"
    _w(ann, "strain\tgene_id\telement_id\tconfidence\nA\tg3\tE1\thigh\nA\tg4\tE1\thigh\nA\tg6\tE2\tlow\n"
            "A\tnope\tE9\thigh\n")
    out = tmp_path / "out"
    sio.main(_file_args(rd) + ["--starship_annotation", str(ann), "--out_dir", str(out), "--n_perm", "5"])
    summ = _summary(out)
    assert summ["n_starship_all_families"] == "3"
    assert summ["n_starship_high_families"] == "2"
    assert summ["n_elements"] == "2"                      # E9 has no placed gene
    ships = {r["element_id"]: r for r in csv.DictReader(open(out / "per_starship.tsv"), delimiter="\t")}
    assert (ships["E1"]["contig"], ships["E1"]["start"], ships["E1"]["end"]) == ("c1", "500", "800")
    checks = {r["check"]: r["value"] for r in csv.DictReader(open(out / "input_checks.tsv"), delimiter="\t")}
    assert checks["elements_gene_level_unplaced"] == "1"


def test_custom_id_sep(tmp_path):
    rd = make_run(tmp_path)
    cl = rd / "cluster" / "tier1_cluster.tsv"
    cl.write_text(cl.read_text().replace("|", "__"))
    ann = tmp_path / "els.tsv"
    _w(ann, ELS)
    out = tmp_path / "out"
    sio.main(_file_args(rd) + ["--starship_annotation", str(ann), "--id_sep", "__",
                               "--out_dir", str(out), "--n_perm", "5"])
    assert _summary(out)["n_starship_all_families"] == "3"
