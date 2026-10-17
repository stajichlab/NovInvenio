"""Tests for bin/pangenome_sm_overlap.py (SM clusters vs Starships and islands)."""
import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent))
import pangenome_sm_overlap as smo  # noqa: E402
from test_pangenome_starship_overlap import ELS, TBL, _w, make_run  # noqa: E402


def _args(rd: Path) -> list[str]:
    return ["--frequency_table", str(rd / "frequency_table.tsv"),
            "--significant_islands", str(rd / "significant_islands.tsv"),
            "--cluster_tsv", str(rd / "cluster" / "tier1_cluster.tsv"),
            "--gene_positions", str(rd / "gene_positions.tsv")]


def _tsv(path: Path) -> list[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# Strain A genes on c1: g1 100-200 .. g6 1100-1200 (families F1 F2 F3 F4 F6 F5).
# E1 = c1:480-820 holds g3, g4; E2 = c1:1050-1250 holds g6.
SM = ("cluster_id\tstrain\tcontig\tstart\tend\ttype\tgenes\tcore_genes\n"
      "K1\tA\tc1\t250\t850\tT1PKS\t\tg3\n"            # coordinates -> g2 g3 g4: 2 of 3 inside E1
      "K2\tA\tc1\t50\t450\tNRPS\tg1,g2,gX\tg1\n"      # explicit list, gX unplaced: 0 inside
      "K3\tA\tc1\t850\t1250\tterpene\tg5,g6\t\n"      # g6 inside E2: 1 of 2 -> both rules
      "K4\tB\tc9\t50\t650\tterpene\t\tx1\n")          # strain B has no elements -> NA


def run(tmp_path, *extra, sm=SM, els=ELS):
    rd = make_run(tmp_path)
    smf = tmp_path / "sm.tsv"
    _w(smf, sm)
    args = _args(rd) + ["--sm_clusters", str(smf), "--out_dir", str(tmp_path / "out"), "--n_perm", "50"]
    if els is not None:
        a = tmp_path / "els.tsv"
        _w(a, els)
        args += ["--starship_annotation", str(a)]
    smo.main(args + list(extra))
    return tmp_path / "out"


def test_cluster_overlap_rules(tmp_path):
    out = run(tmp_path)
    rows = {r["cluster_id"]: r for r in _tsv(out / "sm_clusters.tsv.gz")}
    assert rows["K1"]["genes"] == "g2,g3,g4"
    assert (rows["K1"]["n_genes_in_starship"], rows["K1"]["overlap_ge1gene"], rows["K1"]["overlap_ge50pct"]) == \
        ("2", "yes", "yes")
    assert rows["K1"]["starship_elements"] == "E1"
    assert (rows["K2"]["n_genes"], rows["K2"]["n_genes_in_starship"], rows["K2"]["overlap_ge1gene"]) == ("3", "0", "no")
    assert (rows["K3"]["overlap_ge1gene"], rows["K3"]["overlap_ge50pct"]) == ("yes", "yes")
    assert (rows["K4"]["overlap_ge1gene"], rows["K4"]["n_genes_in_starship"]) == ("NA", "NA")
    assert rows["K4"]["genes"] == "x1,x2,x3"


def test_by_type_and_per_element(tmp_path):
    out = run(tmp_path)
    bt = {r["type"]: r for r in _tsv(out / "sm_by_type.tsv")}
    assert [r["type"] for r in _tsv(out / "sm_by_type.tsv")] == ["terpene", "NRPS", "T1PKS", "ALL"]
    assert (bt["ALL"]["n_clusters"], bt["ALL"]["n_overlap_ge1gene"], bt["ALL"]["n_overlap_ge50pct"]) == ("4", "2", "2")
    assert bt["ALL"]["n_with_captain"] == "NA"          # no captain tblout given
    pe = {r["element_id"]: r for r in _tsv(out / "sm_starship_overlap.tsv")}
    assert list(_tsv(out / "sm_starship_overlap.tsv")[0]) == [
        "element_id", "label", "confidence", "contig", "start", "end", "n_genes_in_element",
        "n_sm_clusters_ge1gene", "n_sm_clusters_ge50pct", "sm_clusters", "n_sm_cluster_genes_in_element",
        "captain_genes_in_element"]
    assert pe["E1"]["sm_clusters"] == "K1(T1PKS;2/3 genes)"
    assert (pe["E1"]["n_genes_in_element"], pe["E1"]["n_sm_cluster_genes_in_element"]) == ("2", "2")
    assert pe["E2"]["n_sm_clusters_ge50pct"] == "1"


def test_summary_and_families(tmp_path):
    out = run(tmp_path)
    s = {r["metric"]: r["value"] for r in _tsv(out / "summary.tsv")}
    assert list(_tsv(out / "summary.tsv")[0]) == ["metric", "value", "source"]
    assert s["sm_clusters"] == "4"
    assert s["sm_clusters_testable"] == "3"
    assert s["sm_cluster_genes_unplaced"] == "1"
    assert s["sm_clusters_overlap_ge1gene"] == "2"
    # gene-level hypergeometric over strain A genes: N=6, K=3 in elements, n=6 placed SM genes, k=3
    assert (s["hypergeom_sm_cluster_genes_in_starships_N"], s["hypergeom_sm_cluster_genes_in_starships_K"],
            s["hypergeom_sm_cluster_genes_in_starships_n"], s["hypergeom_sm_cluster_genes_in_starships_k"]) == \
        ("6", "3", "6", "3")
    fams = {r["family"]: r for r in _tsv(out / "sm_core_families.tsv.gz")}
    assert set(fams) == {"F3", "F1", "F5"}           # core genes g3 (A), g1 (A), x1 (B)
    assert fams["F3"]["in_island"] == "yes" and fams["F1"]["in_island"] == "no"
    perm = {r["statistic"]: r for r in _tsv(out / "sm_permutation.tsv")}
    assert perm["clusters_overlap_ge1gene"]["observed"] == "2"
    assert perm["sm_cluster_genes_in_elements"]["observed"] == "3"
    brk = _tsv(out / "sm_family_class_direction.tsv")
    assert brk and {r["direction"] for r in brk} == {"NA"}


def test_captain_genes(tmp_path):
    tbl = tmp_path / "cap.tblout"
    _w(tbl, TBL + "A|g6 - DUF3435 PF11917.14 1e-30 100 0 1e-30 100 0 1 1 0 0 1 1 1 1 y\n")
    out = run(tmp_path, "--captain_tblout", str(tbl))
    rows = {r["cluster_id"]: r for r in _tsv(out / "sm_clusters.tsv.gz")}
    assert rows["K3"]["captain_genes"] == "g6"
    assert rows["K4"]["captain_genes"] == "x1"        # E-value 1e-50; x3 fails 1e-5
    bt = {r["type"]: r for r in _tsv(out / "sm_by_type.tsv")}
    assert bt["ALL"]["n_with_captain"] == "2"
    pe = {r["element_id"]: r for r in _tsv(out / "sm_starship_overlap.tsv")}
    assert pe["E2"]["captain_genes_in_element"] == "g6"


def test_without_elements(tmp_path):
    out = run(tmp_path, els=None)
    bt = {r["type"]: r for r in _tsv(out / "sm_by_type.tsv")}
    assert bt["ALL"]["n_overlap_ge1gene"] == "NA"
    assert not (out / "sm_starship_overlap.tsv").exists()
    assert not (out / "sm_permutation.tsv").exists()
    s = {r["metric"]: r["value"] for r in _tsv(out / "summary.tsv")}
    assert "hypergeom_sm_cluster_genes_in_starships_k" not in s
    assert s["sm_core_families"] == "3"


def test_permutation_reproducible_and_contig_lengths(tmp_path):
    out = run(tmp_path)
    first = (out / "sm_permutation.tsv").read_text()
    out2 = run(tmp_path / "again")
    assert (out2 / "sm_permutation.tsv").read_text() == first
    cl = tmp_path / "len.tsv"
    _w(cl, "strain\tcontig\tlength\nA\tc1\t100000\n")
    out3 = run(tmp_path / "len", "--contig_lengths", str(cl))
    checks = {r["check"]: r["value"] for r in _tsv(out3 / "input_checks.tsv")}
    assert checks["permutation_contig_length"] == "contig_lengths"
