"""bin/pangenome_island_loci.py, the two DNA check passes (spec section
4b): pass 1 writes the work list (--dna_targets_dir), pass 2 applies the
calls (--dna_check true --dna_calls). Uses the four-strain fixture of
tests/test_pangenome_island_loci.py: S1 and S2 carry the locus (A, B at
ranks 5-6 on c1), S3 has the flanks only, S4 is uninformative. Every gene
at rank r is at bp r*1000+1 .. r*1000+800."""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from pangenome_island_loci import main  # noqa: E402
from test_pangenome_island_loci import LEFT, RIGHT, run, write_fixture  # noqa: E402


def write_genes(d: Path) -> list[str]:
    genes = ["Short\tprotein_id\tcontig\tstart\tend"]
    cluster = []
    orders = {"S1": LEFT + ["A", "B"] + RIGHT, "S2": LEFT + ["A", "B"] + RIGHT,
              "S3": LEFT + RIGHT}
    for s, fams in orders.items():
        for r, f in enumerate(fams):
            genes.append(f"{s}\t{s}_{f}\tc1\t{r * 1000 + 1}\t{r * 1000 + 800}")
            cluster.append(f"{f}\t{s}|{s}_{f}")
    genes.append("S4\tS4_F1\tc9\t1\t800")
    cluster.append("F1\tS4|S4_F1")
    (d / "genes.tsv").write_text("\n".join(genes) + "\n")
    (d / "cluster.tsv").write_text("\n".join(cluster) + "\n")
    return ["--gene_positions", str(d / "genes.tsv"), "--cluster_tsv", str(d / "cluster.tsv")]


def rows_of(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_pass_1_writes_the_query_and_one_target_per_checked_strain(tmp_path):
    write_fixture(tmp_path)
    run(tmp_path, *write_genes(tmp_path), "--dna_targets_dir", str(tmp_path / "t"))
    rows = rows_of(tmp_path / "t" / "batch_001.tsv")
    assert [(r["role"], r["strain"]) for r in rows] == [("query", "S2"), ("target", "S3")]
    q, t = rows
    assert (q["contig"], q["start"], q["end"], q["genes"]) == ("c1", "5001", "6800",
                                                               "5:5001-5800;6:6001-6800")
    assert (t["contig"], t["start"], t["end"]) == ("c1", "4001", "5800")


def test_pass_1_batches_loci(tmp_path):
    write_fixture(tmp_path)
    run(tmp_path, *write_genes(tmp_path), "--dna_targets_dir", str(tmp_path / "t"),
        "--dna_batch", "1")
    assert sorted(p.name for p in (tmp_path / "t").iterdir()) == ["batch_001.tsv"]


def test_pass_1_needs_gene_positions(tmp_path):
    write_fixture(tmp_path)
    with pytest.raises(SystemExit, match="needs --gene_positions"):
        run(tmp_path, "--dna_targets_dir", str(tmp_path / "t"))


def test_pass_1_gives_a_rescue_only_exemplar_column_its_tblastn_span(tmp_path):
    # Spec 4b, Ruling R26: the exemplar S2 has no gene model for B, only a
    # TBLASTN rescue hit. rescue_positions.tsv gives its start; S2's own
    # tblastn output gives the end of the HSP that starts there (minus strand).
    write_fixture(tmp_path)
    args = write_genes(tmp_path)
    genes = tmp_path / "genes.tsv"
    genes.write_text("".join(x for x in genes.read_text().splitlines(True)
                             if not x.startswith("S2\tS2_B\t")))
    run(tmp_path, *args, "--dna_targets_dir", str(tmp_path / "t0"))
    assert rows_of(tmp_path / "t0" / "batch_001.tsv")[0]["genes"] == "5:5001-5800"
    (tmp_path / "rescue.tsv").write_text("Short\tfamily\tcontig\tstart\nS2\tB\tc1\t6101\n")
    (tmp_path / "S2.tblastn.tsv").write_text(
        "B\tS2|c1\t99.0\t200\t0\t0\t1\t200\t6700\t6101\t1e-90\t390\t95\n")
    run(tmp_path, *args, "--dna_targets_dir", str(tmp_path / "t"),
        "--rescue_positions", str(tmp_path / "rescue.tsv"),
        "--rescue_tblastn", str(tmp_path / "S2.tblastn.tsv"))
    q = rows_of(tmp_path / "t" / "batch_001.tsv")[0]
    assert (q["strain"], q["start"], q["end"], q["genes"]) == (
        "S2", "5001", "6700", "5:5001-5800;6:6101-6700")


def write_calls(d: Path, status: str) -> str:
    path = d / "dna_calls_batch_001.tsv"
    path.write_text("locus_id\tstrain\tcol\tstatus\tcoverage\n"
                    f"S1:c1:100-200\tS3\t5\t{status}\t1.000\n"
                    f"S1:c1:100-200\tS3\t6\t{status}\t1.000\n")
    return str(path)


def test_pass_2_dna_present_makes_the_empty_site_a_model_difference(tmp_path):
    write_fixture(tmp_path)
    data = run(tmp_path, "--dna_check", "true", "--dna_calls", write_calls(tmp_path, "present"))
    (locus,) = data["loci"]
    assert locus["counts"] == {"full": 2, "partial": 0, "empty": 0, "model_difference": 1,
                               "uninformative": 1}
    assert [r["codes"] for r in locus["rows"] if r["row_class"] == "model_difference"] == [
        "11111" + "66" + "11111"]
    assert locus["dna"]["empty_to_model_difference"] == 1
    assert data["locus_params"]["dna_check"] is True
    assert data["locus_params"]["dna_min_id"] == 90.0


def test_pass_2_dna_absent_confirms_the_empty_site(tmp_path):
    write_fixture(tmp_path)
    (locus,) = run(tmp_path, "--dna_check", "true",
                   "--dna_calls", write_calls(tmp_path, "absent"))["loci"]
    assert locus["counts"]["empty"] == 1 and locus["dna"]["empty_confirmed"] == 1


def test_without_the_check_the_payload_says_so(tmp_path):
    data = run(tmp_path)
    assert data["locus_params"]["dna_check"] is False
    assert "dna" not in data["loci"][0]


def test_empty_calls_file_leaves_checked_strains_unchecked(tmp_path):
    write_fixture(tmp_path)
    (tmp_path / "none.tsv").write_text("")
    (locus,) = run(tmp_path, "--dna_check", "true",
                   "--dna_calls", str(tmp_path / "none.tsv"))["loci"]
    assert locus["counts"]["empty"] == 1
    assert locus["dna"] == {"checked": 0, "unchecked": 1, "empty_confirmed": 0,
                            "empty_to_model_difference": 0}


def test_cli_accepts_several_calls_files(tmp_path):
    write_fixture(tmp_path)
    first = write_calls(tmp_path, "present")
    (tmp_path / "other.tsv").write_text("locus_id\tstrain\tcol\tstatus\tcoverage\n")
    out = tmp_path / "o.json"
    assert main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
                 "--presence_matrix", str(tmp_path / "matrix.tsv"),
                 "--family_positions", str(tmp_path / "family_positions.tsv"),
                 "--dna_check", "true", "--dna_calls", first, str(tmp_path / "other.tsv"),
                 "--project", "demo", "--output", str(out)]) == 0
