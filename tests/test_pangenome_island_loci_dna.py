"""bin/pangenome_island_loci.py, the two DNA check passes (spec section
4b): pass 1 writes the work list (--dna_targets_dir), pass 2 applies the
calls (--dna_check true --dna_calls). Uses the four-strain fixture of
tests/test_pangenome_island_loci.py: S1 and S2 carry the locus (A, B at
ranks 5-6 on c1), S3 has the flanks only, S4 is uninformative. Every gene
at rank r is at bp r*1000+1 .. r*1000+800."""
import csv
import json
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


# run_two_locus() passes --per_rank 1 (ranks-brief.md, changed 2026-09-27): the
# drawn set is the union of the top --per_rank loci under whole_annot/whole_dna/
# species/within plus a presence-order fill, so with the default --per_rank 20
# both of these two loci would join via whole_annot regardless of --top_loci --
# --per_rank 1 keeps these tests' single-winner intent.
#
# V4: the DNA check covers all candidate loci, not just the pre-check
# --top_loci winner (spec 4b "Ranking", changed 2026-09-27). Two loci A
# (locus_id "S1:cA:100-200") and B ("S1:cB:100-200"), each with 2 full
# carriers (S1, S2) and 12 empty-site candidates (S3..S14): pre-check both
# have informative_score min(12, 2) = 2, so the locus_id tie-break draws A
# with --top_loci 1. The DNA calls make A's empties DNA-present (model
# difference: informative_score -1) and B's DNA-absent (confirmed empty:
# informative_score unchanged at 2) -- B must be drawn instead.
def write_two_locus_fixture(d: Path) -> None:
    left_a, right_a = [f"LA{i}" for i in range(1, 6)], [f"RA{i}" for i in range(1, 6)]
    left_b, right_b = [f"LB{i}" for i in range(1, 6)], [f"RB{i}" for i in range(1, 6)]
    carriers = ["S1", "S2"]
    empties = [f"S{i}" for i in range(3, 15)]
    full_a, full_b = left_a + ["A", "B"] + right_a, left_b + ["C", "D"] + right_b
    empty_a, empty_b = left_a + right_a, left_b + right_b

    (d / "islands.tsv").write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\tlocus_id\tlocus_contig\t"
        "locus_start\tlocus_end\tpfam_domains\n"
        "2\tS1\t2\tA,B\tS1:cA:100-200\tcA\t100\t200\t-\n"
        "2\tS1\t2\tC,D\tS1:cB:100-200\tcB\t100\t200\t-\n")

    pos_rows = ["Short\tfamily\tcontig\trank"]
    for s in carriers:
        pos_rows += [f"{s}\t{f}\tcA\t{r}" for r, f in enumerate(full_a)]
        pos_rows += [f"{s}\t{f}\tcB\t{r}" for r, f in enumerate(full_b)]
    for s in empties:
        pos_rows += [f"{s}\t{f}\tcA\t{r}" for r, f in enumerate(empty_a)]
        pos_rows += [f"{s}\t{f}\tcB\t{r}" for r, f in enumerate(empty_b)]
    (d / "family_positions.tsv").write_text("\n".join(pos_rows) + "\n")

    all_fams = sorted(set(full_a + full_b))
    all_strains = carriers + empties
    mat_lines = ["family\t" + "\t".join(all_strains)]
    for f in all_fams:
        cells = ["present" if (s in carriers or f in empty_a + empty_b) else "absent"
                 for s in all_strains]
        mat_lines.append(f + "\t" + "\t".join(cells))
    (d / "matrix.tsv").write_text("\n".join(mat_lines) + "\n")

    genes = ["Short\tprotein_id\tcontig\tstart\tend"]
    cluster = []
    for s, fams, contig in ([(s, full_a, "cA") for s in carriers] +
                            [(s, full_b, "cB") for s in carriers] +
                            [(s, empty_a, "cA") for s in empties] +
                            [(s, empty_b, "cB") for s in empties]):
        for r, f in enumerate(fams):
            genes.append(f"{s}\t{s}_{f}\t{contig}\t{r * 1000 + 1}\t{r * 1000 + 800}")
            cluster.append(f"{f}\t{s}|{s}_{f}")
    (d / "genes.tsv").write_text("\n".join(genes) + "\n")
    (d / "cluster.tsv").write_text("\n".join(cluster) + "\n")


def run_two_locus(d: Path, *extra) -> dict:
    write_two_locus_fixture(d)
    out = d / "island_loci.json"
    rc = main(["--islands_with_domains", str(d / "islands.tsv"),
               "--presence_matrix", str(d / "matrix.tsv"),
               "--family_positions", str(d / "family_positions.tsv"),
               "--gene_positions", str(d / "genes.tsv"),
               "--cluster_tsv", str(d / "cluster.tsv"),
               "--project", "demo", "--output", str(out), "--top_loci", "1", "--per_rank", "1",
               *extra])
    assert rc == 0
    return json.loads(out.read_text())


def write_two_locus_calls(d: Path) -> str:
    path = d / "dna_calls.tsv"
    lines = ["locus_id\tstrain\tcol\tstatus\tcoverage"]
    for s in (f"S{i}" for i in range(3, 15)):
        for col in (5, 6):
            lines.append(f"S1:cA:100-200\t{s}\t{col}\tpresent\t1.000")
            lines.append(f"S1:cB:100-200\t{s}\t{col}\tabsent\t1.000")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_without_the_dna_check_top_loci_is_chosen_from_the_annotation_ranking(tmp_path):
    data = run_two_locus(tmp_path)
    assert [loc["locus_id"] for loc in data["loci"]] == ["S1:cA:100-200"]


def test_dna_check_can_promote_a_locus_ranked_below_top_loci(tmp_path):
    data = run_two_locus(tmp_path, "--dna_check", "true",
                         "--dna_calls", write_two_locus_calls(tmp_path))
    assert [loc["locus_id"] for loc in data["loci"]] == ["S1:cB:100-200"]
    assert data["loci"][0]["counts"]["empty"] == 12


def test_pass_1_targets_cover_every_candidate_not_just_top_loci(tmp_path):
    write_two_locus_fixture(tmp_path)
    out = tmp_path / "o.json"
    rc = main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
               "--presence_matrix", str(tmp_path / "matrix.tsv"),
               "--family_positions", str(tmp_path / "family_positions.tsv"),
               "--gene_positions", str(tmp_path / "genes.tsv"),
               "--cluster_tsv", str(tmp_path / "cluster.tsv"),
               "--project", "demo", "--output", str(out), "--top_loci", "1",
               "--dna_targets_dir", str(tmp_path / "t")])
    assert rc == 0
    locus_ids = {r["locus_id"] for r in rows_of(tmp_path / "t" / "batch_001.tsv")}
    assert locus_ids == {"S1:cA:100-200", "S1:cB:100-200"}


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
