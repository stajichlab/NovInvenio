"""CLI tests for bin/pangenome_island_loci.py on a four-strain fixture."""
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
from pangenome_island_loci import main, read_n50, scan_family_positions  # noqa: E402

LEFT = ["F1", "F2", "F3", "F4", "F5"]
RIGHT = ["G1", "G2", "G3", "G4", "G5"]


def write_fixture(d: Path) -> dict:
    (d / "islands.tsv").write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\tlocus_id\tlocus_contig\t"
        "locus_start\tlocus_end\tpfam_domains\n"
        "2\tS1\t2\tA,B\tS1:c1:100-200\tc1\t100\t200\t-\n"
        "1\tS2\t1\tA\tS2:c1:100-150\tc1\t100\t150\t-\n")
    rows = ["Short\tfamily\tcontig\trank"]
    full = LEFT + ["A", "B"] + RIGHT
    for s in ("S1", "S2"):
        rows += [f"{s}\t{f}\tc1\t{r}" for r, f in enumerate(full)]
    rows += [f"S3\t{f}\tc1\t{r}" for r, f in enumerate(LEFT + RIGHT)]
    rows += ["S4\tF1\tc9\t0"]
    (d / "family_positions.tsv").write_text("\n".join(rows) + "\n")
    fams = LEFT + ["A", "B"] + RIGHT
    calls = {"S1": "present", "S2": "present", "S3": "present", "S4": "absent"}
    lines = ["family\tS1\tS2\tS3\tS4"]
    for f in fams:
        cells = []
        for s in ("S1", "S2", "S3", "S4"):
            if f in ("A", "B") and s == "S3":
                cells.append("absent")
            elif s == "S4":
                cells.append("present" if f == "F1" else "absent")
            else:
                cells.append(calls[s])
        lines.append(f + "\t" + "\t".join(cells))
    (d / "matrix.tsv").write_text("\n".join(lines) + "\n")
    (d / "aq.tsv").write_text("Short\tn_contigs\tn50\nS1\t10\t500\nS2\t10\t900\n")
    (d / "config.csv").write_text(
        "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
        "IN,Sp one,S1,S1.pep.fa,S1.dna.fa,S1.gff3,S1,t\nIN,Sp one,S2,S2.pep.fa,S2.dna.fa,S2.gff3,S2,t\n"
        "OUT,Sp two,S3,S3.pep.fa,S3.dna.fa,S3.gff3,S3,t\nOUT,Sp two,S4,S4.pep.fa,S4.dna.fa,S4.gff3,S4,t\n")
    (d / "freq.tsv").write_text("family\tfrequency\tstrain_count\tbin\n" +
                                "".join(f"{f}\t1\t4\tcore\n" for f in LEFT + RIGHT) +
                                "A\t0.5\t2\tshell\nB\t0.5\t2\tshell\n")
    return {"islands": d / "islands.tsv", "fp": d / "family_positions.tsv", "matrix": d / "matrix.tsv"}


def run(d: Path, *extra) -> dict:
    write_fixture(d)
    out = d / "island_loci.json"
    rc = main(["--islands_with_domains", str(d / "islands.tsv"),
               "--presence_matrix", str(d / "matrix.tsv"),
               "--family_positions", str(d / "family_positions.tsv"),
               "--frequency_table", str(d / "freq.tsv"),
               "--assembly_quality", str(d / "aq.tsv"),
               "--config", str(d / "config.csv"),
               "--project", "demo", "--output", str(out), *extra])
    assert rc == 0
    return json.loads(out.read_text())


def test_one_locus_with_both_islands_as_variants(tmp_path):
    data = run(tmp_path)
    assert data["n_loci_total"] == 1
    (locus,) = data["loci"]
    assert locus["key"] == "L001"
    assert locus["n_variants"] == 2
    assert locus["locus_id"] == "S1:c1:100-200"


def test_exemplar_is_the_higher_n50_carrier(tmp_path):
    (locus,) = run(tmp_path)["loci"]
    assert locus["exemplar"] == "S2"
    assert locus["tier"] == "full"
    assert locus["families"] == LEFT + ["A", "B"] + RIGHT
    assert (locus["n_left"], locus["n_locus"], locus["n_right"]) == (5, 2, 5)


def test_row_classes_and_counts(tmp_path):
    (locus,) = run(tmp_path)["loci"]
    assert locus["counts"] == {"full": 2, "partial": 0, "empty": 1, "uninformative": 1}
    by_class = {r["row_class"]: r for r in locus["rows"]}
    assert by_class["empty"]["strains"] == ["S3"]
    assert by_class["empty"]["codes"] == "11111" + "00" + "11111"
    assert by_class["full"]["count"] == 2


def test_bins_species_and_params_reach_the_payload(tmp_path):
    data = run(tmp_path)
    (locus,) = data["loci"]
    assert locus["family_bins"][5] == "shell"
    assert locus["counts_by_species"]["Sp two"]["empty"] == 1
    assert data["locus_params"]["flank"] == 5 and data["locus_params"]["k"] == 10
    assert data["n_strains_with_n50"] == 2


def test_top_loci_zero_draws_nothing_and_still_writes_json(tmp_path):
    data = run(tmp_path, "--top_loci", "0")
    assert data["loci"] == [] and data["n_loci_total"] == 1


def test_no_islands_gives_an_empty_valid_payload(tmp_path):
    write_fixture(tmp_path)
    (tmp_path / "islands.tsv").write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\tlocus_id\tlocus_contig\n")
    out = tmp_path / "o.json"
    assert main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
                 "--presence_matrix", str(tmp_path / "matrix.tsv"),
                 "--family_positions", str(tmp_path / "family_positions.tsv"),
                 "--project", "demo", "--output", str(out)]) == 0
    assert json.loads(out.read_text())["loci"] == []


def test_scan_keeps_spans_orders_and_family_copies(tmp_path):
    write_fixture(tmp_path)
    scan = scan_family_positions(str(tmp_path / "family_positions.tsv"), families={"A"},
                                 contigs={("S1", "c1")}, spans=True)
    assert scan.positions[("S1", "A")] == [("c1", 5)]
    assert scan.orders[("S1", "c1")][0] == (0, "F1")
    assert scan.spans[("S3", "c1")] == (0, 9)


def test_reads_zst_family_positions(tmp_path):
    write_fixture(tmp_path)
    subprocess.run(["zstd", "-q", "-f", str(tmp_path / "family_positions.tsv")], check=True)
    scan = scan_family_positions(str(tmp_path / "family_positions.tsv.zst"), families={"B"})
    assert scan.positions[("S2", "B")] == [("c1", 6)]


def test_missing_n50_file_is_empty(tmp_path):
    assert read_n50(str(tmp_path / "nope.tsv")) == {}


def test_cli_runs_as_a_script(tmp_path):
    write_fixture(tmp_path)
    out = tmp_path / "x.json"
    subprocess.run([sys.executable, str(BIN / "pangenome_island_loci.py"),
                    "--islands_with_domains", str(tmp_path / "islands.tsv"),
                    "--presence_matrix", str(tmp_path / "matrix.tsv"),
                    "--family_positions", str(tmp_path / "family_positions.tsv"),
                    "--project", "demo", "--output", str(out)], check=True)
    assert json.loads(out.read_text())["project"] == "demo"


def test_ranks_and_n_species_reach_the_payload(tmp_path):
    # ranks-brief.md: every drawn locus carries the five rank scores plus
    # within_species; locus_params records n_species (2 here: Sp one, Sp two).
    data = run(tmp_path)
    (locus,) = data["loci"]
    assert set(locus["ranks"]) == {"whole_annot", "whole_dna", "presence", "species", "within"}
    assert "within_species" in locus
    assert data["locus_params"]["n_species"] == 2
    assert data["locus_params"]["top_loci"] == 100 and data["locus_params"]["per_rank"] == 20
    assert data["locus_params"]["rank_by"] == "presence"


def test_rank_by_informative_is_normalised_to_presence_in_locus_params(tmp_path):
    data = run(tmp_path, "--rank_by", "informative")
    assert data["locus_params"]["rank_by"] == "presence"


def test_without_config_all_strains_are_one_unknown_species_group(tmp_path):
    # Review Focus 5: no --config means no species map; the breakpoint track
    # then has one "unknown species" group and nothing crashes.
    write_fixture(tmp_path)
    out = tmp_path / "o.json"
    assert main(["--islands_with_domains", str(tmp_path / "islands.tsv"),
                 "--presence_matrix", str(tmp_path / "matrix.tsv"),
                 "--family_positions", str(tmp_path / "family_positions.tsv"),
                 "--project", "demo", "--output", str(out)]) == 0
    (locus,) = json.loads(out.read_text())["loci"]
    assert list(locus["counts_by_species"]) == [""]
    assert all(set(bp["indel"]) <= {""} for bp in locus["breakpoints"])


def test_regions_out_lists_clinker_strains_with_anchors(tmp_path):
    import csv
    run(tmp_path, "--regions_out", str(tmp_path / "regions.tsv"))
    rows = list(csv.DictReader(open(tmp_path / "regions.tsv"), delimiter="\t"))
    assert [(r["strain"], r["reason"]) for r in rows] == [
        ("S2", "exemplar"), ("S1", "best_full"), ("S3", "best_empty")]
    s3 = rows[2]
    assert (s3["contig"], s3["rank_lo"], s3["rank_hi"]) == ("c1", "0", "9")
    assert s3["anchors"].split(";")[0] == "0:F1"


def test_clinker_max_strains_zero_selects_none(tmp_path):
    data = run(tmp_path, "--clinker_max_strains", "0", "--regions_out", str(tmp_path / "r.tsv"))
    assert data["loci"][0]["clinker_strains"] == []
    assert (tmp_path / "r.tsv").read_text().count("\n") == 1
