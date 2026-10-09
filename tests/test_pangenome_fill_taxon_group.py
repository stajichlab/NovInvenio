import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_fill_taxon_group import read_clade_labels, fill_taxon_group

CONFIG_HEADER = "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup"


def _write_config(path, rows):
    with open(path, "w", newline="") as fh:
        fh.write(CONFIG_HEADER + "\n")
        for row in rows:
            fh.write(",".join(row) + "\n")


def test_read_clade_labels_parses_short_and_mash_clade(tmp_path):
    p = tmp_path / "clade_assignments.tsv"
    p.write_text("Short\tmash_clade\tpcoa1\nS1\tclade_0\t0.1\nS2\tclade_1\t-0.2\n")
    labels = read_clade_labels(str(p))
    assert labels == {"S1": "clade_0", "S2": "clade_1"}


def test_fill_taxon_group_only_fills_blank_cells(tmp_path):
    config = tmp_path / "config.csv"
    _write_config(config, [
        ("IN", "sp1", "st1", "p1.fa", "d1.fa", "g1.gff3", "S1", ""),           # blank -> fill
        ("IN", "sp2", "st2", "p2.fa", "d2.fa", "g2.gff3", "S2", "DAPC_1"),      # already set -> keep
        ("OUT", "sp3", "st3", "p3.fa", "d3.fa", "g3.gff3", "S3", ""),           # no clade label available
    ])
    output = tmp_path / "config.filled.csv"
    clade_labels = {"S1": "clade_0", "S2": "clade_5"}

    n_filled, n_total = fill_taxon_group(str(config), clade_labels, str(output))

    assert n_total == 3
    assert n_filled == 1  # only S1 was blank AND had a label available

    with open(output, newline="") as fh:
        rows = {row["Short"]: row for row in csv.DictReader(fh)}
    assert rows["S1"]["TaxonGroup"] == "clade_0"
    assert rows["S2"]["TaxonGroup"] == "DAPC_1"   # untouched, not overwritten
    assert rows["S3"]["TaxonGroup"] == ""         # still blank, no label was available


def test_fill_taxon_group_preserves_column_order(tmp_path):
    config = tmp_path / "config.csv"
    _write_config(config, [("IN", "sp1", "st1", "p1.fa", "d1.fa", "g1.gff3", "S1", "")])
    output = tmp_path / "config.filled.csv"

    fill_taxon_group(str(config), {"S1": "clade_0"}, str(output))

    with open(output, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    assert header == CONFIG_HEADER.split(",")


def test_stratum_from_mash_adds_a_column_and_leaves_taxon_group_alone(tmp_path):
    config = tmp_path / "config.csv"
    _write_config(config, [
        ("IN", "sp1", "st1", "p1.fa", "d1.fa", "g1.gff3", "S1", "Coccidioides immitis"),
        ("IN", "sp1", "st2", "p2.fa", "d2.fa", "g2.gff3", "S2", "Coccidioides immitis"),
        ("OUT", "sp2", "st3", "p3.fa", "d3.fa", "g3.gff3", "S3", "Coccidioides posadasii"),
    ])
    out = tmp_path / "out.csv"
    n_filled, _ = fill_taxon_group(str(config), {"S1": "clade_0", "S2": "clade_3"}, str(out),
                                   stratum_from_mash=True)
    rows = list(csv.DictReader(open(out)))
    assert n_filled == 0                                    # TaxonGroup was already set everywhere
    assert [r["TaxonGroup"] for r in rows] == ["Coccidioides immitis"] * 2 + ["Coccidioides posadasii"]
    assert [r["Stratum"] for r in rows] == ["clade_0", "clade_3", ""]      # outgroup has no Mash clade


def test_stratum_from_mash_keeps_an_existing_stratum(tmp_path):
    config = tmp_path / "config.csv"
    config.write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup,Stratum\n"
                      "IN,sp1,st1,p1.fa,d1.fa,g1.gff3,S1,sp1,hand_curated\n"
                      "IN,sp1,st2,p2.fa,d2.fa,g2.gff3,S2,sp1,\n")
    out = tmp_path / "out.csv"
    fill_taxon_group(str(config), {"S1": "clade_0", "S2": "clade_1"}, str(out), stratum_from_mash=True)
    assert [r["Stratum"] for r in csv.DictReader(open(out))] == ["hand_curated", "clade_1"]


def test_without_the_flag_no_stratum_column_is_added(tmp_path):
    config = tmp_path / "config.csv"
    _write_config(config, [("IN", "sp1", "st1", "p1.fa", "d1.fa", "g1.gff3", "S1", "")])
    out = tmp_path / "out.csv"
    fill_taxon_group(str(config), {"S1": "clade_0"}, str(out))
    assert "Stratum" not in csv.DictReader(open(out)).fieldnames


def test_cooccurrence_stratifies_on_stratum_when_present(tmp_path):
    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from config_parser import parse_config
    cfg = tmp_path / "config.csv"
    cfg.write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup,Stratum\n"
                   "IN,sp1,a,p.fa,d.fa,,S1,sp1,clade_0\n"
                   "IN,sp1,b,p.fa,d.fa,,S2,sp1,clade_1\n"
                   "IN,sp1,c,p.fa,d.fa,,S3,sp1,\n")
    s = {x.short: x.stratum_label for x in parse_config(cfg)}
    assert s == {"S1": "clade_0", "S2": "clade_1", "S3": "sp1"}    # blank Stratum falls back to TaxonGroup
