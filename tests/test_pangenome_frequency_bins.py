import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_frequency_bins import (
    group_class, GroupStrains, compute_group_frequency_table, resolve_group_strains,
)
from pangenome_matrix import PresenceMatrix


def _matrix(tmp_path, text):
    p = tmp_path / "m.tsv"
    p.write_text(text)
    return PresenceMatrix.from_tsv(str(p))


def test_group_class_rules_in_order():
    assert group_class(10, 10, False, False, True) == "core"
    assert group_class(1, 10, False, False, True) == "singleton"
    assert group_class(0, 10, True, True, True) == "nonrep_only"
    assert group_class(0, 10, False, True, True) == "outgroup_only"
    assert group_class(0, 10, False, True, False) == "ingroup_only"
    assert group_class(0, 10, False, False, True) == "absent"


def test_group_frequency_table_classes_both_sides(tmp_path):
    # i1-i3 ingroup reps, i4 ingroup non-rep; o1-o3 outgroup reps, o4 outgroup non-rep.
    m = _matrix(tmp_path,
        "family\ti1\ti2\ti3\ti4\to1\to2\to3\to4\n"
        "fA\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\tpresent\n"
        "fB\tabsent\tabsent\tabsent\tabsent\tpresent\tpresent\tpresent\tabsent\n"
        "fC\tabsent\tabsent\tabsent\tpresent\tabsent\tabsent\tabsent\tabsent\n"
        "fD\tabsent\tabsent\tabsent\tpresent\tpresent\tabsent\tabsent\tabsent\n"
        "fE\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tgenome_only\n"
        "fF\tpresent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], ["i4"], ["o1", "o2", "o3"], ["o4"])
    rows = {r["family"]: r for r in compute_group_frequency_table(m, g)}
    assert (rows["fA"]["bin"], rows["fA"]["bin_out"]) == ("core", "core")
    assert (rows["fB"]["bin"], rows["fB"]["bin_out"]) == ("outgroup_only", "core")
    assert (rows["fC"]["bin"], rows["fC"]["bin_out"]) == ("nonrep_only", "ingroup_only")
    # in an ingroup non-rep AND the outgroup: ingroup side is nonrep_only
    assert (rows["fD"]["bin"], rows["fD"]["bin_out"]) == ("nonrep_only", "singleton")
    assert (rows["fE"]["bin"], rows["fE"]["bin_out"]) == ("outgroup_only", "nonrep_only")
    assert (rows["fF"]["bin"], rows["fF"]["bin_out"]) == ("singleton", "ingroup_only")
    assert rows["fB"]["strain_count"] == 0 and rows["fB"]["strain_count_out"] == 3


def test_outgroup_below_minimum_is_not_binned(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\ti3\to1\to2\nfA\tpresent\tpresent\tpresent\tpresent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], [], ["o1", "o2"], [])
    (row,) = compute_group_frequency_table(m, g, outgroup_min_bin_strains=3)
    assert (row["frequency_out"], row["strain_count_out"], row["bin_out"]) == (None, None, "-")
    assert row["bin"] == "core"


def test_no_outgroup_is_not_binned(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\ti3\nfA\tabsent\tabsent\tabsent\n")
    g = GroupStrains(["i1", "i2", "i3"], [], [], [])
    (row,) = compute_group_frequency_table(m, g)
    assert row["bin_out"] == "-"
    assert row["bin"] == "absent"


def _config(tmp_path, rows):
    p = tmp_path / "c.csv"
    p.write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
                 + "".join(f"{g},sp,{s},{s}.fa,{s}.dna,{s}.gff3,{s},\n" for g, s in rows))
    return str(p)


def test_resolve_group_strains_with_inventory(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\to1\to2\tx1\nfA\tpresent\tpresent\tpresent\tpresent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("OUT", "o1"), ("OUT", "o2"), ("NEAR_INGROUP", "x1")])
    inv = tmp_path / "inv.tsv"
    inv.write_text("Short\tdedup_group\tis_representative\ni1\t0\t1\ni2\t0\t0\no1\t1\t1\no2\t2\t1\nx1\t3\t1\n")
    g = resolve_group_strains(m, cfg, str(inv), "IN", "OUT")
    assert g == GroupStrains(["i1"], ["i2"], ["o1", "o2"], [])


def test_resolve_group_strains_empty_inventory_stub_means_no_dereplication(tmp_path):
    m = _matrix(tmp_path, "family\ti1\ti2\to1\nfA\tpresent\tpresent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("OUT", "o1")])
    stub = tmp_path / "empty_evalues.tsv"
    stub.write_text("")
    g = resolve_group_strains(m, cfg, str(stub), "IN", "OUT")
    assert g == GroupStrains(["i1", "i2"], [], ["o1"], [])


import subprocess


def test_cli_writes_out_columns(tmp_path):
    (tmp_path / "m.tsv").write_text(
        "family\ti1\ti2\ti3\to1\to2\to3\nfA\tpresent\tpresent\tpresent\tpresent\tabsent\tpresent\n")
    cfg = _config(tmp_path, [("IN", "i1"), ("IN", "i2"), ("IN", "i3"),
                             ("OUT", "o1"), ("OUT", "o2"), ("OUT", "o3")])
    out = tmp_path / "ft.tsv"
    subprocess.run([sys.executable, str(Path(__file__).parent.parent / "bin" / "pangenome_frequency_bins.py"),
                    "--matrix", str(tmp_path / "m.tsv"), "--config", cfg,
                    "--outgroup_label", "OUT", "--outgroup_min_bin_strains", "3",
                    "--output", str(out)], check=True)
    lines = out.read_text().splitlines()
    assert lines[0] == "family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out"
    assert lines[1] == "fA\t1.0000\t3\tcore\t0.6667\t2\tshell"
