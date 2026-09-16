import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_report_tables import (
    annotate_islands_with_domains,
    island_size_distribution,
    classification_counts,
    per_strain_summary,
)


def test_annotate_islands_with_domains_unions_member_domains():
    islands_rows = [{"member_families": "famA,famB", "island_size": "2"}]
    family_domains = {"famA": {"PF00001"}, "famB": {"PF00002"}}
    result = annotate_islands_with_domains(islands_rows, family_domains)
    assert result[0]["pfam_domains"] == "PF00001,PF00002"


def test_annotate_islands_with_domains_no_hits_is_dash():
    islands_rows = [{"member_families": "famA", "island_size": "1"}]
    result = annotate_islands_with_domains(islands_rows, {})
    assert result[0]["pfam_domains"] == "-"


def test_island_size_distribution_counts_by_size():
    islands_rows = [{"island_size": "2"}, {"island_size": "2"}, {"island_size": "5"}]
    assert island_size_distribution(islands_rows) == {2: 2, 5: 1}


def test_classification_counts_tallies_column(tmp_path):
    pc = tmp_path / "pair_classification.tsv"
    pc.write_text(
        "family_a\tfamily_b\tclassification\n"
        "f1\tf2\ttrans\n"
        "f3\tf4\ttrans\n"
        "f5\tf6\tambiguous_linkage\n"
    )
    assert classification_counts(str(pc)) == {"trans": 2, "ambiguous_linkage": 1}


def test_per_strain_summary_counts_genes_and_bins(tmp_path):
    pm = tmp_path / "presence_matrix.tsv"
    pm.write_text(
        "family\ts1\ts2\n"
        "famA\tpresent\tpresent\n"    # core
        "famB\tpresent\tabsent\n"     # shell
        "famC\tgenome_only\tabsent\n" # cloud -- genome_only counts as present
    )
    family_bin = {"famA": "core", "famB": "shell", "famC": "cloud"}
    result = per_strain_summary(str(pm), family_bin)
    by_strain = {r["Short"]: r for r in result}
    assert by_strain["s1"] == {"Short": "s1", "n_genes": 3, "core": 1, "soft_core": 0, "shell": 1, "cloud": 1, "singleton": 0}
    assert by_strain["s2"] == {"Short": "s2", "n_genes": 1, "core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0}
