import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

import pangenome_report_tables
from pangenome_report_tables import (
    add_island_locus,
    add_outlier_flags,
    annotate_islands_with_domains,
    island_size_distribution,
    classification_counts,
    marker_summary,
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


def test_add_island_locus_computes_span_from_gene_positions():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "2",
        "member_families": "famA,famB", "n_supporting_pairs": "1",
        "classifications": "starship_explained",
    }]
    # member_to_rep values are always Short<id_sep>protein_id-prefixed in
    # real pipeline data (bin/pangenome_cluster_backend.py:32); gene_positions
    # is keyed on the bare protein_id (bin/pangenome_build_family_positions.py).
    member_to_rep = {"S1|protS1_a": "famA", "S1|protS1_b": "famB"}
    gene_positions = {
        ("S1", "protS1_a"): {"contig": "contig1", "start": 100, "end": 200},
        ("S1", "protS1_b"): {"contig": "contig1", "start": 300, "end": 400},
    }
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")
    assert rows[0]["locus_id"] == "S1:contig1:100-400"
    assert rows[0]["locus_contig"] == "contig1"
    assert rows[0]["n_members_with_coordinates"] == 2
    assert rows[0]["n_contigs_in_locus"] == 1


def test_add_island_locus_flags_multi_contig():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "2",
        "member_families": "famA,famB", "n_supporting_pairs": "1",
        "classifications": "starship_explained",
    }]
    member_to_rep = {"S1|protS1_a": "famA", "S1|protS1_b": "famB"}
    gene_positions = {
        ("S1", "protS1_a"): {"contig": "contig1", "start": 100, "end": 200},
        ("S1", "protS1_b"): {"contig": "contig2", "start": 10, "end": 50},
    }
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")
    # Multi-contig members must not be collapsed into a fabricated
    # cross-contig span -- locus_id/locus_contig/locus_start/locus_end all
    # fall back to the sentinel, while the counts stay real.
    assert rows[0]["n_contigs_in_locus"] == 2
    assert rows[0]["n_members_with_coordinates"] == 2
    assert rows[0]["locus_id"] == "-"
    assert rows[0]["locus_contig"] == "-"
    assert rows[0]["locus_start"] == "-"
    assert rows[0]["locus_end"] == "-"


def test_add_island_locus_sentinel_on_total_failure():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "1",
        "member_families": "famA", "n_supporting_pairs": "0",
        "classifications": "unexplained_physical",
    }]
    rows = add_island_locus(islands_rows, {}, {}, id_sep="|")
    assert rows[0]["locus_id"] == "-"
    assert rows[0]["n_members_with_coordinates"] == 0


def test_add_island_locus_respects_custom_id_sep():
    islands_rows = [{
        "n_strains": "1", "example_strain": "S1", "island_size": "1",
        "member_families": "famA", "n_supporting_pairs": "0",
        "classifications": "unexplained_physical",
    }]
    member_to_rep = {"S1_protS1a": "famA"}
    gene_positions = {("S1", "protS1a"): {"contig": "contig1", "start": 5, "end": 50}}
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="_")
    assert rows[0]["locus_id"] == "S1:contig1:5-50"


def test_add_island_locus_restricts_to_example_strain():
    """Realistic-shape regression test (bugs 1+2 from the task-3 review):
    member_to_rep keys are Short<id_sep>protein_id-prefixed, member_families
    is comma-joined (always -- bin/pangenome_build_islands.py:243 hardcodes
    ','.join(entry['members']), never id_sep), and gene_positions is keyed
    on bare (strain, protein_id) tuples. Also confirms a same-family member
    belonging to a DIFFERENT strain than the island's own example_strain is
    not picked up for this island's locus."""
    islands_rows = [{
        "n_strains": "2", "example_strain": "S1", "island_size": "2",
        "member_families": "famA,famB", "n_supporting_pairs": "1",
        "classifications": "starship_explained",
    }]
    member_to_rep = {
        "S1|protS1_a": "famA", "S1|protS1_b": "famB",
        "S2|protS2_a": "famA",  # same family, other strain -- must be ignored
    }
    gene_positions = {
        ("S1", "protS1_a"): {"contig": "contig1", "start": 100, "end": 200},
        ("S1", "protS1_b"): {"contig": "contig1", "start": 300, "end": 400},
        ("S2", "protS2_a"): {"contig": "contig9", "start": 1, "end": 9},
    }
    rows = add_island_locus(islands_rows, member_to_rep, gene_positions, id_sep="|")
    assert rows[0]["locus_id"] == "S1:contig1:100-400"
    assert rows[0]["n_members_with_coordinates"] == 2
    assert rows[0]["n_contigs_in_locus"] == 1


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
    # Only 2 strains -- add_outlier_flags emits the "-" sentinel (n<3).
    assert by_strain["s1"] == {
        "Short": "s1", "n_families": 3, "core": 1, "soft_core": 0, "shell": 1, "cloud": 1, "singleton": 0,
        "singleton_z": "-", "is_outlier": "-",
    }
    assert by_strain["s2"] == {
        "Short": "s2", "n_families": 1, "core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0,
        "singleton_z": "-", "is_outlier": "-",
    }


def test_marker_summary_computes_cooccurrence_rate():
    fieldnames = [
        "n_strains", "example_strain", "island_size", "member_families",
        "n_supporting_pairs", "classifications", "has_captain", "has_sm_backbone",
    ]
    islands_rows = [
        {"has_captain": "Y", "has_sm_backbone": "N"},
        {"has_captain": "N", "has_sm_backbone": "N"},
        {"has_captain": "Y", "has_sm_backbone": "Y"},
    ]
    result = marker_summary(islands_rows, fieldnames)
    by_name = {r["marker_name"]: r for r in result}
    assert by_name["captain"] == {
        "marker_name": "captain", "n_islands_with_marker": 2,
        "n_islands_total": 3, "pct_islands_with_marker": 66.7,
    }
    assert by_name["sm_backbone"] == {
        "marker_name": "sm_backbone", "n_islands_with_marker": 1,
        "n_islands_total": 3, "pct_islands_with_marker": 33.3,
    }


def test_marker_summary_zero_marker_columns_returns_empty_list():
    fieldnames = ["n_strains", "example_strain", "island_size", "member_families"]
    assert marker_summary([{"n_strains": "1"}], fieldnames) == []


def test_main_writes_marker_summary_tsv(tmp_path, monkeypatch):
    significant_islands = tmp_path / "significant_islands.tsv"
    significant_islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\thas_captain\n"
        "2\ts1\t2\tfamA,famB\t1\ttrans\tY\n"
        "1\ts2\t2\tfamC,famD\t1\ttrans\tN\n"
    )
    island_pfam_enrichment = tmp_path / "island_pfam_enrichment.tsv"
    island_pfam_enrichment.write_text("domain\tfisher_p\tfdr_q\n")
    pair_classification = tmp_path / "pair_classification.tsv"
    pair_classification.write_text("family_a\tfamily_b\tclassification\n")
    presence_matrix = tmp_path / "presence_matrix.tsv"
    presence_matrix.write_text("family\ts1\n")
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text("family\tfrequency\tbin\n")
    domtblout = tmp_path / "test.domtblout"
    domtblout.write_text("")
    cluster_tsv = tmp_path / "cluster.tsv"
    cluster_tsv.write_text("")
    gene_positions = tmp_path / "gene_positions.tsv"
    gene_positions.write_text("Short\tprotein_id\tcontig\tstart\tend\n")
    out_dir = tmp_path / "out"

    argv = [
        "pangenome_report_tables.py",
        "--significant_islands", str(significant_islands),
        "--island_pfam_enrichment", str(island_pfam_enrichment),
        "--pair_classification", str(pair_classification),
        "--presence_matrix", str(presence_matrix),
        "--frequency_table", str(frequency_table),
        "--domtblout", str(domtblout),
        "--cluster_tsv", str(cluster_tsv),
        "--gene_positions", str(gene_positions),
        "--out_dir", str(out_dir),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_report_tables.main()

    marker_bytes = (out_dir / "marker_summary.tsv").read_bytes()
    marker_text = marker_bytes.decode()
    assert marker_text == (
        "marker_name\tn_islands_with_marker\tn_islands_total\tpct_islands_with_marker\n"
        "captain\t1\t2\t50.0\n"
    )
    assert b"\r\n" not in marker_bytes


def test_main_zero_marker_columns_writes_header_only_marker_summary(tmp_path, monkeypatch):
    significant_islands = tmp_path / "significant_islands.tsv"
    significant_islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\n"
    )
    island_pfam_enrichment = tmp_path / "island_pfam_enrichment.tsv"
    island_pfam_enrichment.write_text("domain\tfisher_p\tfdr_q\n")
    pair_classification = tmp_path / "pair_classification.tsv"
    pair_classification.write_text("family_a\tfamily_b\tclassification\n")
    presence_matrix = tmp_path / "presence_matrix.tsv"
    presence_matrix.write_text("family\ts1\n")
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text("family\tfrequency\tbin\n")
    domtblout = tmp_path / "test.domtblout"
    domtblout.write_text("")
    cluster_tsv = tmp_path / "cluster.tsv"
    cluster_tsv.write_text("")
    gene_positions = tmp_path / "gene_positions.tsv"
    gene_positions.write_text("Short\tprotein_id\tcontig\tstart\tend\n")
    out_dir = tmp_path / "out"

    argv = [
        "pangenome_report_tables.py",
        "--significant_islands", str(significant_islands),
        "--island_pfam_enrichment", str(island_pfam_enrichment),
        "--pair_classification", str(pair_classification),
        "--presence_matrix", str(presence_matrix),
        "--frequency_table", str(frequency_table),
        "--domtblout", str(domtblout),
        "--cluster_tsv", str(cluster_tsv),
        "--gene_positions", str(gene_positions),
        "--out_dir", str(out_dir),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_report_tables.main()

    marker_text = (out_dir / "marker_summary.tsv").read_text()
    assert marker_text == "marker_name\tn_islands_with_marker\tn_islands_total\tpct_islands_with_marker\n"


def test_main_zero_islands_writes_fallback_header_and_lf_endings(tmp_path, monkeypatch):
    # Header-only significant_islands.tsv -> annotate_islands_with_domains
    # produces zero rows -> islands_with_domains.tsv must still get a real
    # header (not a blank line from an empty fieldnames list), and every
    # TSV this script writes must use LF, not csv's CRLF default.
    significant_islands = tmp_path / "significant_islands.tsv"
    significant_islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\n"
    )
    island_pfam_enrichment = tmp_path / "island_pfam_enrichment.tsv"
    island_pfam_enrichment.write_text("domain\tfisher_p\tfdr_q\n")
    pair_classification = tmp_path / "pair_classification.tsv"
    pair_classification.write_text("family_a\tfamily_b\tclassification\n")
    presence_matrix = tmp_path / "presence_matrix.tsv"
    presence_matrix.write_text("family\ts1\n")
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text("family\tfrequency\tbin\n")
    domtblout = tmp_path / "test.domtblout"
    domtblout.write_text("")
    cluster_tsv = tmp_path / "cluster.tsv"
    cluster_tsv.write_text("")
    gene_positions = tmp_path / "gene_positions.tsv"
    gene_positions.write_text("Short\tprotein_id\tcontig\tstart\tend\n")
    out_dir = tmp_path / "out"

    argv = [
        "pangenome_report_tables.py",
        "--significant_islands", str(significant_islands),
        "--island_pfam_enrichment", str(island_pfam_enrichment),
        "--pair_classification", str(pair_classification),
        "--presence_matrix", str(presence_matrix),
        "--frequency_table", str(frequency_table),
        "--domtblout", str(domtblout),
        "--cluster_tsv", str(cluster_tsv),
        "--gene_positions", str(gene_positions),
        "--out_dir", str(out_dir),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_report_tables.main()

    islands_bytes = (out_dir / "islands_with_domains.tsv").read_bytes()
    assert islands_bytes == (
        b"n_strains\texample_strain\tisland_size\tmember_families\t"
        b"n_supporting_pairs\tclassifications\tpfam_domains\t"
        b"locus_id\tlocus_contig\tlocus_start\tlocus_end\t"
        b"n_members_with_coordinates\tn_contigs_in_locus\n"
    )
    assert b"\r\n" not in islands_bytes

    per_strain_bytes = (out_dir / "per_strain_summary.tsv").read_bytes()
    assert b"\r\n" not in per_strain_bytes


def test_per_strain_summary_flags_clear_outlier():
    # 5 strains, singleton counts 10,11,9,10,150 -- strain E is the outlier.
    totals = {
        "A": {"Short": "A", "n_families": 100, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 10},
        "B": {"Short": "B", "n_families": 101, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 11},
        "C": {"Short": "C", "n_families": 99, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 9},
        "D": {"Short": "D", "n_families": 100, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 10},
        "E": {"Short": "E", "n_families": 240, "core": 80, "soft_core": 5, "shell": 3, "cloud": 2, "singleton": 150},
    }
    rows = add_outlier_flags(list(totals.values()))
    by_short = {r["Short"]: r for r in rows}
    assert by_short["E"]["is_outlier"] == "Y"
    assert by_short["A"]["is_outlier"] == "N"


def test_per_strain_summary_no_outliers_all_n():
    totals = [
        {"Short": s, "singleton": v} for s, v in
        [("A", 10), ("B", 11), ("C", 9), ("D", 10), ("E", 12)]
    ]
    rows = add_outlier_flags(totals)
    assert all(r["is_outlier"] == "N" for r in rows)


def test_per_strain_summary_below_min_n_emits_sentinel():
    totals = [{"Short": "A", "singleton": 10}, {"Short": "B", "singleton": 500}]
    rows = add_outlier_flags(totals)
    assert all(r["singleton_z"] == "-" and r["is_outlier"] == "-" for r in rows)


def test_per_strain_summary_zero_mad_emits_sentinel():
    totals = [{"Short": s, "singleton": 10} for s in ("A", "B", "C", "D")]
    rows = add_outlier_flags(totals)
    assert all(r["singleton_z"] == "-" and r["is_outlier"] == "-" for r in rows)
