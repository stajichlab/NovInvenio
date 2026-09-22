import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix
from pangenome_assembly_quality_qc import (
    parse_contig_lengths,
    assembly_stats_from_lengths,
    is_at_terminus,
    find_private_family_strain,
    spearman,
    partial_spearman,
    compute_per_strain_content,
    compute_terminus_enrichment,
    build_per_strain_table,
    build_correlation_table,
    render_report_markdown,
    main,
)

CONFIG_HEADER = "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup"


def _write_config(path, rows):
    with open(path, "w", newline="") as fh:
        fh.write(CONFIG_HEADER + "\n")
        for row in rows:
            fh.write(",".join(row) + "\n")


def test_parse_contig_lengths_two_contigs(tmp_path):
    fa = tmp_path / "s1.dna.fa"
    fa.write_text(">contig1\n" + "A" * 100 + "\n>contig2 desc\n" + "C" * 300 + "\n")
    lengths = parse_contig_lengths(fa)
    assert lengths == {"contig1": 100, "contig2": 300}


def test_assembly_stats_from_lengths_n50():
    stats = assembly_stats_from_lengths({"c1": 100, "c2": 300})
    assert stats == {"n_contigs": 2, "n50": 300, "total_length": 400}


def test_is_at_terminus_start_within_window():
    assert is_at_terminus(start=50, end=200, contig_length=10000, window_bp=100)


def test_is_at_terminus_end_within_window():
    assert is_at_terminus(start=9000, end=9950, contig_length=10000, window_bp=100)


def test_is_at_terminus_middle_of_contig_is_false():
    assert not is_at_terminus(start=5000, end=5200, contig_length=10000, window_bp=100)


def test_is_at_terminus_zero_length_contig_is_false():
    assert not is_at_terminus(start=1, end=1, contig_length=0, window_bp=100)


def test_find_private_family_strain_only_single_carrier():
    matrix = PresenceMatrix(families=["famA", "famB"], strains=["S1", "S2"])
    matrix.set_call("famA", "S1", "present")
    matrix.set_call("famB", "S1", "present")
    matrix.set_call("famB", "S2", "present")
    result = find_private_family_strain(matrix, {"famA", "famB"})
    assert result == {"famA": "S1"}


def test_spearman_perfect_positive_correlation():
    assert abs(spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1.0) < 1e-9


def test_spearman_perfect_negative_correlation():
    assert abs(spearman([1, 2, 3, 4], [40, 30, 20, 10]) - (-1.0)) < 1e-9


def test_partial_spearman_removes_confound():
    # y depends only on z; x is independent of both -> raw rho(x,y) may be
    # nonzero by chance in a tiny sample, but controlling for z should not
    # blow up (NaN) and should differ from the uncontrolled rho when z
    # fully explains y.
    x = [5, 1, 4, 2, 3]
    z = [1, 2, 3, 4, 5]
    y = [10, 20, 30, 40, 50]  # y perfectly determined by z
    partial = partial_spearman(x, y, z)
    # With y a perfect function of z, ryz = 1, so partial correlation of x,y
    # controlling for z is undefined (denominator 0) -> NaN.
    import math
    assert math.isnan(partial)


def test_compute_per_strain_content_counts_accessory_and_core_missing():
    matrix = PresenceMatrix(families=["core1", "shell1", "cloud1"], strains=["S1", "S2"])
    matrix.set_call("core1", "S1", "present")
    matrix.set_call("shell1", "S1", "present")
    matrix.set_call("core1", "S2", "absent")
    matrix.set_call("cloud1", "S2", "present")
    freq_bin = {"core1": "core", "shell1": "shell", "cloud1": "cloud"}
    content = compute_per_strain_content(matrix, freq_bin, ["S1", "S2"])
    assert content["S1"] == {"n_families_present": 2, "accessory_present": 1, "core_missing": 0}
    assert content["S2"] == {"n_families_present": 1, "accessory_present": 1, "core_missing": 1}


def test_compute_terminus_enrichment(tmp_path):
    gene_positions = tmp_path / "gene_positions.tsv"
    gene_positions.write_text(
        "Short\tprotein_id\tcontig\tstart\tend\n"
        "S1\tp1\tc1\t1\t100\n"       # at terminus (start<=window)
        "S1\tp2\tc1\t5000\t5100\n"   # middle
    )
    contig_lengths = {"S1": {"c1": 10000}}
    member_to_rep = {"S1|p1": "famA", "S1|p2": "famB"}
    private_families_by_strain = {"S1": {"famA"}}
    counts = compute_terminus_enrichment(
        gene_positions, contig_lengths, member_to_rep, private_families_by_strain,
        id_sep="|", window_bp=200,
    )
    assert counts["S1"]["n_all_proteins"] == 2
    assert counts["S1"]["n_all_at_terminus"] == 1
    assert counts["S1"]["n_private_proteins"] == 1
    assert counts["S1"]["n_private_at_terminus"] == 1


def test_build_per_strain_table_and_correlations():
    strains = ["S1", "S2", "S3"]
    assembly_stats = {
        "S1": {"n_contigs": 2000, "n50": 20000, "total_length": 27_000_000},
        "S2": {"n_contigs": 300, "n50": 300000, "total_length": 27_200_000},
        "S3": {"n_contigs": 100, "n50": 900000, "total_length": 27_100_000},
    }
    content = {
        "S1": {"n_families_present": 9000, "accessory_present": 900, "core_missing": 50},
        "S2": {"n_families_present": 8800, "accessory_present": 700, "core_missing": 10},
        "S3": {"n_families_present": 8700, "accessory_present": 600, "core_missing": 5},
    }
    terminus_counts = {s: {"n_all_proteins": 100, "n_all_at_terminus": 5,
                           "n_private_proteins": 10, "n_private_at_terminus": 3} for s in strains}
    rows = build_per_strain_table(strains, assembly_stats, content, terminus_counts)
    assert len(rows) == 3
    assert rows[0]["Short"] == "S1"
    assert rows[0]["pct_private_at_terminus"] == "30.00"

    corr_rows = build_correlation_table(rows)
    by_comparison = {r["comparison"]: r for r in corr_rows}
    # More fragmentation (S1) has more accessory content -> positive rho vs n_contigs,
    # negative rho vs n50, matching the confound this diagnostic exists to surface.
    assert float(by_comparison["accessory_present_vs_n_contigs"]["rho_raw"]) > 0
    assert float(by_comparison["accessory_present_vs_n50"]["rho_raw"]) < 0


def test_build_correlation_table_below_3_strains_is_dash():
    strains = ["S1", "S2"]
    assembly_stats = {
        "S1": {"n_contigs": 100, "n50": 100000, "total_length": 27_000_000},
        "S2": {"n_contigs": 200, "n50": 50000, "total_length": 27_000_000},
    }
    content = {
        "S1": {"n_families_present": 100, "accessory_present": 10, "core_missing": 1},
        "S2": {"n_families_present": 90, "accessory_present": 20, "core_missing": 2},
    }
    terminus_counts = {s: {"n_all_proteins": 0, "n_all_at_terminus": 0,
                           "n_private_proteins": 0, "n_private_at_terminus": 0} for s in strains}
    rows = build_per_strain_table(strains, assembly_stats, content, terminus_counts)
    corr_rows = build_correlation_table(rows)
    assert all(r["rho_raw"] == "-" for r in corr_rows)


def test_render_report_markdown_warns_above_threshold():
    correlation_rows = [
        {"comparison": "accessory_present_vs_n_contigs", "y": "accessory_present", "x": "n_contigs",
         "rho_raw": "0.53", "rho_partial_length": "0.60", "n_strains": 15},
        {"comparison": "accessory_present_vs_n50", "y": "accessory_present", "x": "n50",
         "rho_raw": "-0.53", "rho_partial_length": "-0.55", "n_strains": 15},
    ]
    per_strain_rows = [
        {"n_all_proteins": "100", "pct_all_at_terminus": "9.38",
         "n_private_proteins": "10", "pct_private_at_terminus": "31.40"},
    ]
    report, warned = render_report_markdown(correlation_rows, per_strain_rows, rho_warn_threshold=0.3, window_bp=1000)
    assert warned is True
    assert "WARNING" in report
    assert "31.40" not in report  # the aggregate stat is recomputed, not just copied verbatim
    assert "contig-terminus" in report.lower() or "Contig-terminus" in report


def test_render_report_markdown_no_warning_below_threshold():
    correlation_rows = [
        {"comparison": "accessory_present_vs_n_contigs", "y": "accessory_present", "x": "n_contigs",
         "rho_raw": "0.05", "rho_partial_length": "0.02", "n_strains": 15},
    ]
    per_strain_rows = [
        {"n_all_proteins": "100", "pct_all_at_terminus": "9.00",
         "n_private_proteins": "10", "pct_private_at_terminus": "10.00"},
    ]
    report, warned = render_report_markdown(correlation_rows, per_strain_rows, rho_warn_threshold=0.3, window_bp=1000)
    assert warned is False
    assert "WARNING" not in report


def test_main_end_to_end(tmp_path, capsys):
    data_dir = tmp_path / "data_dir"
    (data_dir / "dna").mkdir(parents=True)
    # S1: fragmented (many small contigs); S2/S3: well-assembled (few large contigs).
    (data_dir / "dna" / "S1.dna.fa").write_text(
        "".join(f">c{i}\n{'A' * 100}\n" for i in range(1, 21))
    )
    (data_dir / "dna" / "S2.dna.fa").write_text(">c1\n" + "A" * 2000 + "\n")
    (data_dir / "dna" / "S3.dna.fa").write_text(">c1\n" + "A" * 2000 + "\n")

    config = tmp_path / "config.csv"
    _write_config(config, [
        ("IN", "sp1", "st1", "S1.protein.fa", "S1.dna.fa", "", "S1", ""),
        ("IN", "sp1", "st2", "S2.protein.fa", "S2.dna.fa", "", "S2", ""),
        ("IN", "sp1", "st3", "S3.protein.fa", "S3.dna.fa", "", "S3", ""),
    ])

    matrix_path = tmp_path / "matrix.tsv"
    matrix_path.write_text(
        "family\tS1\tS2\tS3\n"
        "core1\tpresent\tpresent\tpresent\n"
        "shared1\tpresent\tpresent\tabsent\n"
        "private_S1\tpresent\tabsent\tabsent\n"
    )

    freq_table = tmp_path / "frequency_table.tsv"
    freq_table.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "core1\t1.0000\t3\tcore\n"
        "shared1\t0.6667\t2\tshell\n"
        "private_S1\t0.3333\t1\tsingleton\n"
    )

    cluster_tsv = tmp_path / "cluster.tsv"
    cluster_tsv.write_text(
        "core1\tS1|p1\n"
        "core1\tS2|p1\n"
        "core1\tS3|p1\n"
        "shared1\tS1|p2\n"
        "shared1\tS2|p2\n"
        "private_S1\tS1|p3\n"
    )

    gene_positions = tmp_path / "gene_positions.tsv"
    gene_positions.write_text(
        "Short\tprotein_id\tcontig\tstart\tend\n"
        "S1\tp1\tc1\t50\t150\n"
        "S1\tp2\tc2\t50\t150\n"
        "S1\tp3\tc3\t1\t100\n"   # private family protein, sits at contig terminus
        "S2\tp1\tc1\t500\t600\n"
        "S2\tp2\tc1\t700\t800\n"
        "S3\tp1\tc1\t500\t600\n"
    )

    out_dir = tmp_path / "out"
    argv = [
        "pangenome_assembly_quality_qc.py",
        "--config", str(config),
        "--data_dir", str(data_dir),
        "--matrix", str(matrix_path),
        "--frequency_table", str(freq_table),
        "--gene_positions", str(gene_positions),
        "--cluster_tsv", str(cluster_tsv),
        "--ingroup_label", "IN",
        "--terminus_window_bp", "100",
        "--out_dir", str(out_dir),
    ]
    sys.argv = argv
    rc = main()
    assert rc == 0

    with open(out_dir / "assembly_quality_vs_content.tsv", newline="") as fh:
        rows = {row["Short"]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"S1", "S2", "S3"}
    assert rows["S1"]["n_contigs"] == "20"
    assert int(rows["S1"]["accessory_present"]) == 2   # shared1 + private_S1
    assert int(rows["S1"]["core_missing"]) == 0
    assert rows["S1"]["n_private_proteins"] == "1"
    assert rows["S1"]["pct_private_at_terminus"] == "100.00"

    with open(out_dir / "assembly_quality_correlations.tsv", newline="") as fh:
        corr_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(corr_rows) == 4

    report_text = (out_dir / "assembly_quality_report.md").read_text()
    assert "Assembly Quality vs Pangenome Content QC" in report_text

    captured = capsys.readouterr()
    assert "pangenome_assembly_quality_qc" in captured.err
