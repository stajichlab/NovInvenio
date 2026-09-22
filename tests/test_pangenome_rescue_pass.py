import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from pangenome_rescue_pass import (
    parse_tblastn_hits,
    apply_rescue,
    main,
    StructuralFilter,
)


def test_parse_tblastn_hits_filters_by_identity_and_coverage():
    # outfmt6 with qcovs appended: qseqid sseqid pident length ... qcovs
    # family "famA" query, subject header encodes strain as "s2|contig1"
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "famA\ts3|contig1\t60.0\t100\t0\t0\t1\t100\t500\t600\t1e-10\t80\t95",  # too low identity
        "famA\ts4|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t50",  # too low coverage
    ]
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {("famA", "s2")}


def test_apply_rescue_only_upgrades_absent_calls():
    pm = PresenceMatrix(families=["famA"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    applied, skipped = apply_rescue(pm, {("famA", "s1"), ("famA", "s2")})
    assert pm.call("famA", "s1") == PRESENT       # unchanged, was already PRESENT
    assert pm.call("famA", "s2") == GENOME_ONLY   # upgraded from ABSENT
    assert pm.call("famA", "s3") == ABSENT         # no rescue hit, stays ABSENT
    assert applied == 1
    assert skipped == 0


def test_parse_tblastn_hits_skips_blank_and_comment_lines():
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "",
        "# This is a comment line",
        "famB\ts3|contig2\t92.0\t100\t0\t0\t1\t100\t500\t600\t1e-40\t200\t85",
    ]
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {("famA", "s2"), ("famB", "s3")}


def test_apply_rescue_skips_unrecognized_strain_or_family():
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", ABSENT)
    applied, skipped = apply_rescue(pm, {("famA", "s1"), ("famA", "s999"), ("famZ", "s1")})
    assert pm.call("famA", "s1") == GENOME_ONLY
    assert applied == 1
    assert skipped == 2


def test_apply_rescue_idempotence_with_genome_only():
    pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", GENOME_ONLY)
    applied, skipped = apply_rescue(pm, {("famA", "s1")})
    assert pm.call("famA", "s1") == GENOME_ONLY
    assert applied == 0
    assert skipped == 0


def test_main_exits_when_all_hits_are_skipped(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text("famA\ts999|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n")

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1

        captured = capsys.readouterr()
        assert "All 1 parsed tblastn hits were skipped" in captured.err


def test_main_succeeds_when_all_hits_recognized_but_nothing_applied(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", PRESENT)
        pm.set_call("famA", "s2", PRESENT)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
            "famA\ts2|contig2\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
        )

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        main()

        assert output_file.exists()
        captured = capsys.readouterr()
        assert "Rescue: 0 ABSENT→GENOME_ONLY, 0 skipped" in captured.err


def test_main_succeeds_with_zero_tblastn_tsv_args(monkeypatch, capsys):
    # Regression test: the workflow can feed a collected channel with
    # .ifEmpty([]) when there are zero ABSENT calls upstream, producing zero
    # --tblastn_tsv args. This must not be an argparse-required error.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--output", str(output_file),
        ])

        main()

        assert output_file.exists()
        result = PresenceMatrix.from_tsv(output_file)
        assert result.call("famA", "s1") == ABSENT  # unchanged, nothing to rescue with
        captured = capsys.readouterr()
        assert "Rescue: 0 ABSENT→GENOME_ONLY, 0 skipped" in captured.err


# --- Issue #126: zero-parsed-rows guard ------------------------------------
#
# The incident this guards against: 530 per-strain tblastn files were all
# present on disk (Nextflow staged them as symlinks) but every single one
# read as empty because `zstd -dc` (no `-f`) refused the symlink and the
# then-current compressed_io ignored the non-zero exit. RESCUE_PASS applied
# zero cells and the pipeline reported SUCCESS. The underlying symlink bug
# is fixed on main; these tests cover the second line of defence: RESCUE_PASS
# itself must notice and refuse to succeed silently when the files it was
# given produced zero hit rows.


def test_main_fails_when_zero_rows_parsed_across_present_but_empty_files(monkeypatch, capsys):
    # Files ARE provided (unlike the zero-args regression test above) but
    # every one of them is empty -- the exact shape of the real incident.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        empty_file_1 = tmpdir / "s1.tblastn.tsv"
        empty_file_1.write_text("")
        empty_file_2 = tmpdir / "s2.tblastn.tsv"
        empty_file_2.write_text("")  # header-only/empty: zero data rows

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(empty_file_1),
            "--tblastn_tsv", str(empty_file_2),
            "--output", str(output_file),
        ])

        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code != 0

        captured = capsys.readouterr()
        assert "zero hit rows parsed" in captured.err.lower()
        # Must name likely causes so the next person isn't debugging blind.
        assert "mis-staged" in captured.err.lower() or "symlink" in captured.err.lower()
        assert "decompress" in captured.err.lower()
        assert "--allow_zero_rescue" in captured.err
        # Output must NOT have been written on a hard failure.
        assert not output_file.exists()


def test_main_reports_funnel_counts_in_normal_case(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.set_call("famB", "s2", ABSENT)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"  # passes, applied
            "famB\ts2|contig2\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"  # passes, applied
            "famA\ts1|contig3\t50.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"  # fails identity
            "\n"  # blank line, not a parsed row
            "# comment\n"  # comment line, not a parsed row
        )

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        main()

        assert output_file.exists()
        captured = capsys.readouterr()
        # 3 real data rows parsed (blank/comment lines excluded), 2 pass threshold.
        assert "3 hit rows parsed" in captured.err
        assert "2 rows passed" in captured.err
        assert "1 input file" in captured.err
        assert "Rescue: 2 ABSENT→GENOME_ONLY, 0 skipped" in captured.err


def test_main_warns_but_succeeds_when_rows_parsed_but_nothing_applied(monkeypatch, capsys):
    # Rows parsed and pass threshold, but every hit's strain/family is
    # already PRESENT (or unrecognized) -- legitimate, but flag it loudly.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", PRESENT)
        pm.set_call("famA", "s2", PRESENT)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
        )

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        main()

        assert output_file.exists()
        captured = capsys.readouterr()
        assert "WARNING" in captured.err
        assert "zero cells" in captured.err.lower() or "0 applied" in captured.err.lower()


def test_main_allow_zero_rescue_flag_warns_and_succeeds(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        empty_file = tmpdir / "s1.tblastn.tsv"
        empty_file.write_text("")

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(empty_file),
            "--output", str(output_file),
            "--allow_zero_rescue",
        ])

        main()  # must NOT raise SystemExit

        assert output_file.exists()
        result = PresenceMatrix.from_tsv(output_file)
        assert result.call("famA", "s1") == ABSENT
        captured = capsys.readouterr()
        assert "WARNING" in captured.err
        assert "zero hit rows parsed" in captured.err.lower()


def test_main_zero_tblastn_args_does_not_trigger_zero_row_guard(monkeypatch, capsys):
    # Regression guard: the pre-existing "no --tblastn_tsv args at all"
    # legitimate case (zero ABSENT calls upstream) must NOT be treated as
    # the "files present but empty" failure case -- it should still just
    # succeed quietly (see test_main_succeeds_with_zero_tblastn_tsv_args).
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--output", str(output_file),
        ])

        main()  # must not raise SystemExit

        assert output_file.exists()
        captured = capsys.readouterr()
        assert "ERROR" not in captured.err


# --- Issue #133: structural rescue criterion -------------------------------
#
# Rescue as configured turns 77.6% of rescuable cells into the same locus
# double-counted (the hit lands on a predicted gene already assigned to a
# DIFFERENT family in that strain) and another chunk into short/repetitive
# noise. StructuralFilter rejects those three shapes; only a hit that lands
# genuinely intergenic, on a plausible-length family rep, outside a
# repeat-hotspot window, survives.


def test_structural_filter_rejects_hit_overlapping_a_different_family_gene():
    # gene_positions: strain s1, contig1, protein p1 spans 400-700
    gene_positions = [("s1", "p1", "contig1", 400, 700)]
    member_to_rep = {"s1|p1": "famB"}  # p1 belongs to famB
    family_lengths = {"famA": 200, "famB": 200}
    sf = StructuralFilter(gene_positions, member_to_rep, family_lengths)
    # famA hit at 500-600 overlaps p1 (famB) -- reject
    assert sf.overlapping_other_family("s1", "contig1", 500, 600, "famA") is True
    # famB hit at 500-600 overlaps p1, but SAME family -- not rejected on this ground
    assert sf.overlapping_other_family("s1", "contig1", 500, 600, "famB") is False
    # A hit on a different contig/strain never overlaps
    assert sf.overlapping_other_family("s1", "contig2", 500, 600, "famA") is False
    assert sf.overlapping_other_family("s2", "contig1", 500, 600, "famA") is False


def test_structural_filter_no_overlap_when_intergenic():
    gene_positions = [("s1", "p1", "contig1", 400, 700)]
    member_to_rep = {"s1|p1": "famB"}
    sf = StructuralFilter(gene_positions, member_to_rep, {})
    # Hit at 1000-1100 does not overlap the 400-700 gene.
    assert sf.overlapping_other_family("s1", "contig1", 1000, 1100, "famA") is False


def test_structural_filter_rep_too_short():
    sf = StructuralFilter([], {}, {"famA": 100, "famB": 201}, min_rep_length=150)
    assert sf.rep_too_short("famA") is True
    assert sf.rep_too_short("famB") is False
    # A family with no recorded length (rep not in the FASTA given) is never
    # rejected on this ground alone -- absence of data isn't evidence.
    assert sf.rep_too_short("famC") is False


def test_parse_tblastn_hits_structural_filter_rejects_overlap_with_other_family():
    # famA hit at 500-600 on s1|contig1 overlaps a gene there assigned to famB.
    gene_positions = [("s1", "p1", "contig1", 400, 700)]
    member_to_rep = {"s1|p1": "famB"}
    sf = StructuralFilter(gene_positions, member_to_rep, {"famA": 200})
    lines = [
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats, structural=sf)
    assert hits == set()
    assert stats["rows_rejected_overlap"] == 1


def test_parse_tblastn_hits_structural_filter_accepts_genuinely_intergenic_hit():
    gene_positions = [("s1", "p1", "contig1", 4000, 4700)]  # far from the hit
    member_to_rep = {"s1|p1": "famB"}
    sf = StructuralFilter(gene_positions, member_to_rep, {"famA": 200})
    lines = [
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats, structural=sf)
    assert hits == {("famA", "s1")}
    assert stats["rows_rejected_overlap"] == 0
    assert stats["rows_rejected_short_rep"] == 0
    assert stats["rows_rejected_hotspot"] == 0


def test_parse_tblastn_hits_structural_filter_rejects_short_rep():
    sf = StructuralFilter([], {}, {"famA": 80}, min_rep_length=150)
    lines = [
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats, structural=sf)
    assert hits == set()
    assert stats["rows_rejected_short_rep"] == 1


def test_parse_tblastn_hits_structural_filter_rejects_repeat_hotspot():
    # Three different families all hit within the same 2kb window on
    # s1|contig1, all genuinely intergenic and long enough -- a repeat
    # hotspot, not three independent real gene losses.
    sf = StructuralFilter(
        [], {}, {"famA": 200, "famB": 200, "famC": 200},
        hotspot_window_bp=2000, hotspot_min_families=3,
    )
    lines = [
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "famB\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t550\t650\t1e-50\t200\t95",
        "famC\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t900\t1000\t1e-50\t200\t95",
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats, structural=sf)
    assert hits == set()
    assert stats["rows_rejected_hotspot"] == 3


def test_parse_tblastn_hits_structural_filter_keeps_two_family_window():
    # Only two distinct families in the window -- below the >=3 hotspot
    # threshold, so both survive (assuming they pass the other criteria).
    sf = StructuralFilter(
        [], {}, {"famA": 200, "famB": 200},
        hotspot_window_bp=2000, hotspot_min_families=3,
    )
    lines = [
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        "famB\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t550\t650\t1e-50\t200\t95",
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats, structural=sf)
    assert hits == {("famA", "s1"), ("famB", "s1")}
    assert stats["rows_rejected_hotspot"] == 0


def test_main_applies_structural_filter_end_to_end(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA", "famB"], strains=["s1"])
        pm.set_call("famA", "s1", ABSENT)
        pm.set_call("famB", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        # famA's hit overlaps an existing famX gene on s1 -- rejected.
        # famB's hit is genuinely intergenic -- applied.
        tblastn_file = tmpdir / "s1.tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
            "famB\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t5000\t5100\t1e-50\t200\t95\n"
        )

        gene_positions_file = tmpdir / "gene_positions.tsv"
        gene_positions_file.write_text(
            "Short\tprotein_id\tcontig\tstart\tend\n"
            "s1\tp1\tcontig1\t400\t700\n"
        )
        cluster_tsv_file = tmpdir / "tier1_cluster.tsv"
        cluster_tsv_file.write_text("famX\ts1|p1\n")
        rep_fasta_file = tmpdir / "tier1_rep_seq.fasta"
        rep_fasta_file.write_text(">famA\n" + "M" * 200 + "\n>famB\n" + "M" * 200 + "\n")

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--gene_positions", str(gene_positions_file),
            "--cluster_tsv", str(cluster_tsv_file),
            "--rep_fasta", str(rep_fasta_file),
            "--output", str(output_file),
        ])

        main()

        result = PresenceMatrix.from_tsv(output_file)
        assert result.call("famA", "s1") == ABSENT       # rejected: overlaps famX
        assert result.call("famB", "s1") == GENOME_ONLY  # applied: intergenic
        captured = capsys.readouterr()
        assert "1 rejected: overlaps a different family's gene" in captured.err


def test_main_requires_all_three_structural_inputs_together(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        output_file = tmpdir / "output.tsv"
        gene_positions_file = tmpdir / "gene_positions.tsv"
        gene_positions_file.write_text("Short\tprotein_id\tcontig\tstart\tend\n")

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--output", str(output_file),
            "--gene_positions", str(gene_positions_file),
        ])

        with pytest.raises(SystemExit):
            main()


# --- Issue #134: machine-readable funnel-stats sidecar ---------------------
#
# The funnel counts (rows parsed, passed threshold, structural rejections,
# cells applied) previously existed only as stderr prints -- this sidecar
# gives bin/pangenome_diagnostics.py real numbers to compute the rescue-
# redundancy diagnostic from, instead of re-parsing stderr.


def test_main_writes_funnel_tsv_when_requested(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA", "famB"], strains=["s1"])
        pm.set_call("famA", "s1", ABSENT)
        pm.set_call("famB", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "s1.tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
            "famB\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t5000\t5100\t1e-50\t200\t95\n"
        )
        gene_positions_file = tmpdir / "gene_positions.tsv"
        gene_positions_file.write_text(
            "Short\tprotein_id\tcontig\tstart\tend\n"
            "s1\tp1\tcontig1\t400\t700\n"
        )
        cluster_tsv_file = tmpdir / "tier1_cluster.tsv"
        cluster_tsv_file.write_text("famX\ts1|p1\n")
        rep_fasta_file = tmpdir / "tier1_rep_seq.fasta"
        rep_fasta_file.write_text(">famA\n" + "M" * 200 + "\n>famB\n" + "M" * 200 + "\n")

        output_file = tmpdir / "output.tsv"
        funnel_file = tmpdir / "funnel.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--gene_positions", str(gene_positions_file),
            "--cluster_tsv", str(cluster_tsv_file),
            "--rep_fasta", str(rep_fasta_file),
            "--output", str(output_file),
            "--funnel_tsv", str(funnel_file),
        ])

        main()

        assert funnel_file.exists()
        rows = dict(
            line.split("\t") for line in funnel_file.read_text().splitlines()[1:]
        )
        assert rows["rows_parsed"] == "2"
        assert rows["rows_passed_threshold"] == "2"
        assert rows["rows_rejected_overlap"] == "1"
        assert rows["rows_rejected_short_rep"] == "0"
        assert rows["rows_rejected_hotspot"] == "0"
        assert rows["applied"] == "1"


def test_main_does_not_write_funnel_tsv_when_not_requested():
    # Default behavior (no --funnel_tsv) must be unchanged: no sidecar file.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)
        output_file = tmpdir / "output.tsv"

        import sys as _sys
        _sys.argv = [
            "pangenome_rescue_pass.py",
            "--matrix", str(matrix_file),
            "--output", str(output_file),
        ]
        main()
        assert not (tmpdir / "funnel.tsv").exists()


def test_parse_tblastn_hits_tracks_row_stats():
    lines = [
        "famA\ts2|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",  # parsed + passes
        "famA\ts3|contig1\t60.0\t100\t0\t0\t1\t100\t500\t600\t1e-10\t80\t95",  # parsed, fails identity
        "",  # blank: not parsed
        "# comment",  # comment: not parsed
    ]
    stats: dict = {}
    hits = parse_tblastn_hits(lines, min_pident=90.0, min_qcov=80.0, stats=stats)
    assert hits == {("famA", "s2")}
    assert stats["rows_parsed"] == 2
    assert stats["rows_passed_threshold"] == 1
