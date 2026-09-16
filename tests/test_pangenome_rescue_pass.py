import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from pangenome_rescue_pass import parse_tblastn_hits, apply_rescue, main


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
