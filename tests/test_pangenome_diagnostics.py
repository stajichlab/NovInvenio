import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_diagnostics import (
    Diagnostic,
    evaluate_rescue_redundancy,
    not_computed_diagnostic,
    render_banner_markdown,
    render_banner_html,
    write_diagnostics_tsv,
    read_funnel_tsv,
    main,
)


# --- rescue redundancy (issue #133 data, issue #134 surfacing) -------------


def test_evaluate_rescue_redundancy_trips_when_overlap_majority():
    # 60 of 100 rescuable (threshold-passing) cells rejected as overlapping
    # a different family's already-annotated gene -- over the 50% bar.
    funnel = {
        "rows_passed_threshold": 100,
        "rows_rejected_overlap": 60,
        "rows_rejected_short_rep": 5,
        "rows_rejected_hotspot": 5,
        "applied": 30,
    }
    d = evaluate_rescue_redundancy(funnel)
    assert d.diagnostic_id == "rescue_redundancy"
    assert d.status == "triggered"
    assert d.would_fail_strict is True
    assert "60" in d.detail and "100" in d.detail
    assert any("do not enable rescue" in opt.lower() for opt in d.prune_options)


def test_evaluate_rescue_redundancy_ok_below_threshold():
    funnel = {
        "rows_passed_threshold": 100,
        "rows_rejected_overlap": 10,
        "rows_rejected_short_rep": 5,
        "rows_rejected_hotspot": 5,
        "applied": 80,
    }
    d = evaluate_rescue_redundancy(funnel)
    assert d.status == "ok"
    assert d.would_fail_strict is False


def test_evaluate_rescue_redundancy_ok_when_no_rescuable_cells():
    # Zero rows passed threshold -- nothing to be redundant about; must not
    # divide by zero or fabricate a fraction.
    funnel = {
        "rows_passed_threshold": 0,
        "rows_rejected_overlap": 0,
        "rows_rejected_short_rep": 0,
        "rows_rejected_hotspot": 0,
        "applied": 0,
    }
    d = evaluate_rescue_redundancy(funnel)
    assert d.status == "ok"
    assert d.would_fail_strict is False


def test_evaluate_rescue_redundancy_custom_threshold():
    funnel = {
        "rows_passed_threshold": 100,
        "rows_rejected_overlap": 40,
        "rows_rejected_short_rep": 0,
        "rows_rejected_hotspot": 0,
        "applied": 60,
    }
    assert evaluate_rescue_redundancy(funnel, threshold=0.5).status == "ok"
    assert evaluate_rescue_redundancy(funnel, threshold=0.3).status == "triggered"


# --- not-yet-computed diagnostics (issues #130, #131, gain/loss) -----------


def test_not_computed_diagnostic_never_fabricates_a_status():
    d = not_computed_diagnostic(
        "assembly_quality_confound", tracking_issue="#130",
        description="abs(rho) accessory-vs-N50 > ~0.3",
    )
    assert d.status == "not_computed"
    assert d.would_fail_strict is False
    assert "#130" in d.detail


# --- read_funnel_tsv --------------------------------------------------------


def test_read_funnel_tsv_round_trips_metrics():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "funnel.tsv"
        path.write_text(
            "metric\tvalue\n"
            "rows_passed_threshold\t100\n"
            "rows_rejected_overlap\t60\n"
            "rows_rejected_short_rep\t5\n"
            "rows_rejected_hotspot\t5\n"
            "applied\t30\n"
        )
        funnel = read_funnel_tsv(path)
        assert funnel["rows_passed_threshold"] == 100
        assert funnel["rows_rejected_overlap"] == 60


# --- rendering ---------------------------------------------------------


def test_render_banner_markdown_lists_all_diagnostics():
    diagnostics = [
        Diagnostic("rescue_redundancy", "triggered", True,
                   "60/100 (60.0%) rescuable cells overlap another family's gene",
                   ["do not enable rescue; see clustering sweep #132"]),
        not_computed_diagnostic("clade_structure_validity", "#131", "PCoA1 variance share < ~50%"),
    ]
    md = render_banner_markdown(diagnostics)
    assert md.startswith("## Pipeline diagnostics")
    assert "rescue_redundancy" in md
    assert "TRIGGERED" in md
    assert "do not enable rescue" in md
    assert "clade_structure_validity" in md
    assert "not computed" in md.lower()
    assert "#131" in md


def test_render_banner_markdown_ok_case_still_visible():
    diagnostics = [Diagnostic("rescue_redundancy", "ok", False, "10/100 (10.0%)", [])]
    md = render_banner_markdown(diagnostics)
    assert "rescue_redundancy" in md
    assert "OK" in md


def test_render_banner_html_escapes_and_lists_diagnostics():
    diagnostics = [
        Diagnostic("rescue_redundancy", "triggered", True,
                   "<script>alert(1)</script>", ["opt & <b>one</b>"]),
    ]
    html = render_banner_html(diagnostics)
    assert "<script>alert(1)</script>" not in html  # must be escaped
    assert "&lt;script&gt;" in html
    assert "rescue_redundancy" in html


# --- write_diagnostics_tsv --------------------------------------------------


def test_write_diagnostics_tsv_writes_all_columns():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "diagnostics.tsv"
        diagnostics = [
            Diagnostic("rescue_redundancy", "triggered", True, "60/100 (60.0%)",
                       ["do not enable rescue"]),
        ]
        write_diagnostics_tsv(diagnostics, path)
        lines = path.read_text().splitlines()
        assert lines[0] == "diagnostic_id\tstatus\twould_fail_strict\tdetail\tprune_options"
        assert lines[1] == "rescue_redundancy\ttriggered\tTrue\t60/100 (60.0%)\tdo not enable rescue"


# --- main() / strict mode ---------------------------------------------------


def test_main_advisory_by_default_does_not_exit(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        funnel = tmpdir / "funnel.tsv"
        funnel.write_text(
            "metric\tvalue\n"
            "rows_passed_threshold\t100\n"
            "rows_rejected_overlap\t60\n"
            "rows_rejected_short_rep\t0\n"
            "rows_rejected_hotspot\t0\n"
            "applied\t40\n"
        )
        out_dir = tmpdir / "out"
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--rescue_funnel", str(funnel),
            "--out_dir", str(out_dir),
        ])
        main()  # must not raise SystemExit
        assert (out_dir / "diagnostics.tsv").exists()
        assert (out_dir / "diagnostics_banner.md").exists()
        assert (out_dir / "diagnostics_banner.html").exists()
        captured = capsys.readouterr()
        assert "rescue_redundancy" in captured.err


def test_main_strict_mode_fails_when_a_diagnostic_triggers(monkeypatch, capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        funnel = tmpdir / "funnel.tsv"
        funnel.write_text(
            "metric\tvalue\n"
            "rows_passed_threshold\t100\n"
            "rows_rejected_overlap\t60\n"
            "rows_rejected_short_rep\t0\n"
            "rows_rejected_hotspot\t0\n"
            "applied\t40\n"
        )
        out_dir = tmpdir / "out"
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--rescue_funnel", str(funnel),
            "--out_dir", str(out_dir),
            "--pangenome_strict",
        ])
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
        # Outputs must still be written on a strict failure -- the whole point
        # is a machine-readable record of WHY it failed.
        assert (out_dir / "diagnostics.tsv").exists()
        captured = capsys.readouterr()
        assert "pangenome_strict" in captured.err.lower()


def test_main_strict_mode_succeeds_when_nothing_triggers(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        funnel = tmpdir / "funnel.tsv"
        funnel.write_text(
            "metric\tvalue\n"
            "rows_passed_threshold\t100\n"
            "rows_rejected_overlap\t10\n"
            "rows_rejected_short_rep\t0\n"
            "rows_rejected_hotspot\t0\n"
            "applied\t90\n"
        )
        out_dir = tmpdir / "out"
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--rescue_funnel", str(funnel),
            "--out_dir", str(out_dir),
            "--pangenome_strict",
        ])
        main()  # must not raise


def test_main_without_rescue_funnel_only_reports_not_computed(monkeypatch):
    # A study that ran with --pangenome_rescue_enable false has no funnel
    # file at all -- must not error, just report rescue_redundancy as
    # not_computed alongside the other not-yet-wired diagnostics.
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir) / "out"
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--out_dir", str(out_dir),
        ])
        main()
        tsv = (out_dir / "diagnostics.tsv").read_text()
        assert "rescue_redundancy\tnot_computed" in tsv
        assert "assembly_quality_confound\tnot_computed" in tsv
        assert "clade_structure_validity\tnot_computed" in tsv
        assert "gain_loss_skew\tnot_computed" in tsv
