import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_diagnostics import (
    Diagnostic,
    evaluate_rescue_redundancy,
    evaluate_assembly_quality_confound,
    read_correlations_tsv,
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


# --- assembly quality confound (issue #130 data, issue #134 surfacing) ------


CORR_HEADER = "comparison\ty\tx\trho_raw\trho_partial_length\tn_strains\n"


def _corr_rows(acc_contigs, acc_n50, core_contigs="0.1000", core_n50="-0.1000",
               acc_contigs_partial="0.0000", acc_n50_partial="0.0000", n=529):
    return [
        {"comparison": "accessory_present_vs_n_contigs", "y": "accessory_present", "x": "n_contigs",
         "rho_raw": acc_contigs, "rho_partial_length": acc_contigs_partial, "n_strains": str(n)},
        {"comparison": "accessory_present_vs_n50", "y": "accessory_present", "x": "n50",
         "rho_raw": acc_n50, "rho_partial_length": acc_n50_partial, "n_strains": str(n)},
        {"comparison": "core_missing_vs_n_contigs", "y": "core_missing", "x": "n_contigs",
         "rho_raw": core_contigs, "rho_partial_length": "0.0000", "n_strains": str(n)},
        {"comparison": "core_missing_vs_n50", "y": "core_missing", "x": "n50",
         "rho_raw": core_n50, "rho_partial_length": "0.0000", "n_strains": str(n)},
    ]


def test_assembly_quality_confound_trips_on_real_cocci_values():
    # The 529-strain Coccidioides values quoted in issue #130.
    d = evaluate_assembly_quality_confound(_corr_rows("0.5300", "-0.5300"))
    assert d.diagnostic_id == "assembly_quality_confound"
    assert d.status == "triggered"
    assert d.would_fail_strict is True
    assert "0.53" in d.detail and "n_contigs" in d.detail
    assert d.prune_options


def test_assembly_quality_confound_ok_below_threshold():
    d = evaluate_assembly_quality_confound(_corr_rows("0.2000", "-0.1000"))
    assert d.status == "ok"
    assert d.would_fail_strict is False


def test_assembly_quality_confound_trips_on_partial_rho():
    # Same rule as pangenome_assembly_quality_qc.py's report warning: raw OR
    # partial |rho| over the threshold.
    d = evaluate_assembly_quality_confound(
        _corr_rows("0.1000", "-0.1000", acc_n50_partial="-0.4500"))
    assert d.status == "triggered"


def test_assembly_quality_confound_ignores_core_missing_rows():
    # core_missing rows are reported, but the #130 warning rule only uses
    # the accessory_present rows.
    d = evaluate_assembly_quality_confound(
        _corr_rows("0.1000", "-0.1000", core_contigs="0.9000"))
    assert d.status == "ok"


def test_assembly_quality_confound_all_undefined_is_ok_not_fabricated():
    # < 3 strains: every rho is "-". Nothing to assess -- ok, and the detail
    # says so rather than printing a made-up number.
    d = evaluate_assembly_quality_confound(_corr_rows("-", "-", acc_contigs_partial="-",
                                                      acc_n50_partial="-", n=2))
    assert d.status == "ok"
    assert "undefined" in d.detail.lower()


def test_assembly_quality_confound_custom_threshold():
    rows = _corr_rows("0.4000", "-0.2000")
    assert evaluate_assembly_quality_confound(rows, threshold=0.5).status == "ok"
    assert evaluate_assembly_quality_confound(rows, threshold=0.3).status == "triggered"


def test_read_correlations_tsv_round_trips():
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "assembly_quality_correlations.tsv"
        p.write_text(CORR_HEADER
                     + "accessory_present_vs_n50\taccessory_present\tn50\t-0.5300\t-0.4100\t529\n")
        rows = read_correlations_tsv(p)
        assert rows == [{"comparison": "accessory_present_vs_n50", "y": "accessory_present",
                         "x": "n50", "rho_raw": "-0.5300", "rho_partial_length": "-0.4100",
                         "n_strains": "529"}]


def test_main_with_assembly_correlations_reports_real_status(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        corr = Path(tmpdir) / "assembly_quality_correlations.tsv"
        corr.write_text(CORR_HEADER
                        + "accessory_present_vs_n_contigs\taccessory_present\tn_contigs\t0.5300\t0.5000\t529\n")
        out_dir = Path(tmpdir) / "out"
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--assembly_correlations", str(corr),
            "--out_dir", str(out_dir),
        ])
        main()
        tsv = (out_dir / "diagnostics.tsv").read_text()
        assert "assembly_quality_confound\ttriggered" in tsv
        assert "assembly_quality_confound\tnot_computed" not in tsv
        assert "clade_structure_validity\tnot_computed" in tsv
        assert "gain_loss_skew\tnot_computed" in tsv


def test_main_strict_fails_on_assembly_confound(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        corr = Path(tmpdir) / "assembly_quality_correlations.tsv"
        corr.write_text(CORR_HEADER
                        + "accessory_present_vs_n50\taccessory_present\tn50\t-0.5300\t-0.5000\t529\n")
        monkeypatch.setattr(sys, "argv", [
            "pangenome_diagnostics.py",
            "--assembly_correlations", str(corr),
            "--pangenome_strict",
            "--out_dir", str(Path(tmpdir) / "out"),
        ])
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 1
