import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

import pangenome_report_render
from pangenome_report_render import render_report_markdown, fit_heaps_law, fit_core_decay, build_presence_bool_array
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT


def test_render_report_markdown_includes_key_sections():
    md = render_report_markdown(
        counts={"core": 100, "soft_core": 10, "shell": 50, "cloud": 200, "singleton": 300},
        size_dist={2: 50, 3: 20, 10: 5},
        classification_counts_dict={"trans": 40, "unexplained_physical": 10},
        top_domains=[{"domain": "PF00001", "fisher_p": 1e-5, "fdr_q": 1e-4}],
        n_islands=75,
        heaps_fit={"kappa": 500.0, "gamma": 0.4, "r_squared": 0.95, "is_open": True},
        core_decay={"core_inf": 90.0, "tau": 20.0, "fit_ok": True},
        strain_family_counts=[8000, 8100, 8050],
    )
    assert "# Pangenome Island + Pfam Enrichment Report" in md
    assert "## Pangenome composition" in md
    assert "## Accessory islands" in md
    assert "## Pfam domain enrichment" in md
    assert "## Pangenome openness" in md
    assert "PF00001" in md
    assert "75" in md  # n_islands appears somewhere
    assert "open" in md.lower()


def test_render_report_markdown_handles_zero_enriched_domains():
    md = render_report_markdown(
        counts={"core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0},
        size_dist={},
        classification_counts_dict={},
        top_domains=[],
        n_islands=0,
        heaps_fit=None,
        core_decay=None,
        strain_family_counts=[],
    )
    assert "No significantly enriched" in md
    # classification_counts_dict is empty -> pair_classification_summary.png
    # is never drawn (plot_classification_counts is only called `if
    # classification_counts_dict:`), so the markdown must not embed it --
    # a broken image link otherwise.
    assert "pair_classification_summary.png" not in md


def test_render_report_markdown_embeds_classification_figure_when_present():
    md = render_report_markdown(
        counts={"core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0},
        size_dist={},
        classification_counts_dict={"trans": 5},
        top_domains=[],
        n_islands=0,
        heaps_fit=None,
        core_decay=None,
        strain_family_counts=[],
    )
    assert "![Classification breakdown](figures/pair_classification_summary.png)" in md


def test_plot_classification_counts_with_non_empty_dict(tmp_path):
    # Regression test for the fig.tick_params -> ax.tick_params bug
    # (todo/report-render-tick-params-bug.md): tick_params is an Axes
    # method, not a Figure method. This raises AttributeError on any real
    # run where classification_counts_dict is non-empty -- which every
    # prior smoke test happened to avoid (zero pair-classification counts
    # at test scale), so the existing tests never called this function
    # directly with real data. First real study to hit it was the
    # coccidioides_pangenome genus_vs_ureesii run (2026-09-18).
    pangenome_report_render.plot_classification_counts(
        {"trans": 5, "unexplained_physical": 2}, tmp_path,
    )
    assert (tmp_path / "figures" / "pair_classification_summary.png").exists()
    assert (tmp_path / "figures_pdf" / "pair_classification_summary.pdf").exists()


def test_render_report_markdown_includes_marker_section_when_markers_present():
    md = render_report_markdown(
        counts={"core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0},
        size_dist={}, classification_counts_dict={}, top_domains=[], n_islands=0,
        heaps_fit=None, core_decay=None, strain_family_counts=[],
        marker_rows=[{
            "marker_name": "captain", "n_islands_with_marker": "2",
            "n_islands_total": "4", "pct_islands_with_marker": "50.0",
        }],
    )
    assert "## Marker co-occurrence" in md
    assert "captain" in md
    assert "50.0%" in md


def test_render_report_markdown_skips_marker_section_when_no_markers():
    md = render_report_markdown(
        counts={"core": 1, "soft_core": 0, "shell": 0, "cloud": 0, "singleton": 0},
        size_dist={}, classification_counts_dict={}, top_domains=[], n_islands=0,
        heaps_fit=None, core_decay=None, strain_family_counts=[], marker_rows=[],
    )
    assert "## Marker co-occurrence" not in md


def test_fit_heaps_law_detects_open_pangenome():
    # A slowly-still-growing pangenome curve (gamma < 1 -> open)
    n = np.arange(1, 21)
    pan_mean = 100 * n**0.3  # gamma=0.3, kappa=100 by construction
    result = fit_heaps_law(pan_mean)
    assert result["is_open"] is True
    assert abs(result["gamma"] - 0.3) < 0.05
    assert result["r_squared"] > 0.99


def test_fit_core_decay_returns_asymptote():
    n = np.arange(1, 21).astype(float)
    core_mean = 50 + 200 * np.exp(-n / 5)  # decays toward asymptote 50
    result = fit_core_decay(core_mean)
    assert result["fit_ok"] is True
    assert abs(result["core_inf"] - 50) < 5


def test_fit_core_decay_falls_back_gracefully_on_tiny_cohort():
    # 2 strains -> curve_fit's 3-parameter model has fewer data points than
    # free parameters, which raises TypeError (not RuntimeError) -- this
    # previously crashed REPORT_RENDER, the pipeline's last step, uncaught.
    core_mean = np.array([50.0, 45.0])
    result = fit_core_decay(core_mean)
    assert result["fit_ok"] is False
    assert result["core_inf"] == core_mean[-1]
    assert np.isnan(result["tau"])


def test_fit_heaps_law_short_circuits_below_3_strains():
    # log-log linear regression degenerates below ~3 strains (2 points fit
    # a line trivially with a meaningless R^2; 1 point can't fit at all) --
    # must short-circuit to the not-available sentinel instead of attempting
    # np.polyfit.
    pan_mean = np.array([10.0, 15.0])
    result = fit_heaps_law(pan_mean)
    assert result["fit_ok"] is False
    assert np.isnan(result["gamma"])
    assert np.isnan(result["kappa"])
    assert np.isnan(result["r_squared"])


def test_domain_table_renders_pfam_url_as_link():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "pfam_accession": "PF13577.9",
                      "pfam_url": "https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/",
                      "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "[SnoaL_2](https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/)" in md


def test_domain_table_falls_back_to_bare_name_without_pfam_url_key():
    # Existing-style row with no pfam_url key at all -- must not KeyError.
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "| SnoaL_2 |" in md


def test_domain_table_falls_back_to_bare_name_when_pfam_url_is_sentinel():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "pfam_url": "-",
                      "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "| SnoaL_2 |" in md


def test_domain_table_adds_go_columns_when_present():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5",
                      "n_go_terms": "2", "go_terms": "GO:0016491;GO:0008152"}],
        0, None, None, [], None,
    )
    assert "GO terms" in md
    assert "GO:0016491;GO:0008152" in md


def test_domain_table_omits_go_columns_when_absent():
    md = render_report_markdown(
        {}, {}, {}, [{"domain": "SnoaL_2", "fisher_p": "1e-5", "fdr_q": "2e-5"}],
        0, None, None, [], None,
    )
    assert "GO terms" not in md


def test_accessory_islands_section_adds_top_islands_table_when_locus_present():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{
            "locus_id": "S1:contig1:100-400", "island_size": "5",
            "n_strains": "2", "pfam_domains": "SnoaL_2",
        }],
    )
    assert "S1:contig1:100-400" in md


def test_accessory_islands_section_omits_top_islands_table_without_locus():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{"locus_id": "-", "island_size": "5",
                                     "n_strains": "2", "pfam_domains": "-"}],
    )
    assert "Top islands" not in md


def test_top_islands_table_excludes_single_strain_islands_by_default():
    # 62% of islands in a real genus-scale run are n_strains=1 (strain-private,
    # often assembly/annotation artefacts) and they dominate the size ranking,
    # so the headline table showed nothing shared between strains.
    md = render_report_markdown(
        {}, {}, {}, [], 2, None, None, [], None,
        islands_with_domains_rows=[
            {"locus_id": "S1:c1:1-90000", "island_size": "69", "n_strains": "1",
             "locus_start": "1", "locus_end": "90000", "pfam_domains": "-"},
            {"locus_id": "S2:c2:100-400", "island_size": "5", "n_strains": "3",
             "locus_start": "100", "locus_end": "400", "pfam_domains": "SnoaL_2"},
        ],
    )
    assert "S1:c1:1-90000" not in md
    assert "S2:c2:100-400" in md


def test_top_islands_table_notes_how_many_islands_were_excluded():
    md = render_report_markdown(
        {}, {}, {}, [], 3, None, None, [], None,
        islands_with_domains_rows=[
            {"locus_id": "S1:c1:1-90000", "island_size": "69", "n_strains": "1",
             "locus_start": "1", "locus_end": "90000", "pfam_domains": "-"},
            {"locus_id": "S3:c3:1-500", "island_size": "8", "n_strains": "1",
             "locus_start": "1", "locus_end": "500", "pfam_domains": "-"},
            {"locus_id": "S2:c2:100-400", "island_size": "5", "n_strains": "3",
             "locus_start": "100", "locus_end": "400", "pfam_domains": "SnoaL_2"},
        ],
    )
    assert "2 single-strain islands excluded" in md


def test_top_islands_min_strains_threshold_is_configurable():
    rows = [
        {"locus_id": "S1:c1:1-900", "island_size": "9", "n_strains": "2",
         "locus_start": "1", "locus_end": "900", "pfam_domains": "-"},
    ]
    assert "S1:c1:1-900" in render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=rows, top_islands_min_strains=2,
    )
    assert "S1:c1:1-900" not in render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=rows, top_islands_min_strains=3,
    )


def test_top_islands_table_reports_span_in_kb():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{
            "locus_id": "UT:scaffold_106:415-51663", "island_size": "69",
            "n_strains": "2", "locus_start": "415", "locus_end": "51663",
            "pfam_domains": "-",
        }],
    )
    assert "51.2" in md


def test_top_islands_table_header_carries_units():
    md = render_report_markdown(
        {}, {}, {}, [], 1, None, None, [], None,
        islands_with_domains_rows=[{
            "locus_id": "S1:c1:100-400", "island_size": "5", "n_strains": "2",
            "locus_start": "100", "locus_end": "400", "pfam_domains": "SnoaL_2",
        }],
    )
    assert "| Locus (strain:contig:start-end) |" in md
    assert "Families (#)" in md
    assert "Span (kb)" in md
    assert "Strains (#)" in md


def test_composition_total_names_the_unit():
    md = render_report_markdown(
        {"core": 2}, {}, {}, [], 0, None, None, [], None,
    )
    assert "Total gene families: 2" in md


def test_flagged_outlier_strains_line_present():
    md = render_report_markdown(
        {}, {}, {}, [], 0, None, None, [100, 101, 99, 100, 240],
        per_strain_rows=[
            {"Short": "A", "is_outlier": "N"}, {"Short": "E", "is_outlier": "Y"},
        ],
    )
    assert "Outlier strains" in md
    assert "E" in md


def test_flagged_outlier_strains_line_absent_when_none_flagged():
    md = render_report_markdown(
        {}, {}, {}, [], 0, None, None, [100, 101],
        per_strain_rows=[{"Short": "A", "is_outlier": "N"}, {"Short": "B", "is_outlier": "N"}],
    )
    assert "Outlier strains" not in md


def test_main_passes_islands_with_domains_rows_to_render(tmp_path, monkeypatch):
    # F1 regression test: main() previously discarded the parsed
    # islands_with_domains.tsv rows (only using them to count n_islands)
    # and never passed islands_with_domains_rows through to
    # render_report_markdown, so the "Top islands" table (and locus_id)
    # could never appear in a real report.md. Build a minimal fixture with
    # a real locus_id and assert it reaches the written report.md.
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text("family\tfrequency\tbin\nfamA\t1.0\tcore\n")
    presence_matrix = tmp_path / "presence_matrix.tsv"
    presence_matrix.write_text("family\ts1\nfamA\tpresent\n")
    islands_with_domains = tmp_path / "islands_with_domains.tsv"
    islands_with_domains.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\t"
        "locus_id\tlocus_contig\tlocus_start\tlocus_end\t"
        "n_members_with_coordinates\tn_contigs_in_locus\n"
        "2\ts1\t2\tfamA,famB\t1\ttrans\tPF00001\t"
        "s1:contig1:100-400\tcontig1\t100\t400\t2\t1\n"
    )
    island_size_distribution = tmp_path / "island_size_distribution.tsv"
    island_size_distribution.write_text("island_size\tcount\n2\t1\n")
    # classification_counts.tsv is left header-only (empty dict) so this
    # test does not exercise plot_classification_counts, a pre-existing and
    # out-of-scope code path unrelated to the F1 fix under test here.
    classification_counts = tmp_path / "classification_counts.tsv"
    classification_counts.write_text("classification\tcount\n")
    island_pfam_enrichment = tmp_path / "island_pfam_enrichment.tsv"
    island_pfam_enrichment.write_text("domain\tfisher_p\tfdr_q\n")
    marker_summary = tmp_path / "marker_summary.tsv"
    marker_summary.write_text("marker_name\tn_islands_with_marker\tn_islands_total\tpct_islands_with_marker\n")
    per_strain_summary = tmp_path / "per_strain_summary.tsv"
    per_strain_summary.write_text("Short\tn_families\tis_outlier\ns1\t1\tN\n")
    out_dir = tmp_path / "out"

    argv = [
        "pangenome_report_render.py",
        "--frequency_table", str(frequency_table),
        "--presence_matrix", str(presence_matrix),
        "--islands_with_domains", str(islands_with_domains),
        "--island_size_distribution", str(island_size_distribution),
        "--classification_counts", str(classification_counts),
        "--island_pfam_enrichment", str(island_pfam_enrichment),
        "--marker_summary", str(marker_summary),
        "--per_strain_summary", str(per_strain_summary),
        "--n_permutations", "1",
        "--out_dir", str(out_dir),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_report_render.main()

    report_md = (out_dir / "report.md").read_text()
    assert "s1:contig1:100-400" in report_md
    assert "Top islands" in report_md


def test_render_report_markdown_prepends_diagnostics_banner_when_given():
    md = render_report_markdown(
        {}, {}, {}, [], 0, None, None, [],
        diagnostics_banner="## Pipeline diagnostics\n\n- **rescue_redundancy** [OK]: fine\n",
    )
    assert md.startswith("## Pipeline diagnostics")
    assert md.index("## Pipeline diagnostics") < md.index("# Pangenome Island + Pfam Enrichment Report")


def test_render_report_markdown_omits_diagnostics_banner_when_absent():
    md = render_report_markdown({}, {}, {}, [], 0, None, None, [])
    assert "Pipeline diagnostics" not in md
    assert md.startswith("# Pangenome Island + Pfam Enrichment Report")


def test_main_prepends_diagnostics_banner_file_when_given(tmp_path, monkeypatch):
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text("family\tfrequency\tbin\nfamA\t1.0\tcore\n")
    presence_matrix = tmp_path / "presence_matrix.tsv"
    presence_matrix.write_text("family\ts1\nfamA\tpresent\n")
    islands_with_domains = tmp_path / "islands_with_domains.tsv"
    islands_with_domains.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\n"
    )
    island_size_distribution = tmp_path / "island_size_distribution.tsv"
    island_size_distribution.write_text("island_size\tcount\n")
    classification_counts = tmp_path / "classification_counts.tsv"
    classification_counts.write_text("classification\tcount\n")
    island_pfam_enrichment = tmp_path / "island_pfam_enrichment.tsv"
    island_pfam_enrichment.write_text("domain\tfisher_p\tfdr_q\n")
    marker_summary = tmp_path / "marker_summary.tsv"
    marker_summary.write_text("marker_name\tn_islands_with_marker\tn_islands_total\tpct_islands_with_marker\n")
    per_strain_summary = tmp_path / "per_strain_summary.tsv"
    per_strain_summary.write_text("Short\tn_families\ts1\t1\n")
    diagnostics_banner = tmp_path / "diagnostics_banner.md"
    diagnostics_banner.write_text("## Pipeline diagnostics\n\n- **rescue_redundancy** [OK]: fine\n")
    out_dir = tmp_path / "out"

    argv = [
        "pangenome_report_render.py",
        "--frequency_table", str(frequency_table),
        "--presence_matrix", str(presence_matrix),
        "--islands_with_domains", str(islands_with_domains),
        "--island_size_distribution", str(island_size_distribution),
        "--classification_counts", str(classification_counts),
        "--island_pfam_enrichment", str(island_pfam_enrichment),
        "--marker_summary", str(marker_summary),
        "--per_strain_summary", str(per_strain_summary),
        "--n_permutations", "1",
        "--diagnostics_banner", str(diagnostics_banner),
        "--out_dir", str(out_dir),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_report_render.main()

    report_md = (out_dir / "report.md").read_text()
    assert report_md.startswith("## Pipeline diagnostics")


def _make_mixed_state_matrix() -> PresenceMatrix:
    """Small PresenceMatrix with PRESENT, GENOME_ONLY, ABSENT-by-explicit-call,
    and ABSENT-by-omission (never called at all) cells -- the four cases
    `build_presence_bool_array` must handle identically to `is_present`."""
    matrix = PresenceMatrix(families=["famA", "famB", "famC"], strains=["s1", "s2", "s3"])
    matrix.set_call("famA", "s1", PRESENT)
    matrix.set_call("famA", "s2", ABSENT)
    matrix.set_call("famB", "s1", GENOME_ONLY)
    matrix.set_call("famB", "s3", PRESENT)
    matrix.set_call("famC", "s2", ABSENT)
    # famA/s3, famB/s2, famC/s1, famC/s3 are never set at all -> default ABSENT.
    return matrix


def test_build_presence_bool_array_matches_hand_written_expected():
    matrix = _make_mixed_state_matrix()
    family_order = ["famA", "famB", "famC"]
    strain_order = ["s1", "s2", "s3"]
    expected = np.array([
        [True, False, False],   # famA: present, absent, absent(default)
        [True, False, True],    # famB: genome_only, absent(default), present
        [False, False, False],  # famC: absent(default), absent, absent(default)
    ], dtype=bool)
    result = build_presence_bool_array(matrix, family_order, strain_order)
    assert result.dtype == bool
    assert result.shape == expected.shape
    assert np.array_equal(result, expected)


def test_build_presence_bool_array_matches_naive_is_present_oracle():
    matrix = _make_mixed_state_matrix()
    # A different (non-identity) ordering, to also exercise permutation.
    family_order = ["famC", "famA", "famB"]
    strain_order = ["s3", "s1", "s2"]

    # Ground-truth oracle: the exact naive double-loop the production code
    # used to run, written fresh here rather than imported/duplicated from
    # production, calling matrix.is_present() directly.
    oracle = np.zeros((len(family_order), len(strain_order)), dtype=bool)
    for i, fam in enumerate(family_order):
        for j, strain in enumerate(strain_order):
            oracle[i, j] = matrix.is_present(fam, strain)

    result = build_presence_bool_array(matrix, family_order, strain_order)
    assert np.array_equal(result, oracle)


def test_build_presence_bool_array_ignores_calls_outside_requested_orders():
    """A (family, strain) pair present in matrix.calls but not included in
    the requested family_order/strain_order must not affect the result --
    matches calling is_present() only for the requested rows/columns."""
    matrix = _make_mixed_state_matrix()
    family_order = ["famA"]
    strain_order = ["s1"]
    result = build_presence_bool_array(matrix, family_order, strain_order)
    assert result.shape == (1, 1)
    assert result[0, 0] == True  # noqa: E712 -- numpy bool identity check reads clearer this way here
