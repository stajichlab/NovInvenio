"""Pairwise pipeline clarity and correctness fixes found in the 2026-10-07 review."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))

from core_report_template import CORE_HTML_TEMPLATE  # noqa: E402
from losses_report_template import LOSSES_HTML_TEMPLATE  # noqa: E402
from report_template import HTML_TEMPLATE  # noqa: E402


def test_evalue_param_reaches_the_presence_matrix_filter():
    # --evalue must set filter 1; before, the script's hard-coded 1e-5 always applied.
    nf = (REPO / "modules" / "build_presence_matrix.nf").read_text()
    assert "--default-evalue ${params.evalue}" in nf


def test_every_report_explains_how_to_read_it():
    for template in (HTML_TEMPLATE, CORE_HTML_TEMPLATE, LOSSES_HTML_TEMPLATE):
        assert "How to read this page" in template
        assert "Qualifying hit" in template


def test_novelty_page_defines_the_tblastn_and_signal_columns():
    assert "does not remove a candidate" in HTML_TEMPLATE        # TBLASTN is evidence only
    assert "Outgroup signal" in HTML_TEMPLATE and "not a validated cutoff" in HTML_TEMPLATE


def test_novelty_page_says_no_qualifying_hit_not_absent():
    assert "absent from every outgroup proteome" not in HTML_TEMPLATE
    assert "No qualifying hit" in HTML_TEMPLATE


def test_core_note_states_the_count_threshold():
    assert "Math.ceil(DATA.core_min_frac * PROTEOMES.length" in CORE_HTML_TEMPLATE
    assert "of every sampled proteome" not in CORE_HTML_TEMPLATE


def test_core_family_tile_counts_only_families_core_rows_belong_to():
    assert "famUsed" in CORE_HTML_TEMPLATE
    assert 'textContent = FAMILIES.length.toLocaleString()' not in CORE_HTML_TEMPLATE


def test_make_novelties_docstring_does_not_claim_tblastn_filters_candidates():
    doc = (REPO / "bin" / "make_novelties.py").read_text().split('"""')[1]
    assert "OFF in the pipeline" in doc and "--skip_tblastn_filter" in doc


def test_core_card_lists_the_hit_in_each_species():
    assert 'field("Best hit in each species", hitsNode(row))' in CORE_HTML_TEMPLATE
    assert "DATA.protein_names" in CORE_HTML_TEMPLATE
