"""Decision 2 (2026-10-08): permutation_p at full precision, a BH q column, optional gating."""
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_cooccurrence import add_permutation_q, format_permutation_p  # noqa: E402
from pangenome_pair_classification import classify_pair  # noqa: E402


def test_small_p_values_are_not_rounded_to_zero():
    assert format_permutation_p(1.2e-9) == "1.2000e-09"
    assert format_permutation_p(0.0213) == "2.1300e-02"
    assert float(format_permutation_p(3e-5)) == pytest.approx(3e-5)
    assert float(format_permutation_p(3e-5)) != 0.0                 # was printed as 0.0000 before


def test_legacy_format_is_available_for_comparison():
    assert format_permutation_p(3e-5, "fixed4") == "0.0000"


def test_bh_q_over_the_exact_tested_pairs():
    rows = [{"permutation_p": 0.01}, {"permutation_p": 0.02}, {"permutation_p": 0.5}]
    add_permutation_q(rows)
    assert [round(r["permutation_q"], 4) for r in rows] == [0.03, 0.03, 0.5]


def test_bh_q_is_empty_safe():
    rows = []
    add_permutation_q(rows)
    assert rows == []


def _linked_positions():
    return {s: {"famA": [("c1", 1)], "famB": [("c1", 500)]} for s in "ABCDEF"}     # not linked


def test_trans_gate_uses_the_value_it_is_given():
    gp = _linked_positions()
    two_clades = {"x": 3, "y": 3}
    # the caller passes the p-value (correction none) or the q-value (correction bh) as `permutation_p`
    assert classify_pair("famA", "famB", 0.02, two_clades, gp, {})[0] == "trans"
    assert classify_pair("famA", "famB", 0.2, two_clades, gp, {})[0] == "trans_unconfirmed"


def test_pair_classification_module_passes_the_correction_param():
    nf = (BIN.parent / "modules" / "pangenome" / "pair_classification.nf").read_text()
    assert "--perm_correction ${params.pangenome_pair_class_perm_correction}" in nf
    cfg = (BIN.parent / "nextflow.config").read_text()
    assert "pangenome_pair_class_perm_correction" in cfg and "'none'" in cfg.split("pangenome_pair_class_perm_correction")[1][:80]


def test_cooccurrence_module_changes_its_command_so_a_resumed_run_recomputes():
    # A Nextflow task hash covers the command text, not bin/ script contents. COOCCURRENCE's
    # output gains a column, so its command must change or -resume would reuse the old output.
    nf = (BIN.parent / "modules" / "pangenome" / "frequency_cooccurrence.nf").read_text()
    assert "--permutation_p_format ${params.pangenome_permutation_p_format}" in nf
