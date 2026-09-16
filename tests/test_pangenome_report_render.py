import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_report_render import render_report_markdown, fit_heaps_law, fit_core_decay
import numpy as np


def test_render_report_markdown_includes_key_sections():
    md = render_report_markdown(
        counts={"core": 100, "soft_core": 10, "shell": 50, "cloud": 200, "singleton": 300},
        size_dist={2: 50, 3: 20, 10: 5},
        classification_counts_dict={"trans": 40, "unexplained_physical": 10},
        top_domains=[{"domain": "PF00001", "fisher_p": 1e-5, "fdr_q": 1e-4}],
        n_islands=75,
        heaps_fit={"kappa": 500.0, "gamma": 0.4, "r_squared": 0.95, "is_open": True},
        core_decay={"core_inf": 90.0, "tau": 20.0, "fit_ok": True},
        strain_gene_counts=[8000, 8100, 8050],
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
        strain_gene_counts=[],
    )
    assert "No significantly enriched" in md


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
