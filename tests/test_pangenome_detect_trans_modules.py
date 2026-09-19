import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_detect_trans_modules import (  # noqa: E402
    build_graph,
    detect_modules,
    load_trans_edges,
)

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin" / "pangenome_detect_trans_modules.py"


def _write_pair_classification(path, rows):
    header = "family_a\tfamily_b\tclassification\tfdr_q\n"
    body = "".join(f"{a}\t{b}\t{cls}\t{q}\n" for a, b, cls, q in rows)
    path.write_text(header + body)


def test_load_trans_edges_filters_to_trans_classification_only(tmp_path):
    pc = tmp_path / "pair_classification.tsv"
    _write_pair_classification(pc, [
        ("famA", "famB", "trans", "1e-3"),
        ("famC", "famD", "unexplained_physical", "1e-3"),
        ("famE", "famF", "trans_unconfirmed", "1e-3"),
    ])
    edges = load_trans_edges(str(pc))
    assert [(a, b) for a, b, _ in edges] == [("famA", "famB")]


def test_load_trans_edges_weight_is_neg_log10_fdr_q(tmp_path):
    import math
    pc = tmp_path / "pair_classification.tsv"
    _write_pair_classification(pc, [("famA", "famB", "trans", "0.01")])
    edges = load_trans_edges(str(pc))
    assert len(edges) == 1
    _, _, weight = edges[0]
    assert abs(weight - (-math.log10(0.01))) < 1e-9


def test_load_trans_edges_clips_weight_for_zero_fdr_q(tmp_path):
    pc = tmp_path / "pair_classification.tsv"
    _write_pair_classification(pc, [("famA", "famB", "trans", "0.0")])
    edges = load_trans_edges(str(pc))
    assert edges[0][2] == 300.0


def test_detect_modules_separates_disconnected_components():
    # Two disconnected triangles -- Leiden must never merge families that
    # share no edge, regardless of resolution/seed.
    edges = [
        ("a1", "a2", 10.0), ("a2", "a3", 10.0), ("a1", "a3", 10.0),
        ("b1", "b2", 10.0), ("b2", "b3", 10.0), ("b1", "b3", 10.0),
    ]
    g = build_graph(edges)
    modules = detect_modules(g, resolution=1.0, seed=0)
    assert modules["a1"] == modules["a2"] == modules["a3"]
    assert modules["b1"] == modules["b2"] == modules["b3"]
    assert modules["a1"] != modules["b1"]


def test_main_writes_empty_valid_output_when_no_trans_edges(tmp_path, capsys):
    # Regression test for the real failure mode this pipeline hit: a
    # single-species ingroup study (immitis_in_posadasii_out,
    # posadasii_in_immitis_out) produces zero "trans"-classified pairs
    # (classify_pair()'s min_clades>=2 requirement), which the original
    # ported script (sys.exit(1) on empty edges) would turn into a hard
    # Nextflow pipeline failure. This must degrade gracefully instead:
    # write valid, empty-but-well-formed output and exit 0.
    pc = tmp_path / "pair_classification.tsv"
    _write_pair_classification(pc, [("famA", "famB", "trans_unconfirmed", "1e-3")])
    families_out = tmp_path / "family_modules.tsv"
    modules_out = tmp_path / "module_summary.tsv"

    result = subprocess.run(
        [sys.executable, str(BIN),
         "--pair_classification", str(pc),
         "--output_families", str(families_out),
         "--output_modules", str(modules_out)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert families_out.read_text() == "family\tmodule_id\tmodule_size\n"
    assert modules_out.read_text() == "module_id\tsize\n"


def test_main_writes_real_modules_end_to_end(tmp_path):
    pc = tmp_path / "pair_classification.tsv"
    _write_pair_classification(pc, [
        ("a1", "a2", "trans", "1e-5"),
        ("a2", "a3", "trans", "1e-5"),
        ("a1", "a3", "trans", "1e-5"),
        ("b1", "b2", "trans", "1e-5"),
    ])
    families_out = tmp_path / "family_modules.tsv"
    modules_out = tmp_path / "module_summary.tsv"

    result = subprocess.run(
        [sys.executable, str(BIN),
         "--pair_classification", str(pc),
         "--resolution", "1.0", "--seed", "0",
         "--output_families", str(families_out),
         "--output_modules", str(modules_out)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    families_text = families_out.read_text()
    assert "a1\t" in families_text and "b1\t" in families_text
    modules_text = modules_out.read_text()
    assert modules_text.startswith("module_id\tsize\n")
