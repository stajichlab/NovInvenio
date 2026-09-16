import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_build_islands import accessory_islands, build_pair_index, find_significant_islands


def test_accessory_islands_merges_consecutive_noncore_same_contig():
    gene_order = [
        ("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1), ("g3", "chr1", 2, 2),
        ("g4", "chr1", 3, 3), ("g5", "chr1", 4, 4),
    ]
    is_core = {"g1": True, "g2": False, "g3": False, "g4": False, "g5": True}
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 1
    assert [g[0] for g in islands[0]] == ["g2", "g3", "g4"]


def test_accessory_islands_splits_at_contig_boundary():
    gene_order = [
        ("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1),
        ("g3", "chr2", 0, 0), ("g4", "chr2", 1, 1),
    ]
    is_core = {"g1": False, "g2": False, "g3": False, "g4": False}
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 2
    assert [g[0] for g in islands[0]] == ["g1", "g2"]
    assert [g[0] for g in islands[1]] == ["g3", "g4"]


def test_accessory_islands_defaults_missing_gene_to_core():
    gene_order = [("g1", "chr1", 0, 0), ("g2", "chr1", 1, 1), ("g3", "chr1", 2, 2)]
    is_core = {"g1": False, "g3": False}  # g2 missing -> defaults to core, splits the run
    islands = accessory_islands(gene_order, is_core)
    assert len(islands) == 2
    assert [g[0] for g in islands[0]] == ["g1"]
    assert [g[0] for g in islands[1]] == ["g3"]


def test_accessory_islands_warns_above_10pct_missing():
    gene_order = [(f"g{i}", "chr1", i, i) for i in range(10)]
    is_core = {"g0": False}  # 9/10 missing -> 90% > 10% threshold
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        accessory_islands(gene_order, is_core)
    assert any("absent from is_core" in str(w.message) for w in caught)


def test_build_pair_index_indexes_by_both_families():
    pairs = {frozenset({"famA", "famB"}): "unexplained_physical"}
    index = build_pair_index(pairs)
    assert index["famA"] == [(frozenset({"famA", "famB"}), "unexplained_physical")]
    assert index["famB"] == [(frozenset({"famA", "famB"}), "unexplained_physical")]


def test_find_significant_islands_keeps_only_islands_with_a_supporting_pair():
    strain_gene_orders = {
        "s1": [
            ("famA", "chr1", 0, 0), ("famB", "chr1", 1, 1),  # supported island
            ("core1", "chr1", 2, 2),
            ("famC", "chr1", 3, 3), ("famD", "chr1", 4, 4),  # unsupported island
        ],
    }
    is_core = {"famA": False, "famB": False, "core1": True, "famC": False, "famD": False}
    significant_pairs = {frozenset({"famA", "famB"}): "unexplained_physical"}
    islands = find_significant_islands(strain_gene_orders, is_core, significant_pairs)
    assert len(islands) == 1
    assert set(islands[0]["members"]) == {"famA", "famB"}
