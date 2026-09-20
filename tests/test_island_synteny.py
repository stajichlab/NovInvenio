import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_synteny import select_islands


def island(locus_id, size, n_strains):
    return {"locus_id": locus_id, "island_size": str(size),
            "n_strains": str(n_strains), "member_families": "a,b",
            "pfam_domains": "-"}


def test_single_strain_islands_are_excluded_by_default():
    # 62% of located islands on a real run are single-strain; their grid is
    # one filled row and N-1 empty ones, so there is no breakpoint to show.
    rows = [island("S1:c1:1-9000", 69, 1), island("S2:c2:1-400", 5, 3)]
    kept = select_islands(rows)
    assert [r["locus_id"] for r in kept] == ["S2:c2:1-400"]


def test_islands_are_ranked_by_size_descending():
    rows = [island("small", 3, 2), island("big", 40, 2), island("mid", 9, 2)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["big", "mid", "small"]


def test_top_n_truncates_after_filtering_not_before():
    rows = [island("private", 99, 1), island("keep1", 8, 2), island("keep2", 7, 2)]
    kept = select_islands(rows, top_n=2)
    assert [r["locus_id"] for r in kept] == ["keep1", "keep2"]


def test_min_strains_is_configurable():
    rows = [island("two", 5, 2), island("four", 4, 4)]
    assert [r["locus_id"] for r in select_islands(rows, min_strains=4)] == ["four"]


def test_rows_without_a_locus_are_excluded():
    rows = [island("-", 40, 5), island("S1:c1:1-400", 5, 5)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["S1:c1:1-400"]


def test_malformed_counts_do_not_crash():
    rows = [island("bad", "", ""), island("ok", 5, 2)]
    assert [r["locus_id"] for r in select_islands(rows)] == ["ok"]


def test_no_islands_returns_empty_list():
    assert select_islands([]) == []


from island_synteny import collapse_haplotypes


def test_identical_rows_collapse_into_one_haplotype():
    # 530 strains x a 20-56 gene island is unreadable as 530 rows, but the
    # number of DISTINCT presence patterns is typically tens.
    rows = {"s1": [True, False], "s2": [True, False], "s3": [False, True]}
    haps = collapse_haplotypes(rows)
    assert len(haps) == 2
    by_pattern = {h["pattern"]: h for h in haps}
    assert by_pattern["10"]["count"] == 2
    assert by_pattern["10"]["strains"] == ["s1", "s2"]
    assert by_pattern["01"]["count"] == 1


def test_haplotypes_sort_lexicographically_on_the_pattern():
    # Lexicographic order on the 0/1 string groups related haplotypes, which
    # is what makes breakpoints line up visually down the grid.
    rows = {"a": [True, True], "b": [False, False], "c": [False, True]}
    assert [h["pattern"] for h in collapse_haplotypes(rows)] == ["00", "01", "11"]


def test_strains_within_a_haplotype_are_sorted():
    rows = {"zeta": [True], "alpha": [True]}
    assert collapse_haplotypes(rows)[0]["strains"] == ["alpha", "zeta"]


def test_counts_sum_to_the_number_of_strains():
    rows = {f"s{i}": [i % 2 == 0] for i in range(10)}
    assert sum(h["count"] for h in collapse_haplotypes(rows)) == 10


def test_no_strains_gives_no_haplotypes():
    assert collapse_haplotypes({}) == []
