import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from island_synteny import (
    select_islands,
    collapse_haplotypes,
    order_families_by_locus,
    build_payload,
)
from pangenome_matrix import PresenceMatrix


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


def test_families_are_ordered_by_rank_in_the_example_strain():
    positions = {("S1", "famC"): 3, ("S1", "famA"): 1, ("S1", "famB"): 2}
    assert order_families_by_locus(["famA", "famB", "famC"], positions, "S1") == \
        ["famA", "famB", "famC"]


def test_input_order_does_not_matter():
    positions = {("S1", "famC"): 3, ("S1", "famA"): 1, ("S1", "famB"): 2}
    assert order_families_by_locus(["famC", "famA", "famB"], positions, "S1") == \
        ["famA", "famB", "famC"]


def test_families_without_a_position_sort_last_in_stable_name_order():
    # A member with no coordinate in this strain still deserves a column --
    # dropping it would silently shrink the island.
    positions = {("S1", "famB"): 2}
    assert order_families_by_locus(["famZ", "famB", "famA"], positions, "S1") == \
        ["famB", "famA", "famZ"]


def test_positions_from_another_strain_are_ignored():
    positions = {("S2", "famA"): 1, ("S1", "famB"): 5}
    assert order_families_by_locus(["famA", "famB"], positions, "S1") == \
        ["famB", "famA"]


def test_empty_family_list_gives_empty_order():
    assert order_families_by_locus([], {}, "S1") == []


def make_matrix(calls):
    """calls: {family: {strain: bool}} -> a PresenceMatrix."""
    families = sorted(calls.keys())
    strains = sorted(set(s for per_strain in calls.values() for s in per_strain.keys()))
    m = PresenceMatrix(families=families, strains=strains)
    for family, per_strain in calls.items():
        for strain, present in per_strain.items():
            m.set_call(family, strain, "present" if present else "absent")
    return m


def test_payload_carries_one_island_with_ordered_columns_and_haplotypes():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famB,famA", "pfam_domains": "NACHT,MFS_1",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True, "S3": False},
        "famB": {"S1": True, "S2": False, "S3": False},
    })
    positions = {("S1", "famA"): 1, ("S1", "famB"): 2}
    payload = build_payload(rows, matrix, positions, project="demo")

    assert payload["project"] == "demo"
    assert len(payload["islands"]) == 1
    island = payload["islands"][0]
    assert island["families"] == ["famA", "famB"]
    assert island["locus_id"] == "S1:c1:1-400"
    patterns = {h["pattern"]: h["count"] for h in island["haplotypes"]}
    assert patterns == {"11": 1, "10": 1, "00": 1}


def test_payload_reports_how_many_islands_were_excluded():
    rows = [
        {"locus_id": "S1:c1:1-9", "island_size": "9", "n_strains": "1",
         "member_families": "famA", "pfam_domains": "-", "example_strain": "S1"},
        {"locus_id": "S2:c2:1-9", "island_size": "2", "n_strains": "2",
         "member_families": "famA", "pfam_domains": "-", "example_strain": "S2"},
    ]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["n_islands_total"] == 2
    assert payload["n_islands_excluded"] == 1
    assert len(payload["islands"]) == 1


def test_each_family_carries_its_pfam_class_for_the_glyph_strip():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "1", "n_strains": "2",
             "member_families": "famA", "pfam_domains": "NACHT",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["islands"][0]["dominant_class"] == "nlr"
    assert "nlr" in payload["classes"]


def test_zero_islands_produces_a_valid_empty_payload():
    # The single-species-ingroup case. A valid page saying so, not a crash.
    payload = build_payload([], make_matrix({}), {}, project="demo")
    assert payload["islands"] == []
    assert payload["n_islands_total"] == 0
    assert payload["project"] == "demo"


def test_families_absent_from_the_matrix_are_scored_absent_not_dropped():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA,ghost", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["islands"][0]["families"] == ["famA", "ghost"]
    for hap in payload["islands"][0]["haplotypes"]:
        assert len(hap["pattern"]) == 2
