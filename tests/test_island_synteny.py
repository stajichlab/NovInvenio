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


# ---- C1: multi-copy families must order by the LOCUS CONTIG's copy, not the
# highest-ranked (often off-contig) paralog. family_positions.tsv has one row
# per (strain, family, COPY), so a family with a paralog elsewhere in the
# genome has more than one (contig, rank) entry.

def test_multicopy_family_orders_by_the_locus_contig_not_a_stray_paralog():
    # famA has a copy at rank 1 on c1 (the island's own locus contig) AND a
    # stray paralog at rank 900 on c9. Naive last-row-wins / max-rank
    # behaviour would sort famA after famB; the correct column order follows
    # the c1 copy.
    positions = {
        ("S1", "famA"): [("c1", 1), ("c9", 900)],
        ("S1", "famB"): [("c1", 2)],
    }
    assert order_families_by_locus(["famB", "famA"], positions, "S1", contig="c1") == \
        ["famA", "famB"]


def test_multicopy_family_falls_back_to_global_minimum_off_contig():
    # No copy of famA on the requested contig c2 -- fall back to its global
    # minimum rank (1, from c1) rather than dropping the column. famB's
    # on-contig rank (5) is higher, so famA still sorts first.
    positions = {
        ("S1", "famA"): [("c1", 1), ("c9", 900)],
        ("S1", "famB"): [("c2", 5)],
    }
    assert order_families_by_locus(["famB", "famA"], positions, "S1", contig="c2") == \
        ["famA", "famB"]


def test_single_copy_list_form_behaves_like_the_legacy_int_form():
    positions_list = {("S1", "famA"): [("c1", 1)], ("S1", "famB"): [("c1", 2)]}
    positions_int = {("S1", "famA"): 1, ("S1", "famB"): 2}
    assert (order_families_by_locus(["famB", "famA"], positions_list, "S1", contig="c1")
            == order_families_by_locus(["famB", "famA"], positions_int, "S1")
            == ["famA", "famB"])


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


def test_payload_orders_columns_by_locus_contig_for_a_multicopy_family():
    # C1 regression: build_payload must pass the island's own locus_contig
    # through to order_families_by_locus so a stray off-contig paralog can't
    # drag a column to the wrong position.
    #
    # The fixture is chosen so it DISCRIMINATES. An earlier version had famA's
    # paralog on "c9" -- alphabetically after the locus contig -- which meant a
    # naive lexicographic sort over the raw [(contig, rank), ...] lists produced
    # the correct answer by coincidence, so the test passed against the very bug
    # it was written to catch. Here famA's paralog sits on "a0" (alphabetically
    # BEFORE the locus contig) while its real on-locus copy is at rank 5, after
    # famB's rank 2. Correct answer: ["famB", "famA"]. A naive sort that ignores
    # locus_contig compares ("a0", 1) against ("c1", 2) and yields
    # ["famA", "famB"] -- the wrong order, so the bug is now actually detected.
    rows = [{"locus_id": "S1:c1:1-400", "locus_contig": "c1", "island_size": "2",
             "n_strains": "2", "member_families": "famB,famA",
             "pfam_domains": "-", "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True},
        "famB": {"S1": True, "S2": False},
    })
    positions = {
        ("S1", "famA"): [("a0", 1), ("c1", 5)],
        ("S1", "famB"): [("c1", 2)],
    }
    payload = build_payload(rows, matrix, positions, project="demo")
    assert payload["islands"][0]["families"] == ["famB", "famA"]


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


# ---- issue #119: "by species" row sort -- payload carries a strain->species map

def test_payload_carries_species_when_provided():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True, "S3": False}})
    species_of = {"S1": "Coccidioides immitis", "S2": "Coccidioides posadasii",
                  "S3": "Coccidioides immitis"}
    payload = build_payload(rows, matrix, {}, project="demo", species_of=species_of)
    assert payload["species"] == species_of


def test_payload_omits_species_when_not_provided():
    # Graceful degradation (issue #119): callers that never pass a config
    # (bin/pangenome_island_synteny.py's --config is optional) must still
    # get a valid payload -- the species sort option simply isn't offered.
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["species"] == {}


def test_payload_species_is_missing_strain_does_not_crash():
    # A Short present in the presence matrix but absent from the config
    # (e.g. a strain added after the config CSV was last touched) must not
    # raise -- it's simply absent from the species map, not an error.
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo",
                            species_of={"S1": "Coccidioides immitis"})
    assert payload["species"] == {"S1": "Coccidioides immitis"}
    assert "S2" not in payload["species"]


def test_families_absent_from_the_matrix_are_scored_absent_not_dropped():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA,ghost", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo")
    assert payload["islands"][0]["families"] == ["famA", "ghost"]
    for hap in payload["islands"][0]["haplotypes"]:
        assert len(hap["pattern"]) == 2


def test_n_islands_excluded_counts_only_min_strains_failures_not_top_n_truncation():
    # Regression: n_islands_excluded used to conflate filter exclusions with
    # top_n truncation. 60 all-multi-strain islands, top_n=50 should report
    # n_islands_excluded=0 (none failed the filter) and n_islands_truncated=10
    # (qualified but beyond top_n).
    rows = [
        {"locus_id": f"S{i}:c{i}:1-400", "island_size": str(100 - i),
         "n_strains": "2", "member_families": "famA", "pfam_domains": "-",
         "example_strain": f"S{i}"}
        for i in range(60)
    ]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo", top_n=50)
    assert payload["n_islands_total"] == 60
    assert payload["n_islands_excluded"] == 0  # all multi-strain, none failed filter
    assert payload["n_islands_truncated"] == 10  # 60 qualified, 50 selected


def test_mixed_exclusion_and_truncation_must_sum_correctly():
    # Both single-strain (excluded) and qualifying (some truncated) islands.
    # n_islands_total = n_islands_excluded + n_islands_truncated + len(islands)
    rows = [
        {"locus_id": f"S{i}:c{i}:1-400", "island_size": str(100 - i),
         "n_strains": str(1 if i < 10 else 2),  # first 10 are single-strain
         "member_families": "famA", "pfam_domains": "-",
         "example_strain": f"S{i}"}
        for i in range(60)  # 10 single-strain + 50 multi-strain
    ]
    matrix = make_matrix({"famA": {"S1": True, "S2": True}})
    payload = build_payload(rows, matrix, {}, project="demo", top_n=30)
    # 60 total, 10 single-strain (excluded), 50 qualifying multi-strain,
    # 30 selected, so 20 truncated
    assert payload["n_islands_total"] == 60
    assert payload["n_islands_excluded"] == 10
    assert payload["n_islands_truncated"] == 20
    assert len(payload["islands"]) == 30
    assert (payload["n_islands_total"] ==
            payload["n_islands_excluded"] + payload["n_islands_truncated"] +
            len(payload["islands"]))


# ---- per-family Pfam class (ruling R8: the glyph strip must be per-column,
# not island-wide) --------------------------------------------------------

def test_per_family_classes_come_from_the_supplied_family_domains_mapping():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famB,famA", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True},
        "famB": {"S1": True, "S2": True},
    })
    positions = {("S1", "famA"): 1, ("S1", "famB"): 2}
    family_domains = {"famA": {"NACHT"}, "famB": {"MFS_1"}}
    payload = build_payload(rows, matrix, positions, project="demo",
                             family_domains=family_domains)
    island = payload["islands"][0]
    assert island["families"] == ["famA", "famB"]
    assert island["family_classes"] == ["nlr", "transporter"]
    assert island["family_domains"] == ["NACHT", "MFS_1"]


def test_a_family_absent_from_family_domains_is_unannotated():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA,famB", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True},
        "famB": {"S1": True, "S2": True},
    })
    positions = {("S1", "famA"): 1, ("S1", "famB"): 2}
    family_domains = {"famA": {"NACHT"}}  # famB not covered by the scan
    payload = build_payload(rows, matrix, positions, project="demo",
                             family_domains=family_domains)
    island = payload["islands"][0]
    assert island["family_classes"] == ["nlr", "unannotated"]
    assert island["family_domains"] == ["NACHT", ""]


def test_omitting_family_domains_degrades_gracefully():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "2", "n_strains": "2",
             "member_families": "famA,famB", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True},
        "famB": {"S1": True, "S2": True},
    })
    positions = {("S1", "famA"): 1, ("S1", "famB"): 2}
    payload = build_payload(rows, matrix, positions, project="demo")
    island = payload["islands"][0]
    assert island["family_classes"] == ["unannotated", "unannotated"]
    assert island["family_domains"] == ["", ""]


def test_family_class_and_domain_arrays_are_the_same_length_as_families():
    rows = [{"locus_id": "S1:c1:1-400", "island_size": "3", "n_strains": "2",
             "member_families": "famA,famB,famC", "pfam_domains": "-",
             "example_strain": "S1"}]
    matrix = make_matrix({
        "famA": {"S1": True, "S2": True},
        "famB": {"S1": True, "S2": True},
        "famC": {"S1": True, "S2": True},
    })
    family_domains = {"famA": {"NACHT"}}
    payload = build_payload(rows, matrix, {}, project="demo",
                             family_domains=family_domains)
    island = payload["islands"][0]
    assert len(island["family_classes"]) == len(island["families"])
    assert len(island["family_domains"]) == len(island["families"])
