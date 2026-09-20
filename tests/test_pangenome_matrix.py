import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import (
    PRESENT, GENOME_ONLY, ABSENT,
    read_cluster_tsv, build_families, PresenceMatrix, copy_number_path,
    iter_tblout_family_hits,
)


def test_read_cluster_tsv_parses_rep_member_pairs(tmp_path):
    p = tmp_path / "clusters.tsv"
    p.write_text("repA\trepA\nrepA\tmemberA2\nrepB\trepB\n")
    result = read_cluster_tsv(p)
    assert result == {"repA": "repA", "memberA2": "repA", "repB": "repB"}


def test_build_families_groups_by_rep_including_singletons():
    member_to_rep = {"repA": "repA", "memberA2": "repA", "repB": "repB"}
    families = build_families(member_to_rep)
    assert families == {"repA": ["memberA2", "repA"], "repB": ["repB"]}


def test_presence_matrix_set_and_query():
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2", "s3"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famA", "s2", GENOME_ONLY)
    # famA absent in s3 (never set) -- default ABSENT
    assert pm.call("famA", "s1") == PRESENT
    assert pm.call("famA", "s3") == ABSENT
    assert pm.is_present("famA", "s1") is True
    assert pm.is_present("famA", "s2") is True   # genome_only still counts as present
    assert pm.is_present("famA", "s3") is False
    assert pm.presence_vector("famA") == [True, True, False]
    assert pm.strain_count("famA") == 2
    assert pm.frequency("famA") == 2 / 3


def test_iter_tblout_family_hits_default_id_sep(tmp_path):
    tblout = tmp_path / "marker.tblout"
    tblout.write_text(
        "# comment\n"
        "s1|proteinA - - - - - - - - - - - - - - - -\n"
        "s2|proteinB - - - - - - - - - - - - - - - -\n"
    )
    member_to_rep = {"s1|proteinA": "famA", "s2|proteinB": "famB"}
    hits = sorted(iter_tblout_family_hits(str(tblout), member_to_rep))
    assert hits == [("s1", "famA"), ("s2", "famB")]


def test_iter_tblout_family_hits_custom_id_sep(tmp_path):
    # Confirms a non-default id_sep is actually threaded through, not
    # hardcoded -- the bug this consolidation fixes (F6): a study configured
    # with a different --pangenome_id_sep previously got zero marker hits
    # from either load_hit_families or load_captain_families.
    tblout = tmp_path / "marker.tblout"
    tblout.write_text("s1__proteinA - - - - - - - - - - - - - - - -\n")
    member_to_rep = {"s1__proteinA": "famA"}
    # Default id_sep "|" finds nothing (no "|" in the target name at all).
    assert list(iter_tblout_family_hits(str(tblout), member_to_rep)) == []
    # Custom id_sep "__" correctly resolves the hit.
    hits = list(iter_tblout_family_hits(str(tblout), member_to_rep, id_sep="__"))
    assert hits == [("s1", "famA")]


def test_presence_matrix_tsv_roundtrip(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ["s1", "s2"]
    assert loaded.families == ["famA", "famB"]
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY
    assert loaded.call("famA", "s2") == ABSENT


def test_tsv_roundtrip_preserves_copy_number(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT, copies=3)
    pm.set_call("famA", "s2", PRESENT, copies=1)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)

    assert copy_number_path(out).exists()
    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.copy_number[("famA", "s1")] == 3
    assert loaded.copy_number[("famA", "s2")] == 1
    assert loaded.copy_number.get(("famB", "s2"), 0) == 0
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY


@pytest.mark.parametrize("suffix", [".gz", ".zst"])
def test_presence_matrix_tsv_roundtrip_compressed(tmp_path, suffix):
    """to_tsv actually compresses when given a .gz/.zst path (not just
    accepts one on read) -- the write-side half of this repo's general
    storage-compression convention (see ~/.claude/CLAUDE.md), matching
    from_tsv's existing transparent read support."""
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / f"matrix.tsv{suffix}"
    pm.to_tsv(out)

    assert out.exists()
    raw = out.read_bytes()
    assert not raw.startswith(b"family\t")

    loaded = PresenceMatrix.from_tsv(out)
    assert loaded.strains == ["s1", "s2"]
    assert loaded.families == ["famA", "famB"]
    assert loaded.call("famA", "s1") == PRESENT
    assert loaded.call("famB", "s2") == GENOME_ONLY
    assert loaded.call("famA", "s2") == ABSENT


def _write_wide_matrix(tmp_path, n_families=1000, n_strains=8, wanted=("famA", "famB")):
    """A matrix with many families -- only `wanted` should be loaded when
    `families=` is passed, mirroring the real island-synteny case where 298
    of 54,421 families are actually referenced by the islands drawn."""
    strains = [f"s{i}" for i in range(n_strains)]
    lines = ["family\t" + "\t".join(strains)]
    for i in range(n_families):
        fam = wanted[i] if i < len(wanted) else f"fam{i:05d}"
        # every strain present, alternating for variety
        row = ["present" if (i + j) % 2 == 0 else "absent" for j in range(n_strains)]
        lines.append(fam + "\t" + "\t".join(row))
    p = tmp_path / "wide_matrix.tsv"
    p.write_text("\n".join(lines) + "\n")
    return p, strains


def test_from_tsv_default_behaviour_is_unchanged_without_families_arg(tmp_path):
    """Pin: from_tsv(path) with no `families` arg must be byte-identical to
    today's behaviour -- six other consumers depend on the unfiltered load."""
    pm = PresenceMatrix(families=["famA", "famB"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT)
    pm.set_call("famB", "s2", GENOME_ONLY)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)

    loaded_default = PresenceMatrix.from_tsv(out)
    loaded_explicit_none = PresenceMatrix.from_tsv(out, families=None)
    assert loaded_default.families == loaded_explicit_none.families == ["famA", "famB"]
    assert loaded_default.strains == loaded_explicit_none.strains == ["s1", "s2"]
    assert loaded_default.calls == loaded_explicit_none.calls
    assert loaded_default.copy_number == loaded_explicit_none.copy_number


def test_from_tsv_with_families_loads_only_requested_families(tmp_path):
    p, strains = _write_wide_matrix(tmp_path, n_families=1000)
    loaded = PresenceMatrix.from_tsv(p, families={"famA", "famB"})
    assert sorted(loaded.families) == ["famA", "famB"]
    # Memory characteristic: an order of magnitude fewer calls than a full load.
    assert len(loaded.calls) == 2 * len(strains)


def test_from_tsv_with_families_strains_is_the_full_list(tmp_path):
    """matrix.strains comes from the header row and must stay the FULL
    strain list even when only a handful of families are loaded -- the
    grid's rows are strains, so truncating strains would be wrong."""
    p, strains = _write_wide_matrix(tmp_path, n_families=500, n_strains=12)
    loaded = PresenceMatrix.from_tsv(p, families={"famA"})
    assert loaded.strains == strains


def test_from_tsv_with_families_requested_family_absent_does_not_crash(tmp_path):
    p, strains = _write_wide_matrix(tmp_path, n_families=50)
    loaded = PresenceMatrix.from_tsv(p, families={"famA", "does_not_exist"})
    assert loaded.families == ["famA"]
    # An absent family must still score absent (default) rather than raising.
    for s in strains:
        assert loaded.call("does_not_exist", s) == ABSENT


def test_from_tsv_with_families_matches_full_load_for_those_families(tmp_path):
    """The filtered load is an optimisation, not a behaviour change: for the
    families it does load, results must match a full unfiltered load exactly."""
    p, strains = _write_wide_matrix(tmp_path, n_families=300)
    full = PresenceMatrix.from_tsv(p)
    filtered = PresenceMatrix.from_tsv(p, families={"famA", "famB"})
    for fam in ("famA", "famB"):
        for s in strains:
            assert filtered.call(fam, s) == full.call(fam, s)


def test_from_tsv_with_families_filters_copy_number_sidecar(tmp_path):
    pm = PresenceMatrix(families=["famA", "famB", "famC"], strains=["s1", "s2"])
    pm.set_call("famA", "s1", PRESENT, copies=3)
    pm.set_call("famB", "s1", PRESENT, copies=2)
    pm.set_call("famC", "s1", PRESENT, copies=5)
    out = tmp_path / "matrix.tsv"
    pm.to_tsv(out)

    loaded = PresenceMatrix.from_tsv(out, families={"famA"})
    assert loaded.copy_number == {("famA", "s1"): 3}


def test_from_tsv_rejects_unknown_state(tmp_path):
    p = tmp_path / "corrupt.tsv"
    p.write_text("family\ts1\ts2\nfamA\tpresent\tPRESENT\n")
    with pytest.raises(ValueError, match="invalid presence state"):
        PresenceMatrix.from_tsv(p)


def test_set_call_rejects_invalid_state():
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    with pytest.raises(ValueError):
        pm.set_call("famA", "s1", "not_a_real_state")
