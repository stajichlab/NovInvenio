import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import (
    PRESENT, GENOME_ONLY, ABSENT,
    read_cluster_tsv, build_families, PresenceMatrix, copy_number_path,
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


def test_from_tsv_rejects_unknown_state(tmp_path):
    p = tmp_path / "corrupt.tsv"
    p.write_text("family\ts1\ts2\nfamA\tpresent\tPRESENT\n")
    with pytest.raises(ValueError, match="invalid presence state"):
        PresenceMatrix.from_tsv(p)


def test_set_call_rejects_invalid_state():
    pm = PresenceMatrix(families=["famA"], strains=["s1"])
    with pytest.raises(ValueError):
        pm.set_call("famA", "s1", "not_a_real_state")
