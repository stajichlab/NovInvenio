"""Unit tests for bin/candidate_families.py (ADR-0002 Q7 family-as-cluster, issue #11).

Verifies the profile pathway's gene families (restricted to candidate-containing families)
reproduce the CLUSTER workflow's three-artifact contract — cluster_tsv, representatives,
candidates_fa — without a second mmseqs run.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'candidate_families.py'

_spec = importlib.util.spec_from_file_location('candidate_families', BIN)
sys.path.insert(0, str(REPO / 'lib'))
cf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cf)

# Families: A (rep pA1: pA1,pA2), B (rep pB1: pB1,pB2), C (rep pC1: pC1,pC2).
FAMILY_CLUSTER = "pA1\tpA1\npA1\tpA2\npB1\tpB1\npB1\tpB2\npC1\tpC1\npC1\tpC2\n"
FAMILY_REPS = ">pA1\nMAAA\n>pB1\nMBBB\n>pC1\nMCCC\n"
INGROUP_FA = ">pA1\nMAAA\n>pA2\nMAAX\n>pB1\nMBBB\n>pB2\nMBBX\n>pC1\nMCCC\n>pC2\nMCCX\n"
# Candidates drawn from families A and B only (C is not a candidate).
CANDIDATES = "Ncra::pA1\nAfum::pA2\nNcra::pB1\n"


def _setup(tmp_path):
    (tmp_path / 'family_cluster.tsv').write_text(FAMILY_CLUSTER)
    (tmp_path / 'family_reps.fasta').write_text(FAMILY_REPS)
    (tmp_path / 'ingroup.fasta').write_text(INGROUP_FA)
    (tmp_path / 'candidates.txt').write_text(CANDIDATES)


def test_load_candidate_ids_strips_source_prefix(tmp_path):
    (tmp_path / 'c.txt').write_text(CANDIDATES)
    assert cf.load_candidate_ids(tmp_path / 'c.txt') == {'pA1', 'pA2', 'pB1'}


def _write(tmp_path, text, name='fc.tsv'):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_candidate_family_reps_only_families_with_a_candidate(tmp_path):
    m2r, _ = cf.load_family_membership(_write(tmp_path, FAMILY_CLUSTER))
    reps = cf.candidate_family_reps({'pA1', 'pA2', 'pB1'}, m2r)
    assert reps == {'pA1', 'pB1'}   # family C (pC1) excluded — no candidate


def _run(tmp_path):
    subprocess.run([
        sys.executable, str(BIN),
        '--candidates', str(tmp_path / 'candidates.txt'),
        '--family-cluster-tsv', str(tmp_path / 'family_cluster.tsv'),
        '--family-reps', str(tmp_path / 'family_reps.fasta'),
        '--ingroup-fasta', str(tmp_path / 'ingroup.fasta'),
        '--out-cluster-tsv', str(tmp_path / 'out_cluster.tsv'),
        '--out-representatives', str(tmp_path / 'out_reps.fasta'),
        '--out-candidates-fa', str(tmp_path / 'out_candidates.fa'),
    ], check=True, capture_output=True, text=True)


def test_end_to_end_emits_candidate_family_artifacts(tmp_path):
    _setup(tmp_path)
    _run(tmp_path)

    # cluster_tsv: only candidate families (A, B) with full membership; C excluded.
    cluster_lines = set(
        (tmp_path / 'out_cluster.tsv').read_text().splitlines())
    assert cluster_lines == {'pA1\tpA1', 'pA1\tpA2', 'pB1\tpB1', 'pB1\tpB2'}

    # representatives: reps of the candidate families only.
    from fasta import extract_ids
    assert extract_ids(tmp_path / 'out_reps.fasta') == {'pA1', 'pB1'}

    # candidates_fa: the candidate protein sequences (pA1, pA2, pB1), from the ingroup.
    assert extract_ids(tmp_path / 'out_candidates.fa') == {'pA1', 'pA2', 'pB1'}


def test_candidate_outside_any_family_still_gets_a_sequence(tmp_path):
    _setup(tmp_path)
    # pLone is a candidate that is in no profiled family (singleton edge case).
    (tmp_path / 'ingroup.fasta').write_text(INGROUP_FA + '>pLone\nMLONE\n')
    (tmp_path / 'candidates.txt').write_text(CANDIDATES + 'Ncra::pLone\n')
    _run(tmp_path)
    from fasta import extract_ids
    # pLone has no family → absent from cluster_tsv/reps, but its sequence is still emitted.
    assert 'pLone' in extract_ids(tmp_path / 'out_candidates.fa')
    assert cf.load_family_membership(tmp_path / 'family_cluster.tsv')[0].get('pLone') is None
