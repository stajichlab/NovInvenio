"""Unit tests for bin/busco_presence_recovery.py (2026-09-03 follow-up to ADR-0002 Q8).

Unlike bin/busco_family_recovery.py (clustering quality), this measures hmmsearch
presence-calling quality: does a family correctly show "present" for an OUTGROUP/
DISCOVERY_OUT species that BUSCO independently says carries a Complete copy? A seed-group
species' own presence is trivial by clustering construction and must be excluded from
scoring -- that's the main behavioural contract under test here.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'busco_presence_recovery.py'

_spec = importlib.util.spec_from_file_location('busco_presence_recovery', BIN)
bpr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bpr)

CLUSTER = "pA1\tpA1\npA1\tpA2\n"
FAMILIES = "family_index\trepresentative_id\tn_members\nfam_000001\tpA1\t2\n"
MATRIX_HEADER = "protein_id\tsource_proteome\tIn1\tIn2\tOut1\n"


def _setup(tmp_path, out1_present):
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    (tmp_path / 'families.tsv').write_text(FAMILIES)
    (tmp_path / 'In1.full_table.tsv').write_text(
        "# BUSCO full table\nB1\tComplete\tpA1\t1\t300\t+\t500\t250\n")
    (tmp_path / 'In2.full_table.tsv').write_text(
        "# BUSCO full table\nB1\tComplete\tpA2\t1\t300\t+\t500\t250\n")
    (tmp_path / 'Out1.full_table.tsv').write_text(
        "# BUSCO full table\nB1\tComplete\tpOut1_x\t1\t300\t+\t500\t250\n")
    (tmp_path / 'matrix.tsv').write_text(
        MATRIX_HEADER + f"pA1\tIn1\t1\t1\t{int(out1_present)}\n"
                       + f"pA2\tIn2\t1\t1\t{int(out1_present)}\n"
    )


def _tables(tmp_path, *shorts):
    return [f'{s}={tmp_path / f"{s}.full_table.tsv"}' for s in shorts]


def test_seed_group_only_species_yield_zero_scorable_pairs(tmp_path):
    # Both In1/In2 are family members (seed group) -- their own presence is trivial by
    # construction, so with no outgroup table there is nothing left to test.
    _setup(tmp_path, out1_present=True)
    busco_map, lengths = bpr.load_busco_map(tables=_tables(tmp_path, 'In1', 'In2'))
    reps = bpr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bpr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    matrix_rows = bpr.read_matrix_rows(tmp_path / 'matrix.tsv')
    per_pair, summary, _ = bpr.score_presence_recovery(busco_map, lengths, m2r, matrix_rows,
                                                        min_species=2)
    assert summary['species_pairs_scored'] == 0
    assert summary['recovered_buscos_scored'] == 1


def test_outgroup_species_present_scores_a_hit(tmp_path):
    _setup(tmp_path, out1_present=True)
    busco_map, lengths = bpr.load_busco_map(tables=_tables(tmp_path, 'In1', 'In2', 'Out1'))
    reps = bpr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bpr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    matrix_rows = bpr.read_matrix_rows(tmp_path / 'matrix.tsv')
    per_pair, summary, _ = bpr.score_presence_recovery(busco_map, lengths, m2r, matrix_rows,
                                                        min_species=2)
    assert summary['species_pairs_scored'] == 1
    assert summary['species_pairs_present'] == 1
    assert summary['presence_recovery_rate'] == 1.0
    assert per_pair == [{
        'busco_id': 'B1', 'species': 'Out1', 'family_rep': 'pA1',
        'length': 300, 'length_bucket': '300-600aa', 'verdict': 'present',
    }]


def test_outgroup_species_absent_scores_a_miss(tmp_path):
    # Regression scenario for the coverage-floor bug: BUSCO independently confirms Out1
    # carries a Complete copy, but the presence matrix (e.g. hmm_presence_cov too strict)
    # says absent -- this is exactly the false-negative this metric exists to surface.
    _setup(tmp_path, out1_present=False)
    busco_map, lengths = bpr.load_busco_map(tables=_tables(tmp_path, 'In1', 'In2', 'Out1'))
    reps = bpr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bpr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    matrix_rows = bpr.read_matrix_rows(tmp_path / 'matrix.tsv')
    _, summary, _ = bpr.score_presence_recovery(busco_map, lengths, m2r, matrix_rows,
                                                 min_species=2)
    assert summary['species_pairs_absent'] == 1
    assert summary['presence_recovery_rate'] == 0.0


def test_busco_split_across_families_is_excluded(tmp_path):
    # Mirrors busco_family_recovery.py's "split" verdict -- ambiguous family row, must not
    # be scored here (that failure mode belongs to clustering, not presence-calling).
    (tmp_path / 'cluster.tsv').write_text("pA1\tpA1\npB1\tpB1\n")
    (tmp_path / 'families.tsv').write_text(
        "family_index\trepresentative_id\tn_members\n"
        "fam_000001\tpA1\t1\nfam_000002\tpB1\t1\n"
    )
    busco_map = {'B1': {'In1': 'pA1', 'In2': 'pB1'}}  # same BUSCO, two different families
    reps = bpr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bpr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    result = bpr.recovered_buscos(busco_map, m2r, min_species=2)
    assert result == {}


def test_bucket_for_boundaries():
    assert bpr.bucket_for(149) == '<150aa'
    assert bpr.bucket_for(150) == '150-300aa'
    assert bpr.bucket_for(300) == '300-600aa'
    assert bpr.bucket_for(600) == '>600aa'
    assert bpr.bucket_for(None) == 'unknown'


def test_length_representativeness_compares_eligible_vs_reference(tmp_path):
    dom = tmp_path / 'ref.domtblout'
    lines = []
    for i, qlen in enumerate([100, 200, 700, 800]):
        f = ['-'] * 23
        f[0] = f'target{i}'
        f[3] = f'fam{i}'
        f[5] = str(qlen)
        lines.append(' '.join(f))
    dom.write_text("# hmmsearch domtblout\n" + "\n".join(lines) + "\n")

    rep = bpr.length_representativeness([100, 200], str(dom))
    assert rep['eligible_busco']['n'] == 2
    assert rep['reference']['n'] == 4


def test_end_to_end_cli(tmp_path):
    _setup(tmp_path, out1_present=True)
    out = tmp_path / 'recovery.tsv'
    subprocess.run([
        sys.executable, str(BIN),
        '--tables', *_tables(tmp_path, 'In1', 'In2', 'Out1'),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--families', str(tmp_path / 'families.tsv'),
        '--matrix', str(tmp_path / 'matrix.tsv'),
        '--output', str(out),
    ], check=True, capture_output=True, text=True)
    summary = dict(
        line.split('\t')
        for line in (tmp_path / 'recovery.summary.tsv').read_text().splitlines()[1:]
    )
    assert summary['presence_recovery_rate'] == '1.0'
