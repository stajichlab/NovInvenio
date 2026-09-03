"""Unit tests for bin/busco_family_recovery.py (ADR-0002 Q8 BUSCO metric, issue #6).

A single-copy BUSCO should recover as exactly one gene family across the species that
carry it. Covers recovered / split / partial verdicts, the >= min-species eligibility
gate, over-clustering detection, ignoring non-Complete BUSCO rows, and both input paths
(raw full_table.tsv and a consolidated --busco-map).
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'busco_family_recovery.py'

_spec = importlib.util.spec_from_file_location('busco_family_recovery', BIN)
bfr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bfr)


# Family A (rep pA1): pA1, pA2, pA1b, pA2b   — two BUSCOs land here (over-merge).
# Family B1 (rep pB1) and B2 (rep pB2): the two copies of BUSCO B2 split apart.
CLUSTER = (
    "pA1\tpA1\npA1\tpA2\npA1\tpA1b\npA1\tpA2b\n"
    "pB1\tpB1\npB1\tpBfill\n"
    "pB2\tpB2\npB2\tpB2fill\n"
)
FAMILIES = (
    "family_index\trepresentative_id\tn_members\n"
    "fam_000001\tpA1\t4\n"
    "fam_000002\tpB1\t2\n"
    "fam_000003\tpB2\t2\n"
)

# BUSCO full_table.tsv (only cols 0,1,2 read). Non-Complete rows must be ignored.
NCRA_TABLE = (
    "# BUSCO full table\n"
    "B1\tComplete\tpA1\t1\t100\t+\t500\t250\n"
    "B2\tComplete\tpB1\t1\t100\t+\t500\t250\n"
    "B3\tComplete\tpOnly\t1\t100\t+\t500\t250\n"   # only in Ncra → not eligible
    "B4\tComplete\tpA1b\t1\t100\t+\t500\t250\n"
    "B5\tDuplicated\tpDup\t1\t100\t+\t500\t250\n"  # ignored (not single-copy)
)
AFUM_TABLE = (
    "# BUSCO full table\n"
    "B1\tComplete\tpA2:10-500\t1\t100\t+\t500\t250\n"  # :coords suffix must be stripped
    "B2\tComplete\tpB2\t1\t100\t+\t500\t250\n"
    "B4\tComplete\tpA2b\t1\t100\t+\t500\t250\n"
    "B6\tMissing\t\t\t\t\t\t\n"                        # ignored
)


def _setup(tmp_path):
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    (tmp_path / 'families.tsv').write_text(FAMILIES)
    (tmp_path / 'Ncra.full_table.tsv').write_text(NCRA_TABLE)
    (tmp_path / 'Afum.full_table.tsv').write_text(AFUM_TABLE)


def test_parse_busco_full_table_keeps_only_complete_and_strips_coords(tmp_path):
    (tmp_path / 't.tsv').write_text(AFUM_TABLE)
    rows = list(bfr.parse_busco_full_table(tmp_path / 't.tsv', 'Afum'))
    assert ('B1', 'Afum', 'pA2', 100) in rows      # :10-500 stripped, length=100 kept
    assert all(r[0] != 'B6' for r in rows)         # Missing dropped


def test_score_recovery_verdicts_and_overmerge(tmp_path):
    _setup(tmp_path)
    busco_map, _lengths = bfr.load_busco_map(
        tables=[f'Ncra={tmp_path/"Ncra.full_table.tsv"}',
                f'Afum={tmp_path/"Afum.full_table.tsv"}'])
    reps = bfr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bfr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    per_busco, summary = bfr.score_recovery(busco_map, m2r, min_species=2)

    verdicts = {r['busco_id']: r['verdict'] for r in per_busco}
    assert verdicts == {'B1': 'recovered', 'B2': 'split', 'B4': 'recovered'}  # B3 not eligible
    assert summary['eligible_buscos'] == 3
    assert summary['recovered_single_family'] == 2
    assert summary['split_across_families'] == 1
    assert summary['recovery_rate'] == 2 / 3
    # Family pA1 carries both B1 and B4 → one over-merged family.
    assert summary['overmerged_families'] == 1


def test_partial_when_member_dropped_as_singleton(tmp_path):
    _setup(tmp_path)
    # B7: one copy in a profiled family (pA1), one copy in an unprofiled singleton (pLone).
    busco_map = {'B7': {'Ncra': 'pA1', 'Afum': 'pLone'}}
    reps = bfr.load_profiled_reps(tmp_path / 'families.tsv')
    m2r = bfr.load_member_to_rep(tmp_path / 'cluster.tsv', reps)
    _, summary = bfr.score_recovery(busco_map, m2r, min_species=2)
    assert summary['partial_unmapped'] == 1
    assert summary['recovered_single_family'] == 0


def test_end_to_end_and_busco_map_input(tmp_path):
    _setup(tmp_path)
    out = tmp_path / 'recovery.tsv'
    subprocess.run([
        sys.executable, str(BIN),
        '--tables', f'Ncra={tmp_path/"Ncra.full_table.tsv"}',
                    f'Afum={tmp_path/"Afum.full_table.tsv"}',
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--families', str(tmp_path / 'families.tsv'),
        '--output', str(out),
    ], check=True, capture_output=True, text=True)
    summary = dict(
        line.split('\t')
        for line in (tmp_path / 'recovery.summary.tsv').read_text().splitlines()[1:]
    )
    assert summary['recovery_rate'] == str(2 / 3)

    # Same result from a consolidated --busco-map.
    (tmp_path / 'map.tsv').write_text(
        "B1\tNcra\tpA1\tComplete\nB1\tAfum\tpA2\tComplete\n"
        "B2\tNcra\tpB1\tComplete\nB2\tAfum\tpB2\tComplete\n"
        "B4\tNcra\tpA1b\tComplete\nB4\tAfum\tpA2b\tComplete\n"
    )
    bmap, _lengths = bfr.load_busco_map(map_path=tmp_path / 'map.tsv')
    assert set(bmap['B1']) == {'Ncra', 'Afum'}
    assert bmap['B1']['Afum'] == 'pA2'
