"""Unit tests for bin/profile_to_matrix.py.

Verifies the family-profile pathway reproduces build_presence_matrix.py's output
contract: the same `protein_id, source_proteome, <sorted shorts>` matrix columns and the
same `source::protein` candidates format, driven by hmmsearch domtblout presence +
cluster membership. See ADR-0002.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parent.parent / 'bin' / 'profile_to_matrix.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Species one,,In1.faa,In1.fna,In1,Pezizomycotina
IN,Species two,,In2.faa,In2.fna,In2,Pezizomycotina
OUT,Species three,,Out1.faa,Out1.fna,Out1,Saccharomycotina
"""


def _domtblout_line(target, query, qlen=100, full_e=1e-10, hmm_from=1, hmm_to=60):
    """One hmmsearch --domtblout data row (23 whitespace fields); only cols
    0,3,5,6,15,16 are read by the parser."""
    f = ['-'] * 23
    f[0] = target          # target name
    f[3] = query           # query (HMM) name = family rep id
    f[5] = str(qlen)       # qlen (HMM length)
    f[6] = f"{full_e:g}"   # full-sequence E-value
    f[15] = str(hmm_from)  # hmm coord from
    f[16] = str(hmm_to)    # hmm coord to
    return ' '.join(f)


def _write_domtblout(path, hits):
    """hits: list of (target, query) tuples → a domtblout file with a header comment."""
    with open(path, 'w') as fh:
        fh.write("# hmmsearch domtblout\n")
        for target, query in hits:
            fh.write(_domtblout_line(target, query) + "\n")


def _setup(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)

    # Two clade families (A: pA1/pA2, B: pB1/pB2) + a dropped singleton (pC1).
    (tmp_path / 'cluster.tsv').write_text(
        "pA1\tpA1\npA1\tpA2\npB1\tpB1\npB1\tpB2\npC1\tpC1\n"
    )
    # families.tsv lists only the profiled (>=2 member) families — singleton excluded.
    (tmp_path / 'families.tsv').write_text(
        "family_index\trepresentative_id\tn_members\n"
        "fam_000001\tpA1\t2\n"
        "fam_000002\tpB1\t2\n"
    )
    (tmp_path / 'protein_map.tsv').write_text(
        "pA1\tIn1\npA2\tIn2\npB1\tIn1\npB2\tIn2\npC1\tIn1\n"
    )

    # Family A present in both ingroup, absent from the outgroup (→ novelty candidate).
    # Family B present in both ingroup AND the outgroup (→ not a candidate).
    _write_domtblout(tmp_path / 'In1.family.domtblout', [('pA1', 'pA1'), ('pB1', 'pB1')])
    _write_domtblout(tmp_path / 'In2.family.domtblout', [('pA2', 'pA1'), ('pB2', 'pB1')])
    _write_domtblout(tmp_path / 'Out1.family.domtblout', [('oX', 'pB1')])


def _run(tmp_path, **extra):
    args = [
        sys.executable, str(BIN),
        '--domtblout',
        str(tmp_path / 'In1.family.domtblout'),
        str(tmp_path / 'In2.family.domtblout'),
        str(tmp_path / 'Out1.family.domtblout'),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--families', str(tmp_path / 'families.tsv'),
        '--protein-map', str(tmp_path / 'protein_map.tsv'),
        '--config', str(tmp_path / 'config.csv'),
        '--output-matrix', str(tmp_path / 'matrix.tsv'),
        '--output-candidates', str(tmp_path / 'candidates.txt'),
    ]
    for k, v in extra.items():
        args += [f'--{k}', str(v)]
    subprocess.run(args, check=True)
    matrix = pd.read_csv(tmp_path / 'matrix.tsv', sep='\t')
    cands = [c for c in (tmp_path / 'candidates.txt').read_text().splitlines() if c]
    return matrix, cands


def test_matrix_columns_match_contract(tmp_path):
    _setup(tmp_path)
    matrix, _ = _run(tmp_path)
    assert list(matrix.columns) == ['protein_id', 'source_proteome', 'In1', 'In2', 'Out1']


def test_singleton_family_dropped(tmp_path):
    _setup(tmp_path)
    matrix, _ = _run(tmp_path)
    # pC1 is a singleton not in families.tsv → no row.
    assert 'pC1' not in set(matrix['protein_id'])
    assert set(matrix['protein_id']) == {'pA1', 'pA2', 'pB1', 'pB2'}


def test_presence_values_and_source_present(tmp_path):
    _setup(tmp_path)
    matrix, _ = _run(tmp_path)
    m = matrix.set_index('protein_id')
    # Family A: present in both ingroup, absent from outgroup.
    assert m.loc['pA1', ['In1', 'In2', 'Out1']].tolist() == [1, 1, 0]
    assert m.loc['pA2', ['In1', 'In2', 'Out1']].tolist() == [1, 1, 0]
    # Family B: present everywhere including the outgroup.
    assert m.loc['pB1', ['In1', 'In2', 'Out1']].tolist() == [1, 1, 1]


def test_candidates_are_novelty_families_only(tmp_path):
    _setup(tmp_path)
    _, cands = _run(tmp_path)
    # Family A members are candidates (conserved ingroup, absent outgroup); B members are not.
    assert set(cands) == {'In1::pA1', 'In2::pA2'}


def test_coverage_filter_excludes_short_hits(tmp_path):
    _setup(tmp_path)
    # Rewrite In1/In2 family-A hits with only 30% HMM coverage → family A no longer
    # present in the ingroup, so its members drop out of the candidate set.
    for short, member in (('In1', 'pA1'), ('In2', 'pA2')):
        with open(tmp_path / f'{short}.family.domtblout', 'w') as fh:
            fh.write("# hmmsearch domtblout\n")
            fh.write(_domtblout_line(member, 'pA1', hmm_from=1, hmm_to=30) + "\n")
            fh.write(_domtblout_line(f'{member}b', 'pB1') + "\n")
    _, cands = _run(tmp_path)
    assert cands == []


def test_min_covered_residues_rescues_a_long_partial_hit(tmp_path):
    # Same 30%-coverage scenario as above, but on a long (2000 aa) HMM where the aligned
    # span (600 residues) is substantial in absolute terms despite the low fraction --
    # --min-covered-residues rescues it as real ingroup presence.
    _setup(tmp_path)
    for short, member in (('In1', 'pA1'), ('In2', 'pA2')):
        with open(tmp_path / f'{short}.family.domtblout', 'w') as fh:
            fh.write("# hmmsearch domtblout\n")
            fh.write(_domtblout_line(member, 'pA1', qlen=2000, hmm_from=1, hmm_to=600) + "\n")
            fh.write(_domtblout_line(f'{member}b', 'pB1') + "\n")

    _, cands = _run(tmp_path)
    assert cands == []  # fraction-only (30%): family A absent from ingroup, no candidates

    _, cands = _run(tmp_path, **{'min-covered-residues': 100})
    assert set(cands) == {'In1::pA1', 'In2::pA2'}
