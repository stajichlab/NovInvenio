"""Unit tests for bin/novelty_presence_matrix.py."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parent.parent / 'bin' / 'novelty_presence_matrix.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
DISCOVERY_TARGET,Species one,,target1.pep.fa,target1.dna.fa,T1,X
DISCOVERY_TARGET,Species two,,target2.pep.fa,target2.dna.fa,T2,X
DISCOVERY_OUT,Outgroup one,,disc1.pep.fa,disc1.dna.fa,D1,Y
DISCOVERY_OUT,Outgroup two,,disc2.pep.fa,disc2.dna.fa,D2,Y
"""


def _domtblout_line(target, query, full_e=1e-10, qlen=100, hmm_from=1, hmm_to=60):
    f = ['-'] * 23
    f[0] = target
    f[3] = query
    f[5] = str(qlen)
    f[6] = f'{full_e:g}'
    f[15] = str(hmm_from)
    f[16] = str(hmm_to)
    return ' '.join(f)


def _write_domtblout(path, hits):
    with open(path, 'w') as fh:
        fh.write("# hmmsearch domtblout\n")
        for target, query, *rest in hits:
            evalue = rest[0] if rest else 1e-10
            fh.write(_domtblout_line(target, query, full_e=evalue) + "\n")


def _setup(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)

    # Cluster TSV: famA has 2 members (T1+pA1, T2+pA2), pS1 is a singleton from T1.
    (tmp_path / 'cluster.tsv').write_text(
        "pA1\tpA1\npA1\tpA2\npS1\tpS1\n"
    )
    # Protein map: protein_id → proteome_short
    (tmp_path / 'protein_map.tsv').write_text(
        "pA1\tT1\npA2\tT2\npS1\tT1\n"
    )
    # Families.tsv (only multi-member families)
    (tmp_path / 'families.tsv').write_text(
        "family_index\trepresentative_id\tn_members\n"
        "fam_000001\tpA1\t2\n"
    )


def _run(tmp_path, family_domtblouts=None, output_evalues=False,
        singleton_hits=None, paralog_cutoffs=None, **extra):
    args = [
        sys.executable, str(BIN),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--protein-map', str(tmp_path / 'protein_map.tsv'),
        '--config', str(tmp_path / 'config.csv'),
        '--output-matrix', str(tmp_path / 'matrix.tsv'),
        '--output-candidates', str(tmp_path / 'candidates.txt'),
    ]
    if family_domtblouts:
        args += ['--family-domtblout'] + [str(d) for d in family_domtblouts]
    if output_evalues:
        args += ['--output-evalues', str(tmp_path / 'evalues.tsv')]
    if singleton_hits:
        args += ['--singleton-hits'] + [str(p) for p in singleton_hits]
    if paralog_cutoffs:
        args += ['--paralog-cutoffs'] + [str(p) for p in paralog_cutoffs]
    for k, v in extra.items():
        args += [f'--{k}', str(v)]
    subprocess.run(args, check=True)
    matrix = pd.read_csv(tmp_path / 'matrix.tsv', sep='\t')
    cands = [c for c in (tmp_path / 'candidates.txt').read_text().splitlines() if c]
    if output_evalues:
        evalues = pd.read_csv(tmp_path / 'evalues.tsv', sep='\t', dtype=str, keep_default_na=False)
        return matrix, cands, evalues
    return matrix, cands


def test_matrix_columns_include_all_proteomes(tmp_path):
    _setup(tmp_path)
    matrix, _ = _run(tmp_path)
    assert list(matrix.columns) == ['protein_id', 'source_proteome', 'D1', 'D2', 'T1', 'T2']


def test_family_member_presence_does_not_leak_across_siblings(tmp_path):
    # Regression test (found via issue #29's integration test): with no hmmsearch hits at
    # all for family A, each member should be present ONLY in its own source proteome --
    # not also in its sibling's source proteome. A shared, mutated-in-place presence set
    # previously leaked pA1's source (T1) onto pA2's row once pA2 was processed second.
    _setup(tmp_path)
    matrix, _ = _run(tmp_path, family_domtblouts=[])
    row_a1 = matrix[matrix['protein_id'] == 'pA1'].iloc[0]
    row_a2 = matrix[matrix['protein_id'] == 'pA2'].iloc[0]
    assert row_a1['T1'] == 1 and row_a1['T2'] == 0
    assert row_a2['T2'] == 1 and row_a2['T1'] == 0


def test_family_present_in_targets_absent_from_disc_out(tmp_path):
    _setup(tmp_path)
    # Family A present in both targets (domtblout hits), absent from DISCOVERY_OUT
    dom_t1 = tmp_path / 'T1.domtblout'
    dom_t2 = tmp_path / 'T2.domtblout'
    _write_domtblout(dom_t1, [('pA1', 'pA1')])
    _write_domtblout(dom_t2, [('pA2', 'pA1')])

    matrix, cands = _run(tmp_path, family_domtblouts=[dom_t1, dom_t2])
    # pA1 should be present in T1 (source + domtblout) and T2 (domtblout)
    row = matrix[matrix['protein_id'] == 'pA1'].iloc[0]
    assert row['T1'] == 1 and row['T2'] == 1
    assert row['D1'] == 0 and row['D2'] == 0
    # Should be a novelty candidate
    assert 'T1::pA1' in cands


def test_output_evalues_sidecar_carries_family_hit_evalues(tmp_path):
    # issue #44: the family-HMM hit's full-seq e-value survives into the sidecar for
    # proteomes where the family was called present, and stays empty for the source
    # proteome (presence there is definitional, not a hit) and for absent proteomes.
    _setup(tmp_path)
    dom_t1 = tmp_path / 'T1.domtblout'
    dom_t2 = tmp_path / 'T2.domtblout'
    _write_domtblout(dom_t1, [('pA1', 'pA1', 1e-10)])
    _write_domtblout(dom_t2, [('pA2', 'pA1', 3.5e-8)])

    matrix, cands, evalues = _run(tmp_path, family_domtblouts=[dom_t1, dom_t2],
                                  output_evalues=True)
    row = evalues[evalues['protein_id'] == 'pA1'].iloc[0]
    assert row['T1'] == ''       # source proteome — not a hit
    assert row['T2'] == '3.5e-08'
    assert row['D1'] == '' and row['D2'] == ''  # absent — no hit


# --- Paralog-aware singleton filtering (issue #52) -------------------------------

HIT_HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'
PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def _write_hits(path, lines):
    path.write_text(HIT_HEADER + lines)


def test_singleton_present_via_flat_fallback_without_paralog_data(tmp_path):
    # No --paralog-cutoffs at all: falls back to --singleton-evalue for everyone,
    # identical to the pre-issue-#52 behaviour.
    _setup(tmp_path)
    hits = tmp_path / 'singletons_vs_D1.parsed.tsv'
    _write_hits(hits, 'pS1\td1_x\t1e-10\t100\tT1\tD1\n')
    matrix, cands = _run(tmp_path, singleton_hits=[hits])
    row = matrix[matrix['protein_id'] == 'pS1'].iloc[0]
    assert row['D1'] == 1


def test_singleton_paralog_cutoff_rejects_a_weak_hit(tmp_path):
    _setup(tmp_path)
    hits = tmp_path / 'singletons_vs_D1.parsed.tsv'
    _write_hits(hits, 'pS1\td1_x\t1e-3\t20\tT1\tD1\n')
    paralog = tmp_path / 'paralog_cutoffs.tsv'
    paralog.write_text(PARALOG_HEADER + 'pS1\tpEif1\t42\t1e-8\n')  # tighter than 1e-3
    matrix, cands = _run(tmp_path, singleton_hits=[hits], paralog_cutoffs=[paralog])
    row = matrix[matrix['protein_id'] == 'pS1'].iloc[0]
    assert row['D1'] == 0


def test_singleton_only_true_singletons_are_emitted_not_their_paralogs(tmp_path):
    # pEif1 is searched (as pS1's paralog, for the competition check) but is not itself
    # a singleton in cluster.tsv, so it must never appear as a matrix row.
    _setup(tmp_path)
    hits = tmp_path / 'singletons_vs_D1.parsed.tsv'
    _write_hits(hits, 'pS1\td1_x\t1e-10\t100\tT1\tD1\npEif1\td1_y\t1e-30\t150\tT1\tD1\n')
    matrix, cands = _run(tmp_path, singleton_hits=[hits])
    assert 'pEif1' not in set(matrix['protein_id'])


def test_singleton_competition_target_scope_keeps_hexa_like_ortholog(tmp_path):
    # Mirrors tests/test_build_presence_matrix.py's HEX-1/eIF5A case (and the real
    # NCU08332 finding from job 26997324): pS1's own hit (1e-69) beats its paralog
    # pEif1's hit to the *same* target (5e-12), so 'target' scope keeps the call even
    # though pEif1 hits harder elsewhere in D1 (1e-70).
    _setup(tmp_path)
    hits = tmp_path / 'singletons_vs_D1.parsed.tsv'
    _write_hits(hits, (
        'pS1\thex1\t1e-69\t230\tT1\tD1\n'
        'pEif1\teif2\t1e-70\t233\tT1\tD1\n'
        'pEif1\thex1\t5e-12\t45\tT1\tD1\n'
    ))
    paralog = tmp_path / 'paralog_cutoffs.tsv'
    paralog.write_text(PARALOG_HEADER + 'pS1\tpEif1\t42\t4.2e-11\n')
    matrix, cands = _run(tmp_path, singleton_hits=[hits], paralog_cutoffs=[paralog],
                         **{'paralog-competition-scope': 'target'})
    row = matrix[matrix['protein_id'] == 'pS1'].iloc[0]
    assert row['D1'] == 1


def test_singleton_competition_proteome_scope_drops_hexa_like_ortholog(tmp_path):
    _setup(tmp_path)
    hits = tmp_path / 'singletons_vs_D1.parsed.tsv'
    _write_hits(hits, (
        'pS1\thex1\t1e-69\t230\tT1\tD1\n'
        'pEif1\teif2\t1e-70\t233\tT1\tD1\n'
        'pEif1\thex1\t5e-12\t45\tT1\tD1\n'
    ))
    paralog = tmp_path / 'paralog_cutoffs.tsv'
    paralog.write_text(PARALOG_HEADER + 'pS1\tpEif1\t42\t4.2e-11\n')
    matrix, cands = _run(tmp_path, singleton_hits=[hits], paralog_cutoffs=[paralog],
                         **{'paralog-competition-scope': 'proteome'})
    row = matrix[matrix['protein_id'] == 'pS1'].iloc[0]
    assert row['D1'] == 0
