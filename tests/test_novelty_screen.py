"""Unit tests for bin/novelty_screen.py (issue #27 — novelty_screen three-category
refinement, todo/novelty-discovery-screen.md)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parent.parent / 'bin' / 'novelty_screen.py'

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
DISCOVERY_TARGET,Species one,,target1.pep.fa,target1.dna.fa,T1,X
DISCOVERY_TARGET,Species two,,target2.pep.fa,target2.dna.fa,T2,X
DISCOVERY_OUT,Outgroup one,,disc1.pep.fa,disc1.dna.fa,D1,Y
NEAR_INGROUP,Near relative,,near1.pep.fa,near1.dna.fa,N1,X
BROAD_OUTGROUP,Distant one,,broad1.pep.fa,broad1.dna.fa,B1,Z
BROAD_OUTGROUP,Distant two,,broad2.pep.fa,broad2.dna.fa,B2,Z
"""

# Discovery matrix: three phase-1 candidates (famA -> pA1/pA2, famB -> pB1/pB2, famC -> pC1)
# plus one non-candidate row (other) that should be carried through with an empty category.
DISCOVERY_MATRIX = """\
protein_id\tsource_proteome\tD1\tT1\tT2
pA1\tT1\t0\t1\t1
pA2\tT2\t0\t1\t1
pB1\tT1\t0\t1\t1
pB2\tT2\t0\t1\t1
pC1\tT1\t0\t1\t0
other\tT1\t1\t1\t0
"""

DISCOVERY_CANDIDATES = "T1::pA1\nT1::pB1\nT1::pC1\n"

CLUSTER_TSV = (
    "pA1\tpA1\npA1\tpA2\n"
    "pB1\tpB1\npB1\tpB2\n"
    "pC1\tpC1\n"
)


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
    (tmp_path / 'discovery_matrix.tsv').write_text(DISCOVERY_MATRIX)
    (tmp_path / 'discovery_candidates.txt').write_text(DISCOVERY_CANDIDATES)
    (tmp_path / 'cluster.tsv').write_text(CLUSTER_TSV)


def _run(tmp_path, near_in_domtblouts=None, broad_out_domtblouts=None,
        near_in_singleton_hits=None, broad_out_singleton_hits=None,
        paralog_cutoffs=None, **extra):
    args = [
        sys.executable, str(BIN),
        '--discovery-matrix', str(tmp_path / 'discovery_matrix.tsv'),
        '--discovery-candidates', str(tmp_path / 'discovery_candidates.txt'),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--config', str(tmp_path / 'config.csv'),
        '--output-matrix', str(tmp_path / 'screened_matrix.tsv'),
        '--output-candidates', str(tmp_path / 'screened_candidates.txt'),
    ]
    if near_in_domtblouts:
        args += ['--near-in-domtblout'] + [str(d) for d in near_in_domtblouts]
    if broad_out_domtblouts:
        args += ['--broad-out-domtblout'] + [str(d) for d in broad_out_domtblouts]
    if near_in_singleton_hits:
        args += ['--near-in-singleton-hits'] + [str(p) for p in near_in_singleton_hits]
    if broad_out_singleton_hits:
        args += ['--broad-out-singleton-hits'] + [str(p) for p in broad_out_singleton_hits]
    if paralog_cutoffs:
        args += ['--paralog-cutoffs'] + [str(p) for p in paralog_cutoffs]
    for k, v in extra.items():
        args += [f'--{k}', str(v)]
    subprocess.run(args, check=True, capture_output=True, text=True)
    # keep_default_na=False: an empty novelty_category is a real '' value, not a missing one.
    matrix = pd.read_csv(tmp_path / 'screened_matrix.tsv', sep='\t', keep_default_na=False)
    cands = [c for c in (tmp_path / 'screened_candidates.txt').read_text().splitlines() if c]
    return matrix, cands


def test_matrix_adds_near_in_broad_out_columns_and_category(tmp_path):
    _setup(tmp_path)
    matrix, _ = _run(tmp_path)
    assert list(matrix.columns) == [
        'protein_id', 'source_proteome', 'D1', 'T1', 'T2', 'N1', 'B1', 'B2', 'novelty_category',
    ]


def test_no_hits_anywhere_is_target_specific(tmp_path):
    _setup(tmp_path)
    matrix, cands = _run(tmp_path)
    row = matrix[matrix['protein_id'] == 'pC1'].iloc[0]
    assert row['novelty_category'] == 'target_specific'
    assert row['N1'] == 0 and row['B1'] == 0 and row['B2'] == 0
    assert 'T1::pC1' in cands


def test_hit_in_near_in_only_is_clade_specific(tmp_path):
    _setup(tmp_path)
    near_dom = tmp_path / 'N1.family.domtblout'
    _write_domtblout(near_dom, [('protN1', 'pA1')])

    matrix, cands = _run(tmp_path, near_in_domtblouts=[near_dom])
    row = matrix[matrix['protein_id'] == 'pA1'].iloc[0]
    assert row['novelty_category'] == 'clade_specific'
    assert row['N1'] == 1
    assert row['B1'] == 0 and row['B2'] == 0
    assert 'T1::pA1' in cands
    # A different family member (pA2) is not itself the family rep but shares pA1's category
    # via the shared family HMM -- both fall under the same clade_specific classification.
    row2 = matrix[matrix['protein_id'] == 'pA2'].iloc[0]
    assert row2['novelty_category'] == ''  # pA2 was never itself a phase-1 candidate


def test_hit_in_broad_out_is_false_novelty_and_removed_from_candidates(tmp_path):
    _setup(tmp_path)
    broad_dom = tmp_path / 'B1.family.domtblout'
    _write_domtblout(broad_dom, [('protB1', 'pB1')])

    matrix, cands = _run(tmp_path, broad_out_domtblouts=[broad_dom])
    row = matrix[matrix['protein_id'] == 'pB1'].iloc[0]
    assert row['novelty_category'] == 'false_novelty'
    assert row['B1'] == 1
    # Removed from the screened candidate list even though it's still in the matrix.
    assert 'T1::pB1' not in cands
    assert 'T1::pA1' in cands and 'T1::pC1' in cands


def test_hit_in_both_near_in_and_broad_out_is_still_false_novelty(tmp_path):
    _setup(tmp_path)
    near_dom = tmp_path / 'N1.family.domtblout'
    broad_dom = tmp_path / 'B2.family.domtblout'
    _write_domtblout(near_dom, [('protN1', 'pB1')])
    _write_domtblout(broad_dom, [('protB2', 'pB1')])

    matrix, cands = _run(tmp_path, near_in_domtblouts=[near_dom], broad_out_domtblouts=[broad_dom])
    row = matrix[matrix['protein_id'] == 'pB1'].iloc[0]
    assert row['novelty_category'] == 'false_novelty'
    assert 'T1::pB1' not in cands


def test_non_candidate_rows_carry_empty_category_but_real_presence(tmp_path):
    _setup(tmp_path)
    near_dom = tmp_path / 'N1.family.domtblout'
    _write_domtblout(near_dom, [('protN1', 'other')])

    matrix, cands = _run(tmp_path, near_in_domtblouts=[near_dom])
    row = matrix[matrix['protein_id'] == 'other'].iloc[0]
    assert row['novelty_category'] == ''
    assert row['N1'] == 1
    assert not any(c.endswith('::other') for c in cands)


def test_preexisting_near_in_broad_out_columns_are_updated_not_duplicated(tmp_path):
    # bin/novelty_presence_matrix.py's own header spans every proteome short ID in the
    # config (not just DISCOVERY_TARGET/DISCOVERY_OUT), so a discovery matrix commonly
    # already carries NEAR_INGROUP/BROAD_OUTGROUP columns, always 0 (phase 1 never searches
    # them). novelty_screen.py must overwrite those columns in place rather than appending
    # duplicates.
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'discovery_matrix.tsv').write_text(
        "protein_id\tsource_proteome\tD1\tT1\tT2\tN1\tB1\tB2\n"
        "pA1\tT1\t0\t1\t1\t0\t0\t0\n"
    )
    (tmp_path / 'discovery_candidates.txt').write_text("T1::pA1\n")
    (tmp_path / 'cluster.tsv').write_text("pA1\tpA1\n")

    broad_dom = tmp_path / 'B1.family.domtblout'
    _write_domtblout(broad_dom, [('protB1', 'pA1')])

    matrix, cands = _run(tmp_path, broad_out_domtblouts=[broad_dom])
    assert list(matrix.columns) == [
        'protein_id', 'source_proteome', 'D1', 'T1', 'T2', 'N1', 'B1', 'B2', 'novelty_category',
    ]
    row = matrix[matrix['protein_id'] == 'pA1'].iloc[0]
    assert row['B1'] == 1
    assert row['novelty_category'] == 'false_novelty'
    assert 'T1::pA1' not in cands


def test_family_threshold_gates_presence(tmp_path):
    _setup(tmp_path)
    (tmp_path / 'thresholds.tsv').write_text(
        "rep_id\tthreshold_evalue\npA1\t1e-20\n"
    )
    # Hit E-value (1e-10) is worse (higher) than the calibrated threshold (1e-20) -> not present.
    broad_dom = tmp_path / 'B1.family.domtblout'
    _write_domtblout(broad_dom, [('protB1', 'pA1', 1e-10)])

    matrix, cands = _run(tmp_path, broad_out_domtblouts=[broad_dom],
                         **{'family-thresholds': tmp_path / 'thresholds.tsv'})
    row = matrix[matrix['protein_id'] == 'pA1'].iloc[0]
    assert row['B1'] == 0
    assert row['novelty_category'] == 'target_specific'
    assert 'T1::pA1' in cands


# --- Singleton screening (issue #52) ----------------------------------------------
# pC1 is a singleton in CLUSTER_TSV ('pC1\tpC1', no family) -- family HMM search can
# never see it, so before this feature it always defaulted to target_specific
# regardless of true NEAR_INGROUP/BROAD_OUTGROUP presence. These tests search it
# directly, the same way the pairwise singleton search does.

HIT_HEADER = 'query_id\ttarget_id\tevalue\tbitscore\tquery_proteome\ttarget_proteome\n'
PARALOG_HEADER = 'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'


def _write_hits(path, lines):
    path.write_text(HIT_HEADER + lines)


def test_singleton_hit_in_near_in_reclassifies_clade_specific(tmp_path):
    _setup(tmp_path)
    near_hits = tmp_path / 'singletons_vs_N1.parsed.tsv'
    _write_hits(near_hits, 'pC1\tn1_x\t1e-10\t100\tT1\tN1\n')

    matrix, cands = _run(tmp_path, near_in_singleton_hits=[near_hits])
    row = matrix[matrix['protein_id'] == 'pC1'].iloc[0]
    assert row['N1'] == 1
    assert row['novelty_category'] == 'clade_specific'
    assert 'T1::pC1' in cands


def test_singleton_hit_in_broad_out_reclassifies_false_novelty(tmp_path):
    _setup(tmp_path)
    broad_hits = tmp_path / 'singletons_vs_B1.parsed.tsv'
    _write_hits(broad_hits, 'pC1\tb1_x\t1e-10\t100\tT1\tB1\n')

    matrix, cands = _run(tmp_path, broad_out_singleton_hits=[broad_hits])
    row = matrix[matrix['protein_id'] == 'pC1'].iloc[0]
    assert row['B1'] == 1
    assert row['novelty_category'] == 'false_novelty'
    assert 'T1::pC1' not in cands


def test_singleton_paralog_competition_prevents_a_false_novelty_call(tmp_path):
    # The NCU08332/HEX-1-vs-eIF5A pattern: pC1's own hit to a real (weak) ortholog
    # target loses to its paralog's much stronger hit in the SAME broad-outgroup
    # proteome -- 'proteome' scope disqualifies it, so pC1 is correctly NOT demoted to
    # false_novelty by what is actually cross-reactivity with a conserved paralog.
    _setup(tmp_path)
    broad_hits = tmp_path / 'singletons_vs_B1.parsed.tsv'
    _write_hits(broad_hits, (
        'pC1\thex1\t5e-12\t45\tT1\tB1\n'
        'pParalog\teif2\t1e-70\t233\tT1\tB1\n'
    ))
    paralog = tmp_path / 'paralog_cutoffs.tsv'
    paralog.write_text(PARALOG_HEADER + 'pC1\tpParalog\t42\t4.2e-11\n')

    matrix, cands = _run(tmp_path, broad_out_singleton_hits=[broad_hits],
                         paralog_cutoffs=[paralog])
    row = matrix[matrix['protein_id'] == 'pC1'].iloc[0]
    assert row['B1'] == 0
    assert row['novelty_category'] == 'target_specific'
    assert 'T1::pC1' in cands


def test_singleton_paralog_added_to_search_is_never_itself_a_row(tmp_path):
    # pParalog is searched (as pC1's paralog, for the competition check) but is not
    # itself a singleton in cluster.tsv, so it must never appear as a matrix row.
    _setup(tmp_path)
    broad_hits = tmp_path / 'singletons_vs_B1.parsed.tsv'
    _write_hits(broad_hits, (
        'pC1\thex1\t1e-69\t230\tT1\tB1\n'
        'pParalog\teif2\t1e-70\t233\tT1\tB1\n'
    ))
    matrix, cands = _run(tmp_path, broad_out_singleton_hits=[broad_hits])
    assert 'pParalog' not in set(matrix['protein_id'])
