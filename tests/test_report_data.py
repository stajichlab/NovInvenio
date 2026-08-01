import json
import subprocess
import sys
from pathlib import Path

import pytest

from config_parser import parse_config
from report_data import (
    build_core_payload,
    build_losses_payload,
    build_payload,
    derive_novelties,
    gene_id_from_protein_id,
    read_fasta_seqs,
    read_tblastn_summary,
)

REPO = Path(__file__).resolve().parent.parent

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
OUT,Saccharomyces cerevisiae,S288c,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina
"""

# n1/n2 are novelties (present across the whole ingroup, absent from the outgroup).
# shared is present in an outgroup proteome; lonely fails the ingroup fraction.
MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions\tPfam_Evalues
n1\tNcra\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\tbZIP_1\tPF00170.27\t4.5e-09
n2\tAfum\t1\t1\t0\t0\t\t\t\t\t\t\t
shared\tNcra\t1\t1\t1\t0\t\tconserved thing\tPfam\tsp|P12345|TEST_YEAST Some protein\tAAA\tPF00004.31\t1e-20
lonely\tNcra\t1\t0\t0\t0\t\t\t\t\t\t\t
"""

# core1/core1b are present in every proteome (ingroup + outgroup) — independent
# ingroup-side discoveries of the same universally-conserved gene. near clears
# 3/4 columns (not core at the 0.95 default). low clears only 1/4. an_outgroup_row
# is present everywhere too, but sourced from an outgroup proteome — it must never
# appear in a CORE payload, since core genes are reported from the ingroup side only
# (mirrors presence_matrix.tsv's own scope restriction).
CORE_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions\tPfam_Evalues
core1\tNcra\t1\t1\t1\t1\tcore-gene\tessential enzyme\tPfam\t\tPF00001.1\tPF00001.1\t1e-30
core1b\tAfum\t1\t1\t1\t1\t\t\t\t\t\t\t
near\tNcra\t1\t1\t1\t0\t\t\t\t\t\t\t
low\tNcra\t1\t0\t0\t0\t\t\t\t\t\t\t
an_outgroup_row\tSpom\t1\t1\t1\t1\t\t\t\t\t\t\t
"""

# The loss matrix is the *full* presence table (LOSS_ANNOTATE annotates it, not the
# filtered candidate list), so build_losses_payload re-applies the loss predicate:
#   present in >= outgroup_min_frac of the outgroup AND <= loss_ingroup_max_frac of ingroup.
#   loss1 (Spom) + loss1b (Scer): the same gene recovered from two outgroup species,
#     conserved across the whole outgroup, absent from the ingroup — the strong candidate,
#     and a two-member gene family once clustered.
#   weak (Scer): present in only 1 of 2 outgroup species (frac 0.5) — dropped at the 0.75
#     default, kept only when outgroup_min_frac is lowered.
#   ingpres (Spom): conserved in the outgroup but still present in 1/2 of the ingroup —
#     dropped at the strict 0.0 default, kept only when loss_ingroup_max_frac is raised.
#   an_ingroup_row (Ncra): sourced from an ingroup proteome — must never appear (mirrors
#     CORE_MATRIX's opposite-direction check).
LOSSES_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions
loss1\tSpom\t0\t0\t1\t1\tERG-like\tsterol biosynthesis\tPfam\t\tp450\tPF00067.1
loss1b\tScer\t0\t0\t1\t1\t\t\t\t\t\t
weak\tScer\t0\t0\t0\t1\t\t\t\t\t\t
ingpres\tSpom\t1\t0\t1\t1\t\t\t\t\t\t
an_ingroup_row\tNcra\t1\t1\t1\t1\t\t\t\t\t\t
"""

# TBLASTN of the loss candidates against ingroup genomic DNA (Ncra, Afum columns) —
# loss1b gets a hit in Ncra, which should flag but not remove it.
LOSSES_TBLASTN = """\
protein_id\tNcra\tAfum
loss1\t0\t0
loss1b\t1\t0
"""

TBLASTN = """\
protein_id\tSpom\tScer
n1\t0\t0
n2\t0\t1
"""

FASTA = """\
>n1 some description here
MKVLAA
GGWT
>n2
MPPQQ
"""

# A second pathway's (mmseqs) presence matrix for the cross-method support column.
# n1 stays a novelty (concordant with the pairwise MATRIX); n2 gains an outgroup
# hit (Spom) so mmseqs does NOT call it novel → pairwise-only support.
SUPPORT_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer
n1\tNcra\t1\t1\t0\t0
n2\tAfum\t1\t1\t1\t0
"""


@pytest.fixture
def run_dir(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'tblastn.tsv').write_text(TBLASTN)
    (tmp_path / 'candidates.fa').write_text(FASTA)
    return tmp_path


@pytest.fixture
def samples(run_dir):
    return parse_config(run_dir / 'config.csv')


def payload_for(run_dir, samples, **kw):
    kw.setdefault('tblastn_path', run_dir / 'tblastn.tsv')
    kw.setdefault('candidates_fa', run_dir / 'candidates.fa')
    return build_payload(run_dir / 'matrix.tsv', samples, **kw)


def rows_by_id(payload):
    idx = payload['fields'].index('id')
    return {r[idx]: r for r in payload['rows']}


def test_gene_id_strips_transcript_suffixes():
    assert gene_id_from_protein_id('Afu1g01620-T-p1') == 'Afu1g01620'
    assert gene_id_from_protein_id('NCU00499-t26_1-p1') == 'NCU00499'
    assert gene_id_from_protein_id('AOBS_000123') == 'AOBS_000123'


def test_read_fasta_seqs_joins_wrapped_lines_and_cuts_id_at_whitespace(run_dir):
    seqs = read_fasta_seqs(run_dir / 'candidates.fa')
    assert seqs == {'n1': 'MKVLAAGGWT', 'n2': 'MPPQQ'}


def test_read_tblastn_summary_reports_only_genomes_with_hits(run_dir):
    genomes, hits = read_tblastn_summary(run_dir / 'tblastn.tsv')
    assert genomes == ['Spom', 'Scer']
    assert hits == {'n2': {'Scer'}}


def test_read_tblastn_summary_tolerates_missing_file():
    assert read_tblastn_summary(None) == ([], {})
    assert read_tblastn_summary('/nonexistent/tblastn.tsv') == ([], {})


def test_derive_novelties_applies_ingroup_fraction_and_outgroup_absence(run_dir):
    from report_data import read_matrix
    _, rows = read_matrix(run_dir / 'matrix.tsv')
    novel = derive_novelties(rows, ['Ncra', 'Afum'], ['Spom', 'Scer'], 0.75)
    # shared is present in Spom; lonely only hits 1/2 of the ingroup.
    assert novel == {'n1', 'n2'}


def test_derive_novelties_honours_a_looser_ingroup_fraction(run_dir):
    from report_data import read_matrix
    _, rows = read_matrix(run_dir / 'matrix.tsv')
    novel = derive_novelties(rows, ['Ncra', 'Afum'], ['Spom', 'Scer'], 0.5)
    assert novel == {'n1', 'n2', 'lonely'}


def test_payload_orders_proteomes_ingroup_first(run_dir, samples):
    payload = payload_for(run_dir, samples)
    assert [p['short'] for p in payload['proteomes']] == ['Ncra', 'Afum', 'Spom', 'Scer']
    assert [p['group'] for p in payload['proteomes']] == ['IN', 'IN', 'OUT', 'OUT']


def test_payload_presence_bitstring_follows_proteome_order(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert rows['shared'][F['pres']] == '1110'  # Ncra, Afum, Spom present; Scer absent
    assert rows['lonely'][F['pres']] == '1000'


def test_payload_evalues_align_with_presence_bitstring(run_dir, samples):
    # issue #44: 'ev' is a comma-separated list aligned position-for-position with
    # payload['proteomes'] / 'pres', empty where there's no e-value evidence.
    (run_dir / 'evalues.tsv').write_text(
        'protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\n'
        'n1\tNcra\t\t3.2e-40\t\t\n'
    )
    payload = payload_for(run_dir, samples, evalues_path=run_dir / 'evalues.tsv')
    assert payload['has_evalues'] is True
    row = rows_by_id(payload)['n1']
    ev_idx = payload['fields'].index('ev')
    shorts = [p['short'] for p in payload['proteomes']]
    assert row[ev_idx].split(',')[shorts.index('Afum')] == '3.2e-40'
    assert row[ev_idx].split(',')[shorts.index('Ncra')] == ''


def test_payload_has_evalues_false_without_evalues_path(run_dir, samples):
    payload = payload_for(run_dir, samples)
    assert payload['has_evalues'] is False
    row = rows_by_id(payload)['n1']
    ev_idx = payload['fields'].index('ev')
    assert row[ev_idx] == ',' * (len(payload['proteomes']) - 1)


def test_payload_appends_context_columns_without_affecting_novelty_stats(run_dir, samples):
    # issue #48: NEAR_INGROUP/BROAD_OUTGROUP context columns are appended after the
    # scored ingroup+outgroup columns, tagged {'context': True}, and 'pres'/'ev' extend
    # to cover them -- but they must never be counted as novelty-determining evidence
    # (that's still purely the strict IN/OUT matrix read_evalues/derive_novelties uses).
    context_samples_csv = CONFIG + (
        'NEAR_INGROUP,Near one,,Near1.pep.fa,Near1.dna.fa,Near1,Pezizomycotina\n'
        'BROAD_OUTGROUP,Broad one,,Broad1.pep.fa,Broad1.dna.fa,Broad1,Basidiomycota\n'
    )
    (run_dir / 'config.csv').write_text(context_samples_csv)
    context_samples = parse_config(run_dir / 'config.csv')
    (run_dir / 'context_matrix.tsv').write_text(
        'protein_id\tsource_proteome\tNear1\tBroad1\n'
        'n1\tNcra\t1\t0\n'
    )
    (run_dir / 'context_evalues.tsv').write_text(
        'protein_id\tsource_proteome\tNear1\tBroad1\n'
        'n1\tNcra\t2.1e-30\t\n'
    )
    payload = build_payload(
        run_dir / 'matrix.tsv', context_samples,
        tblastn_path=run_dir / 'tblastn.tsv', candidates_fa=run_dir / 'candidates.fa',
        context_matrix_path=run_dir / 'context_matrix.tsv',
        context_evalues_path=run_dir / 'context_evalues.tsv',
    )
    assert payload['has_context'] is True
    shorts = [p['short'] for p in payload['proteomes']]
    assert shorts[-2:] == ['Near1', 'Broad1']
    assert payload['proteomes'][-1]['context'] is True
    assert 'context' not in payload['proteomes'][0]  # scored proteomes stay untagged

    row = rows_by_id(payload)['n1']
    F = {n: i for i, n in enumerate(payload['fields'])}
    # 4 scored columns (Ncra, Afum, Spom, Scer) + 2 context columns (Near1, Broad1)
    assert row[F['pres']] == '1100' + '10'
    assert row[F['ev']].split(',')[-2:] == ['2.1e-30', '']


def test_payload_has_context_false_without_context_path(run_dir, samples):
    payload = payload_for(run_dir, samples)
    assert payload['has_context'] is False
    assert all('context' not in p for p in payload['proteomes'])


def test_payload_tblastn_bitstring_follows_genome_order(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert payload['tblastn_genomes'] == ['Spom', 'Scer']
    assert rows['n2'][F['tb']] == '01'
    assert rows['n1'][F['tb']] == '00'
    # A protein absent from tblastn_summary still gets a zero-filled bitstring.
    assert rows['shared'][F['tb']] == '00'


def test_payload_marks_novelties_when_no_novelty_tables_given(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert rows['n1'][F['nov']] == 1
    assert rows['n2'][F['nov']] == 1
    assert rows['shared'][F['nov']] == 0
    assert rows['lonely'][F['nov']] == 0


def test_payload_prefers_novelty_tables_over_derivation(run_dir, samples):
    # The table is authoritative even though `shared` would not be derived as a novelty.
    (run_dir / 'novelties.Ncra.tsv').write_text(
        'protein_id\tsource_proteome\tprotein_sequence\n'
        'shared\tNcra\tMMMM\n'
    )
    payload = payload_for(run_dir, samples, novelty_paths=[run_dir / 'novelties.Ncra.tsv'])
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert rows['shared'][F['nov']] == 1
    assert rows['n1'][F['nov']] == 0
    # Sequence comes from the novelties table, not candidates.fa.
    assert rows['shared'][F['seq']] == 'MMMM'


def test_support_single_method_labels_novelties_with_the_one_method(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert payload['methods'] == ['pairwise']
    assert payload['support_method'] is None
    # Novelty rows carry the sole method; non-novelties carry no support.
    assert rows['n1'][F['support']] == 'pairwise'
    assert rows['n2'][F['support']] == 'pairwise'
    assert rows['shared'][F['support']] == ''
    assert rows['lonely'][F['support']] == ''


def test_support_cross_method_concordance(run_dir, samples):
    (run_dir / 'support.tsv').write_text(SUPPORT_MATRIX)
    payload = payload_for(run_dir, samples,
                          support_matrix=run_dir / 'support.tsv',
                          support_method='mmseqs')
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert payload['methods'] == ['pairwise', 'mmseqs']
    assert payload['support_method'] == 'mmseqs'
    # n1 is novel in both pathways → concordant; n2 only in pairwise (mmseqs sees
    # an outgroup hit) → method-specific.
    assert rows['n1'][F['support']] == 'pairwise+mmseqs'
    assert rows['n2'][F['support']] == 'pairwise'
    assert rows['shared'][F['support']] == ''


def test_support_custom_primary_method_label(run_dir, samples):
    payload = payload_for(run_dir, samples, method='mmseqs')
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert payload['methods'] == ['mmseqs']
    assert rows['n1'][F['support']] == 'mmseqs'


def test_sequences_novelties_only_by_default(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert rows['n1'][F['seq']] == 'MKVLAAGGWT'
    assert rows['shared'][F['seq']] == ''  # in candidates.fa? no — and not a novelty anyway


def test_sequences_all_includes_non_novelties(run_dir, samples):
    (run_dir / 'candidates.fa').write_text(FASTA + '>shared\nMWWWW\n')
    payload = payload_for(run_dir, samples, sequences='all')
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert rows_by_id(payload)['shared'][F['seq']] == 'MWWWW'


def test_sequences_none_embeds_nothing(run_dir, samples):
    payload = payload_for(run_dir, samples, sequences='none')
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert all(r[F['seq']] == '' for r in payload['rows'])


def test_function_sources_are_interned(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)
    assert payload['fsources'] == ['ModelOrg_Ncra', 'Pfam']
    assert payload['fsources'][rows['n1'][F['fsrc']]] == 'ModelOrg_Ncra'
    assert rows['n2'][F['fsrc']] == -1  # no annotation


def test_payload_has_no_families_without_cluster_tsv(run_dir, samples):
    payload = payload_for(run_dir, samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert payload['families'] == []
    assert all(r[F['fam']] == -1 for r in payload['rows'])


def test_payload_groups_candidates_into_gene_families(run_dir, samples):
    # n1 (Ncra) and n2 (Afum) are independent hits of the same gene family
    # recovered in two different ingroup species; mmseqs clusters them together.
    (run_dir / 'clusters_cluster.tsv').write_text('n1\tn1\nn1\tn2\n')
    payload = payload_for(run_dir, samples, cluster_tsv=run_dir / 'clusters_cluster.tsv')
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = rows_by_id(payload)

    assert payload['families'] == [{'rep': 'n1', 'size': 2, 'species': ['Afum', 'Ncra']}]
    assert rows['n1'][F['fam']] == 0
    assert rows['n2'][F['fam']] == 0
    # shared/lonely never appear in candidates.fa's cluster, so they carry no family.
    assert rows['shared'][F['fam']] == -1
    assert rows['lonely'][F['fam']] == -1


def test_payload_ignores_singleton_clusters(run_dir, samples):
    # A cluster where the rep is its own only member has nothing to collapse.
    (run_dir / 'clusters_cluster.tsv').write_text('n1\tn1\n')
    payload = payload_for(run_dir, samples, cluster_tsv=run_dir / 'clusters_cluster.tsv')
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert payload['families'] == []
    assert rows_by_id(payload)['n1'][F['fam']] == -1


@pytest.fixture
def core_run_dir(run_dir):
    (run_dir / 'core_matrix.tsv').write_text(CORE_MATRIX)
    return run_dir


def core_rows_by_id(payload):
    idx = payload['fields'].index('id')
    return {r[idx]: r for r in payload['rows']}


def test_core_payload_keeps_only_rows_at_or_above_the_threshold(core_run_dir, samples):
    payload = build_core_payload(core_run_dir / 'core_matrix.tsv', samples)
    rows = core_rows_by_id(payload)
    assert set(rows) == {'core1', 'core1b'}


def test_core_payload_records_presence_fraction(core_run_dir, samples):
    payload = build_core_payload(core_run_dir / 'core_matrix.tsv', samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = core_rows_by_id(payload)
    assert rows['core1'][F['frac']] == 1.0


def test_core_payload_honours_a_looser_threshold(core_run_dir, samples):
    payload = build_core_payload(core_run_dir / 'core_matrix.tsv', samples, core_min_frac=0.7)
    assert set(core_rows_by_id(payload)) == {'core1', 'core1b', 'near'}


def test_core_payload_excludes_outgroup_sourced_rows_even_when_universal(core_run_dir, samples):
    # an_outgroup_row is present in every column but sourced from Spom (outgroup) —
    # core genes are only ever reported from the ingroup side.
    payload = build_core_payload(core_run_dir / 'core_matrix.tsv', samples, core_min_frac=0.5)
    assert 'an_outgroup_row' not in core_rows_by_id(payload)


def test_core_payload_groups_candidates_into_gene_families(core_run_dir, samples):
    (core_run_dir / 'clusters_cluster.tsv').write_text('core1\tcore1\ncore1\tcore1b\n')
    payload = build_core_payload(
        core_run_dir / 'core_matrix.tsv', samples,
        cluster_tsv=core_run_dir / 'clusters_cluster.tsv',
    )
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = core_rows_by_id(payload)
    assert payload['families'] == [{'rep': 'core1', 'size': 2, 'species': ['Afum', 'Ncra']}]
    assert rows['core1'][F['fam']] == 0
    assert rows['core1b'][F['fam']] == 0


def test_core_payload_carries_annotation_columns(core_run_dir, samples):
    payload = build_core_payload(core_run_dir / 'core_matrix.tsv', samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = core_rows_by_id(payload)
    assert rows['core1'][F['gene']] == 'core-gene'
    assert rows['core1'][F['prod']] == 'essential enzyme'
    assert payload['fsources'][rows['core1'][F['fsrc']]] == 'Pfam'
    assert rows['core1'][F['pfam_n']] == 'PF00001.1'


def test_core_payload_rejects_a_matrix_with_no_config_columns(run_dir, samples):
    (run_dir / 'other.tsv').write_text('protein_id\tsource_proteome\tXxxx\nfoo\tXxxx\t1\n')
    with pytest.raises(ValueError, match='No proteome columns'):
        build_core_payload(run_dir / 'other.tsv', samples)


def test_make_core_report_writes_self_contained_html(core_run_dir):
    out = core_run_dir / 'core.html'
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_core_report.py'),
         '--matrix', str(core_run_dir / 'core_matrix.tsv'),
         '--config', str(core_run_dir / 'config.csv'),
         '--project', 'test_run',
         '--output', str(out)],
        check=True, capture_output=True,
    )
    doc = out.read_text()
    assert doc.startswith('<!doctype html>')
    assert 'test_run' in doc
    for pattern in ('src="http', 'href="http://cdn', '<link rel="stylesheet"', 'import('):
        assert pattern not in doc

    payload = json.loads(doc.split('id="payload">')[1].split('</script>')[0])
    assert payload['project'] == 'test_run'
    assert len(payload['rows']) == 2


@pytest.fixture
def losses_run_dir(run_dir):
    (run_dir / 'losses_matrix.tsv').write_text(LOSSES_MATRIX)
    (run_dir / 'losses_tblastn.tsv').write_text(LOSSES_TBLASTN)
    return run_dir


def losses_rows_by_id(payload):
    idx = payload['fields'].index('id')
    return {r[idx]: r for r in payload['rows']}


def test_losses_payload_keeps_only_conserved_ingroup_absent_rows(losses_run_dir, samples):
    # Default thresholds (outgroup_min_frac 0.75, loss_ingroup_max_frac 0.0): the two
    # whole-outgroup, ingroup-absent rows survive; the single-outgroup 'weak', the
    # ingroup-present 'ingpres', and the ingroup-sourced row are all filtered out.
    payload = build_losses_payload(losses_run_dir / 'losses_matrix.tsv', samples)
    assert set(losses_rows_by_id(payload)) == {'loss1', 'loss1b'}


def test_losses_payload_outgroup_min_frac_admits_narrower_candidates(losses_run_dir, samples):
    # Lowering the outgroup conservation threshold lets 'weak' (1 of 2 outgroup species) in.
    payload = build_losses_payload(
        losses_run_dir / 'losses_matrix.tsv', samples, outgroup_min_frac=0.5,
    )
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = losses_rows_by_id(payload)
    assert 'weak' in rows
    assert rows['weak'][F['frac']] == 0.5
    assert rows['loss1'][F['frac']] == 1.0


def test_losses_payload_loss_ingroup_max_frac_admits_nearly_missing(losses_run_dir, samples):
    # 'ingpres' still survives in 1 of 2 ingroup species (frac 0.5): excluded at the
    # strict 0.0 default, admitted once loss_ingroup_max_frac allows half the ingroup.
    strict = build_losses_payload(losses_run_dir / 'losses_matrix.tsv', samples)
    assert 'ingpres' not in losses_rows_by_id(strict)

    relaxed = build_losses_payload(
        losses_run_dir / 'losses_matrix.tsv', samples, loss_ingroup_max_frac=0.5,
    )
    F = {n: i for i, n in enumerate(relaxed['fields'])}
    rows = losses_rows_by_id(relaxed)
    assert 'ingpres' in rows
    assert rows['ingpres'][F['in_retained']] == 1   # retained in one ingroup species


def test_losses_payload_flags_but_keeps_tblastn_hits(losses_run_dir, samples):
    # Reporting-only, mirroring make_novelties.py's --skip_tblastn_filter default:
    # a TBLASTN hit downgrades priority, it never removes the row.
    payload = build_losses_payload(
        losses_run_dir / 'losses_matrix.tsv', samples,
        tblastn_path=losses_run_dir / 'losses_tblastn.tsv',
    )
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = losses_rows_by_id(payload)
    assert rows['loss1'][F['tb_hit']] == 0
    assert rows['loss1'][F['tb_genomes']] == ''
    assert rows['loss1b'][F['tb_hit']] == 1
    assert rows['loss1b'][F['tb_genomes']] == 'Ncra'


def test_losses_payload_without_tblastn_path_has_no_hits(losses_run_dir, samples):
    payload = build_losses_payload(losses_run_dir / 'losses_matrix.tsv', samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert all(r[F['tb_hit']] == 0 for r in payload['rows'])


def test_losses_payload_singleton_breadth_is_its_own_presence(losses_run_dir, samples):
    # With no cluster file every kept row is a singleton: out_breadth is the row's own
    # outgroup presence count, in_retained its own ingroup presence count.
    payload = build_losses_payload(losses_run_dir / 'losses_matrix.tsv', samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = losses_rows_by_id(payload)
    assert rows['loss1'][F['out_breadth']] == 2   # present in both outgroup species
    assert rows['loss1'][F['in_retained']] == 0


def test_losses_payload_groups_candidates_into_gene_families(losses_run_dir, samples):
    # loss1 + loss1b are the same gene recovered from two outgroup species.
    (losses_run_dir / 'loss_clusters_cluster.tsv').write_text('loss1\tloss1\nloss1\tloss1b\n')
    payload = build_losses_payload(
        losses_run_dir / 'losses_matrix.tsv', samples,
        cluster_tsv=losses_run_dir / 'loss_clusters_cluster.tsv',
    )
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = losses_rows_by_id(payload)
    assert payload['families'] == [{'rep': 'loss1', 'size': 2, 'species': ['Scer', 'Spom']}]
    assert rows['loss1'][F['fam']] == 0
    assert rows['loss1b'][F['fam']] == 0
    # Family-level aggregates span every member: both outgroup species, no ingroup retention.
    assert rows['loss1'][F['out_breadth']] == 2
    assert rows['loss1b'][F['out_breadth']] == 2
    assert rows['loss1'][F['in_retained']] == 0


def test_losses_payload_carries_annotation_columns(losses_run_dir, samples):
    payload = build_losses_payload(losses_run_dir / 'losses_matrix.tsv', samples)
    F = {n: i for i, n in enumerate(payload['fields'])}
    rows = losses_rows_by_id(payload)
    assert rows['loss1'][F['gene']] == 'ERG-like'
    assert rows['loss1'][F['prod']] == 'sterol biosynthesis'
    assert payload['fsources'][rows['loss1'][F['fsrc']]] == 'Pfam'
    assert rows['loss1'][F['pfam_n']] == 'p450'


def test_losses_payload_rejects_a_matrix_with_no_config_columns(run_dir, samples):
    (run_dir / 'other.tsv').write_text('protein_id\tsource_proteome\tXxxx\nfoo\tXxxx\t1\n')
    with pytest.raises(ValueError, match='No proteome columns'):
        build_losses_payload(run_dir / 'other.tsv', samples)


def test_make_losses_report_writes_self_contained_html(losses_run_dir):
    out = losses_run_dir / 'losses.html'
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_losses_report.py'),
         '--matrix', str(losses_run_dir / 'losses_matrix.tsv'),
         '--config', str(losses_run_dir / 'config.csv'),
         '--tblastn_summary', str(losses_run_dir / 'losses_tblastn.tsv'),
         '--project', 'test_run',
         '--output', str(out)],
        check=True, capture_output=True,
    )
    doc = out.read_text()
    assert doc.startswith('<!doctype html>')
    assert 'test_run' in doc
    for pattern in ('src="http', 'href="http://cdn', '<link rel="stylesheet"', 'import('):
        assert pattern not in doc

    payload = json.loads(doc.split('id="payload">')[1].split('</script>')[0])
    assert payload['project'] == 'test_run'
    assert len(payload['rows']) == 2


def test_build_payload_rejects_a_matrix_with_no_config_columns(run_dir, samples):
    (run_dir / 'other.tsv').write_text('protein_id\tsource_proteome\tXxxx\nfoo\tXxxx\t1\n')
    with pytest.raises(ValueError, match='No proteome columns'):
        build_payload(run_dir / 'other.tsv', samples)


def test_make_report_writes_self_contained_html(run_dir):
    out = run_dir / 'novelties.html'
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_report.py'),
         '--matrix', str(run_dir / 'matrix.tsv'),
         '--config', str(run_dir / 'config.csv'),
         '--tblastn_summary', str(run_dir / 'tblastn.tsv'),
         '--candidates_fa', str(run_dir / 'candidates.fa'),
         '--project', 'test_run',
         '--output', str(out)],
        check=True, capture_output=True,
    )
    doc = out.read_text()
    assert doc.startswith('<!doctype html>')
    assert 'test_run' in doc
    # Self-contained: nothing may be fetched at open time.
    for pattern in ('src="http', 'href="http://cdn', '<link rel="stylesheet"', 'import('):
        assert pattern not in doc

    payload = json.loads(doc.split('id="payload">')[1].split('</script>')[0])
    assert payload['project'] == 'test_run'
    assert len(payload['rows']) == 4


def test_make_report_escapes_a_script_tag_hiding_in_an_annotation(run_dir):
    # Product descriptions come from SwissProt/Pfam text and are untrusted.
    (run_dir / 'matrix.tsv').write_text(
        MATRIX + 'evil\tNcra\t1\t1\t0\t0\t\t</script><script>alert(1)</script>\tPfam\t\t\t\t\n'
    )
    out = run_dir / 'novelties.html'
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_report.py'),
         '--matrix', str(run_dir / 'matrix.tsv'),
         '--config', str(run_dir / 'config.csv'),
         '--output', str(out)],
        check=True, capture_output=True,
    )
    doc = out.read_text()
    payload_block = doc.split('id="payload">')[1].split('</script>')[0]
    # The only sequence that can terminate a <script> block early is "</script";
    # a bare "<script>" inside it is inert. So the invariant is that no literal
    # closing tag survives into the block, and the payload still parses whole.
    assert '</script' not in payload_block.lower()
    payload = json.loads(payload_block)
    assert len(payload['rows']) == 5

    # The description must survive intact as data — the browser inserts it with
    # textContent, so it renders as literal text rather than markup.
    prod = payload['fields'].index('prod')
    descriptions = [r[prod] for r in payload['rows']]
    assert '</script><script>alert(1)</script>' in descriptions


# novelty_discovery (todo/novelty-discovery-screen.md) configs use DISCOVERY_TARGET/
# DISCOVERY_OUT instead of IN/OUT — the report payload builders must treat them identically
# (issue #33; labels renamed for clarity in todo/rename-novelty-discovery-group-labels.md).
TARGET_DISC_OUT_CONFIG = (
    CONFIG.replace('IN,', 'DISCOVERY_TARGET,').replace('OUT,', 'DISCOVERY_OUT,')
)


def test_payload_treats_target_disc_out_roles_like_in_out(tmp_path):
    (tmp_path / 'config.csv').write_text(TARGET_DISC_OUT_CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    samples = parse_config(tmp_path / 'config.csv')

    payload = build_payload(tmp_path / 'matrix.tsv', samples, sequences='none')

    assert [p['short'] for p in payload['proteomes']] == ['Ncra', 'Afum', 'Spom', 'Scer']
    assert [p['group'] for p in payload['proteomes']] == [
        'DISCOVERY_TARGET', 'DISCOVERY_TARGET', 'DISCOVERY_OUT', 'DISCOVERY_OUT',
    ]
    rows = rows_by_id(payload)
    F = {n: i for i, n in enumerate(payload['fields'])}
    # Same novelty calls as the IN/OUT config: n1/n2 novel, shared/lonely not.
    assert {pid for pid, r in rows.items() if r[F['nov']] == 1} == {'n1', 'n2'}


def test_core_payload_treats_target_disc_out_roles_like_in_out(tmp_path):
    (tmp_path / 'config.csv').write_text(TARGET_DISC_OUT_CONFIG)
    (tmp_path / 'matrix.tsv').write_text(CORE_MATRIX)
    samples = parse_config(tmp_path / 'config.csv')

    payload = build_core_payload(tmp_path / 'matrix.tsv', samples, core_min_frac=0.95)

    ids = {r[payload['fields'].index('id')] for r in payload['rows']}
    assert ids == {'core1', 'core1b'}


# novelty_category (issue #27's bin/novelty_screen.py output, issue #28's report rendering).
# Rows without the column (pairwise/mmseqs runs) must still work -- 'category' defaults to ''.
CATEGORY_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tnovelty_category
n1\tNcra\t1\t1\t0\t0\ttarget_specific
n2\tAfum\t1\t1\t0\t0\tclade_specific
shared\tNcra\t1\t1\t1\t0\tfalse_novelty
lonely\tNcra\t1\t0\t0\t0\t
"""


def test_payload_includes_novelty_category_field(run_dir, samples):
    (run_dir / 'matrix.tsv').write_text(CATEGORY_MATRIX)
    payload = payload_for(run_dir, samples)

    assert 'category' in payload['fields']
    rows = rows_by_id(payload)
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert rows['n1'][F['category']] == 'target_specific'
    assert rows['n2'][F['category']] == 'clade_specific'
    assert rows['shared'][F['category']] == 'false_novelty'
    assert rows['lonely'][F['category']] == ''
    assert payload['novelty_categories'] == ['clade_specific', 'false_novelty', 'target_specific']


def test_payload_defaults_category_to_empty_when_column_absent(run_dir, samples):
    # MATRIX (the module-level default fixture) has no novelty_category column at all --
    # pairwise/mmseqs runs never populate it.
    payload = payload_for(run_dir, samples)

    assert 'category' in payload['fields']
    F = {n: i for i, n in enumerate(payload['fields'])}
    assert all(r[F['category']] == '' for r in payload['rows'])
    assert payload['novelty_categories'] == []


def test_make_report_escapes_a_script_tag_hiding_in_novelty_category(run_dir):
    # novelty_category is always one of three fixed enum strings written by
    # bin/novelty_screen.py, never free text from an external source -- but the payload's
    # escaping is applied to the whole embedded JSON block (see
    # test_make_report_escapes_a_script_tag_hiding_in_an_annotation), so an adversarial value
    # here must be neutralised the same way any other field's would be.
    matrix_with_category = (
        "protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description"
        "\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions\tPfam_Evalues"
        "\tnovelty_category\n"
        "n1\tNcra\t1\t1\t0\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\tbZIP_1"
        "\tPF00170.27\t4.5e-09\t</script><script>alert(1)</script>\n"
    )
    (run_dir / 'matrix.tsv').write_text(matrix_with_category)
    out = run_dir / 'novelties.html'
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_report.py'),
         '--matrix', str(run_dir / 'matrix.tsv'),
         '--config', str(run_dir / 'config.csv'),
         '--output', str(out)],
        check=True, capture_output=True,
    )
    doc = out.read_text()
    payload_block = doc.split('id="payload">')[1].split('</script>')[0]
    assert '</script' not in payload_block.lower()
    payload = json.loads(payload_block)
    cat_idx = payload['fields'].index('category')
    categories = [r[cat_idx] for r in payload['rows']]
    assert '</script><script>alert(1)</script>' in categories
