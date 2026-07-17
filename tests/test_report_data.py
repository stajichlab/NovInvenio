import json
import subprocess
import sys
from pathlib import Path

import pytest

from config_parser import parse_config
from report_data import (
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


def test_build_payload_rejects_a_matrix_with_no_config_columns(run_dir, samples):
    (run_dir / 'other.tsv').write_text('protein_id\tsource_proteome\tXxxx\nfoo\tXxxx\t1\n')
    with pytest.raises(ValueError, match='No proteome columns'):
        build_payload(run_dir / 'other.tsv', samples)


def test_make_report_writes_self_contained_html(run_dir):
    out = run_dir / 'report.html'
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
    out = run_dir / 'report.html'
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
