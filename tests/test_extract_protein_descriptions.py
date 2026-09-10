import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from extract_protein_descriptions import build_descriptions, parse_fasta_headers  # noqa: E402


def test_parses_uniprot_style_header_with_gene_name():
    lines = iter([
        '>tr|A7UWL5|A7UWL5_NEUCR NRS/ER OS=Neurospora crassa (strain X) OX=367110 GN=NCU10683 PE=4 SV=1\n',
        'MKVLA\n',
    ])
    [(pid, gene_name, description)] = list(parse_fasta_headers(lines))
    assert pid == 'tr|A7UWL5|A7UWL5_NEUCR'
    assert gene_name == 'NCU10683'
    assert description == 'NRS/ER'


def test_parses_header_with_no_gene_name():
    lines = iter(['>tr|Q1|Q1_FOO Uncharacterized protein OS=Foo bar OX=1 PE=4 SV=1\n'])
    [(pid, gene_name, description)] = list(parse_fasta_headers(lines))
    assert pid == 'tr|Q1|Q1_FOO'
    assert gene_name == ''
    assert description == 'Uncharacterized protein'


def test_non_uniprot_header_falls_back_to_whole_remainder():
    lines = iter(['>contig_1 some plain description\n'])
    [(pid, gene_name, description)] = list(parse_fasta_headers(lines))
    assert pid == 'contig_1'
    assert gene_name == ''
    assert description == 'some plain description'


def test_header_with_no_description_at_all():
    lines = iter(['>bareID\n'])
    [(pid, gene_name, description)] = list(parse_fasta_headers(lines))
    assert pid == 'bareID'
    assert gene_name == ''
    assert description == ''


def test_build_descriptions_first_occurrence_wins_across_files(tmp_path):
    fa1 = tmp_path / 'a.pep.fa'
    fa1.write_text('>id1 First description OS=Foo OX=1\nMKV\n')
    fa2 = tmp_path / 'b.pep.fa'
    fa2.write_text('>id1 Second description OS=Foo OX=1\nMKV\n>id2 Other OS=Bar OX=2\nMKV\n')

    descriptions = build_descriptions([fa1, fa2])

    assert descriptions['id1'] == ('', 'First description')
    assert descriptions['id2'] == ('', 'Other')
