import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))
from collapse_isoforms import collapse  # noqa: E402

FEATURE_TABLE = """\
# feature\tclass\tassembly\tassembly_unit\tseq_type\tchromosome\tgenomic_accession\tstart\tend\tstrand\tproduct_accession\tnon-redundant_refseq\trelated_accession\tname\tsymbol\tGeneID\tlocus_tag\tfeature_interval_length\tproduct_length\tattributes
gene\tprotein_coding\tGCF_1\tPrimary Assembly\tchromosome\tI\tNC_1\t1\t100\t+\t\t\t\thomt-1\thomt-1\t171590\tCELE_1\t100\t\t
CDS\twith_protein\tGCF_1\tPrimary Assembly\tchromosome\tI\tNC_1\t1\t50\t+\tNP_001.1\t\tNM_001.1\thomt-1\thomt-1\t171590\tCELE_1\t50\t100\t
CDS\twith_protein\tGCF_1\tPrimary Assembly\tchromosome\tI\tNC_1\t1\t80\t+\tNP_002.1\t\tNM_002.1\thomt-1 isoform B\thomt-1\t171590\tCELE_1\t80\t160\t
CDS\twith_protein\tGCF_1\tPrimary Assembly\tchromosome\tI\tNC_1\t200\t300\t+\tNP_003.1\t\tNM_003.1\tnlp-40\tnlp-40\t171591\tCELE_2\t100\t123\t
"""

PROTEIN_FASTA = """\
>NP_001.1 homt-1 [Species]
MAAA
>NP_002.1 homt-1 isoform B [Species]
MAAAAAAA
>NP_003.1 nlp-40 [Species]
MBBB
>NP_999.1 orphan protein with no feature-table entry [Species]
MCCC
"""


def test_collapse_keeps_longest_isoform_per_gene(tmp_path):
    ft = tmp_path / 'feature_table.txt'
    ft.write_text(FEATURE_TABLE)
    fasta = tmp_path / 'proteins.fa'
    fasta.write_text(PROTEIN_FASTA)
    output = tmp_path / 'collapsed.fa'

    stats = collapse(ft, fasta, output)

    text = output.read_text()
    assert '>NP_002.1' in text  # longer isoform (product_length 160) of gene 171590
    assert '>NP_001.1' not in text  # shorter isoform of the same gene, dropped
    assert '>NP_003.1' in text  # the only isoform of gene 171591
    assert '>NP_999.1' not in text  # no feature-table entry -- excluded, not fatal

    assert stats == {
        'input_records': 4,
        'genes': 2,
        'kept': 2,
        'dropped_isoforms': 1,
        'unmapped_no_gene': 1,
    }


def test_collapse_handles_gzipped_feature_table(tmp_path):
    import gzip
    ft = tmp_path / 'feature_table.txt.gz'
    with gzip.open(ft, 'wt') as fh:
        fh.write(FEATURE_TABLE)
    fasta = tmp_path / 'proteins.fa'
    fasta.write_text(PROTEIN_FASTA)
    output = tmp_path / 'collapsed.fa'

    stats = collapse(ft, fasta, output)
    assert stats['genes'] == 2
