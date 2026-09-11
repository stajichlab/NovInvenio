import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from annotate_presence_matrix import UNIPROT_XREF_COLS, load_uniprot_xrefs, parse_pfam_tblout  # noqa: E402

# A real hmmsearch --tblout line (target = candidate protein, query = Pfam-A HMM). Column
# positions that matter to the parser:
#   parts[0] protein/target id, parts[2] domain name, parts[3] accession, parts[4] E-value
TBLOUT = """\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   query name           accession   E-value  score  bias   E-value  score  bias  hmm from   hmm to  ali from   ali to   env from   env to  acc description of target
candidate1           -           Protein_kinase        PF00069.28  1.2e-30  105.3   0.1   3.4e-31  104.9   0.1     1      260      10      265       8      270  0.95 protein kinase domain
candidate2           -           Protein_kinase        PF00069.28  2.0e-25   90.1   0.2   5.1e-26   89.7   0.2     1      258      12      260      10      263  0.93 protein kinase domain
candidate1           -           zf-C2H2               PF00096.28  4.5e-08   30.2   1.1   6.0e-08   29.8   1.1     1       23       5       27       3       28  0.88 zinc finger
"""


def test_parse_pfam_tblout_uses_hmmsearch_column_order(tmp_path):
    tblout = tmp_path / 'candidates.pfam.tblout'
    tblout.write_text(TBLOUT)

    hits = parse_pfam_tblout(str(tblout))

    assert set(hits) == {'candidate1', 'candidate2'}
    # candidate1 has two distinct domains; order preserved, no duplicates.
    names = [name for name, acc, evalue in hits['candidate1']]
    assert names == ['Protein_kinase', 'zf-C2H2']
    accs = [acc for name, acc, evalue in hits['candidate1']]
    assert accs == ['PF00069.28', 'PF00096.28']

    assert hits['candidate2'] == [('Protein_kinase', 'PF00069.28', '2.0e-25')]


def test_parse_pfam_tblout_dedupes_repeated_domain(tmp_path):
    tblout = tmp_path / 'candidates.pfam.tblout'
    # Same protein, same domain hit twice (e.g. two envelopes) -- only the first is kept.
    tblout.write_text(
        "candidate1 - Protein_kinase PF00069.28 1.2e-30 105.3 0.1 3.4e-31 104.9 0.1 "
        "1 1 1 260 10 265 8 270 0.95 first\n"
        "candidate1 - Protein_kinase PF00069.28 9.9e-05  20.0 0.1 9.9e-05  20.0 0.1 "
        "1 1 1 260 280 300 275 305 0.80 second\n"
    )

    hits = parse_pfam_tblout(str(tblout))

    assert hits['candidate1'] == [('Protein_kinase', 'PF00069.28', '1.2e-30')]


def _write_xref_tsv(path, rows):
    with open(path, 'w', newline='') as fh:
        fieldnames = ['protein_id'] + UNIPROT_XREF_COLS
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow(row)


def test_load_uniprot_xrefs_merges_multiple_species_files_by_protein_id(tmp_path):
    f1 = tmp_path / 'Ccin.tsv'
    f2 = tmp_path / 'Agbi.tsv'
    _write_xref_tsv(f1, [{'protein_id': 'XP_001.1', 'uniprot_accession': 'A8MZR5',
                          'uniprot_xrefs': 'RefSeq:XP_001.1|GeneID:123'}])
    _write_xref_tsv(f2, [{'protein_id': 'XP_002.1', 'uniprot_accession': 'K5VHC4',
                          'uniprot_xrefs': 'RefSeq:XP_002.1|GeneID:456'}])

    merged = load_uniprot_xrefs([str(f1), str(f2)])

    assert set(merged) == {'XP_001.1', 'XP_002.1'}
    assert merged['XP_001.1']['uniprot_accession'] == 'A8MZR5'
    assert merged['XP_002.1']['uniprot_xrefs'] == 'RefSeq:XP_002.1|GeneID:456'


def test_annotate_presence_matrix_cli_adds_uniprot_columns(tmp_path):
    matrix = tmp_path / 'presence_matrix.tsv'
    matrix.write_text(
        "protein_id\tsource_proteome\tSpeciesA\n"
        "XP_001.1\tSpeciesA\t1\n"
        "XP_999.1\tSpeciesA\t1\n"
    )
    xref_tsv = tmp_path / 'SpeciesA.uniprot_xref.tsv'
    _write_xref_tsv(xref_tsv, [{'protein_id': 'XP_001.1', 'uniprot_accession': 'A8MZR5',
                                'uniprot_gene_name': 'CC1G_13429',
                                'uniprot_xrefs': 'RefSeq:XP_001.1|GeneID:123'}])
    output = tmp_path / 'presence_matrix.function.tsv'

    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'annotate_presence_matrix.py'),
         '--matrix', str(matrix), '--uniprot_xref_files', str(xref_tsv),
         '--output', str(output)],
        check=True,
    )

    rows = {r['protein_id']: r for r in csv.DictReader(open(output), delimiter='\t')}
    assert rows['XP_001.1']['uniprot_accession'] == 'A8MZR5'
    assert rows['XP_001.1']['uniprot_xrefs'] == 'RefSeq:XP_001.1|GeneID:123'
    # A protein with no crosswalk match gets empty uniprot_* columns, not an error.
    assert rows['XP_999.1']['uniprot_accession'] == ''
    assert rows['XP_999.1']['uniprot_xrefs'] == ''


def test_annotate_presence_matrix_cli_without_uniprot_flag_omits_columns(tmp_path):
    matrix = tmp_path / 'presence_matrix.tsv'
    matrix.write_text("protein_id\tsource_proteome\tSpeciesA\nXP_001.1\tSpeciesA\t1\n")
    output = tmp_path / 'presence_matrix.function.tsv'

    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'annotate_presence_matrix.py'),
         '--matrix', str(matrix), '--output', str(output)],
        check=True,
    )

    header = open(output).readline().rstrip('\n').split('\t')
    assert 'uniprot_xrefs' not in header
