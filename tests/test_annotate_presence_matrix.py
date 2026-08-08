import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from annotate_presence_matrix import parse_pfam_tblout  # noqa: E402

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
