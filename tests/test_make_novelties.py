import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
"""

# n1 (Ncra) and n2 (Afum) are the same gene family, recovered independently
# in both ingroup species; mmseqs clusters them together downstream.
MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions\tPfam_Evalues
n1\tNcra\t1\t1\t0\tada-1\tall development altered-1\tModelOrg_Ncra\t\tbZIP_1\tPF00170.27\t4.5e-09
n2\tAfum\t1\t1\t0\t\t\t\t\t\t\t
"""

CLUSTER_TSV = "n1\tn1\nn1\tn2\n"


@pytest.fixture
def run_dir(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'clusters_cluster.tsv').write_text(CLUSTER_TSV)
    return tmp_path


def run_make_novelties(run_dir, extra_args=()):
    subprocess.run(
        [sys.executable, str(REPO / 'bin' / 'make_novelties.py'),
         '--matrix', str(run_dir / 'matrix.tsv'),
         '--config', str(run_dir / 'config.csv'),
         '--ingroup_min', '0.5',
         '--output_dir', str(run_dir),
         *extra_args],
        check=True, capture_output=True, text=True,
    )


def read_novelties(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def test_family_columns_absent_without_cluster_tsv(run_dir):
    run_make_novelties(run_dir)
    rows = read_novelties(run_dir / 'novelties.Ncra.tsv')
    assert rows[0]['family_id'] == ''
    assert rows[0]['family_size'] == '1'
    assert rows[0]['family_members'] == ''


def test_family_columns_link_the_same_gene_across_ingroup_species(run_dir):
    run_make_novelties(run_dir, extra_args=['--cluster_tsv', str(run_dir / 'clusters_cluster.tsv')])

    ncra_rows = read_novelties(run_dir / 'novelties.Ncra.tsv')
    afum_rows = read_novelties(run_dir / 'novelties.Afum.tsv')
    assert len(ncra_rows) == 1
    assert len(afum_rows) == 1

    n1, n2 = ncra_rows[0], afum_rows[0]
    assert n1['family_id'] == 'n1'
    assert n2['family_id'] == 'n1'
    assert n1['family_size'] == n2['family_size'] == '2'
    assert n1['family_members'] == 'n2'
    assert n2['family_members'] == 'n1'


# Regression test: bin/make_novelties.py's load_groups() previously read row['GROUP']
# directly (bypassing lib/config_parser.py's GROUP_ALIASES normalization), so a config
# still using the pre-rename TARGET/DISC_OUT labels
# (todo/rename-novelty-discovery-group-labels.md) silently produced zero ingroup/outgroup
# matches -- caught when run_sordario.sh's real config (configs/sordario.csv, untracked,
# still uses TARGET/DISC_OUT) failed with "0 total novelties across 0 ingroup species".
OLD_ALIAS_CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
TARGET,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
TARGET,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
DISC_OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
"""


def test_old_group_aliases_still_produce_novelties(tmp_path):
    (tmp_path / 'config.csv').write_text(OLD_ALIAS_CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    run_make_novelties(tmp_path)
    rows = read_novelties(tmp_path / 'novelties.Ncra.tsv')
    assert len(rows) == 1
    assert rows[0]['protein_id'] == 'n1'


# Regression test (issue #48): a --cluster_tool pairwise config may legitimately list
# NEAR_INGROUP/BROAD_OUTGROUP rows now that CONTEXT_SEARCH exists to search them
# separately -- bin/build_presence_matrix.py deliberately never gives them a matrix
# column (see its own comment), so load_groups() must skip samples with no matching
# column instead of hard-erroring "columns not found in matrix".
CONTEXT_ROLE_CONFIG = CONFIG + (
    'NEAR_INGROUP,Near species,,Near1.pep.fa,Near1.dna.fa,Near1,Pezizomycotina\n'
    'BROAD_OUTGROUP,Broad species,,Broad1.pep.fa,Broad1.dna.fa,Broad1,Basidiomycota\n'
)


def test_near_ingroup_broad_outgroup_rows_without_a_matrix_column_are_skipped(tmp_path):
    (tmp_path / 'config.csv').write_text(CONTEXT_ROLE_CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)  # no Near1/Broad1 columns
    run_make_novelties(tmp_path)
    rows = read_novelties(tmp_path / 'novelties.Ncra.tsv')
    assert len(rows) == 1
    assert rows[0]['protein_id'] == 'n1'
