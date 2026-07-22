"""Unit tests for lib/figures.py + bin/make_pdf_report.py (PDF summary report, issue #19).

Builds the figure set from a fixture payload (reusing report_data.build_payload, the same
path the pipeline uses) and asserts each figure renders; plus an end-to-end run of
make_pdf_report.py writing a real multi-page PDF. No display needed (Agg backend).
"""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))

from config_parser import parse_config  # noqa: E402
from report_data import build_payload, build_losses_payload  # noqa: E402
import figures  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
OUT,Saccharomyces cerevisiae,S288c,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina
"""

MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions\tPfam_Evalues
n1\tNcra\t1\t1\t0\t0\tada-1\tthing\tModelOrg_Ncra\t\tbZIP_1\tPF00170.27\t4e-9
n2\tAfum\t1\t1\t0\t0\t\t\tPfam\t\tp450\tPF00067.1\t1e-8
shared\tNcra\t1\t1\t1\t0\t\tconserved\tPfam\t\tAAA\tPF00004.31\t1e-20
lonely\tNcra\t1\t0\t0\t0\t\t\t\t\t\t\t
"""

TBLASTN = """\
protein_id\tSpom\tScer
n2\t0\t1
"""

CLUSTER = "n1\tn1\nn1\tn2\n"

LOSSES_MATRIX = """\
protein_id\tsource_proteome\tNcra\tAfum\tSpom\tScer\tgene_name\tproduct_description\tfunction_source\tBest_Swissprot\tPfam_Names\tPfam_Accessions
loss1\tSpom\t0\t0\t1\t1\tERG3\tsterol\tPfam\t\tp450\tPF00067.1
loss1b\tScer\t0\t0\t1\t1\t\t\t\t\t\t
"""


@pytest.fixture
def payload(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'tblastn.tsv').write_text(TBLASTN)
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    samples = parse_config(tmp_path / 'config.csv')
    return build_payload(tmp_path / 'matrix.tsv', samples,
                         tblastn_path=tmp_path / 'tblastn.tsv',
                         cluster_tsv=tmp_path / 'cluster.tsv', sequences='none')


@pytest.fixture
def losses_payload(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'loss_matrix.tsv').write_text(LOSSES_MATRIX)
    samples = parse_config(tmp_path / 'config.csv')
    return build_losses_payload(tmp_path / 'loss_matrix.tsv', samples)


def test_each_figure_renders(payload):
    for fn in (figures.fig_summary, figures.fig_novelty_per_species,
               figures.fig_presence_heatmap, figures.fig_family_sizes,
               figures.fig_pfam_frequency, figures.fig_annotation_source):
        fig = fn(payload)
        assert isinstance(fig, Figure)
        assert len(fig.axes) >= 0  # rendered without raising


def test_loss_scatter_renders(losses_payload):
    fig = figures.fig_loss_scatter(losses_payload)
    assert isinstance(fig, Figure)


def test_build_all_includes_loss_page_only_when_given(payload, losses_payload):
    assert len(figures.build_all(payload)) == 6
    assert len(figures.build_all(payload, losses_payload)) == 7


def test_palette_is_the_validated_set():
    # guard against silent palette drift (the CVD-validated categorical set)
    assert figures.CATEGORICAL == ['#2a78d6', '#e69f00', '#cc79a7', '#008300']


def test_end_to_end_writes_pdf(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'tblastn.tsv').write_text(TBLASTN)
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    (tmp_path / 'loss_matrix.tsv').write_text(LOSSES_MATRIX)
    out = tmp_path / 'summary.pdf'
    subprocess.run([
        sys.executable, str(REPO / 'bin' / 'make_pdf_report.py'),
        '--matrix', str(tmp_path / 'matrix.tsv'),
        '--config', str(tmp_path / 'config.csv'),
        '--tblastn_summary', str(tmp_path / 'tblastn.tsv'),
        '--cluster_tsv', str(tmp_path / 'cluster.tsv'),
        '--loss_matrix', str(tmp_path / 'loss_matrix.tsv'),
        '--output', str(out),
    ], check=True, capture_output=True, text=True)
    assert out.exists()
    assert out.stat().st_size > 1000
    assert out.read_bytes()[:5] == b'%PDF-'  # a real PDF
