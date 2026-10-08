"""Issue #208: other-group signal columns in novelties.<SHORT>.tsv and the class boundaries."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))
from other_evidence import candidate_signal, signal_class  # noqa: E402

CONFIG = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,P
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,P
OUT,Schizosaccharomyces pombe,972h,Spom.pep.fa,Spom.dna.fa,Spom,T
"""
MATRIX = "protein_id\tsource_proteome\tNcra\tAfum\tSpom\nn1\tNcra\t1\t1\t0\nn2\tNcra\t1\t1\t0\nn3\tNcra\t1\t1\t0\n"


def test_class_boundaries():
    t = 50
    assert signal_class(49.99, t, True) == 'domain_only'
    assert signal_class(50, t, True) == 'broad'
    assert signal_class(None, t, False) == 'none'
    assert signal_class(None, t, True) == ''


def test_candidate_signal_uses_max_of_protein_and_scaled_tblastn():
    prot = [{'max_qcov': '20', 'status': 'paralog_filtered', 'paralog_evalue': '1e-100'},
            {'max_qcov': '30', 'status': 'kept', 'paralog_evalue': ''}]
    tb = [{'genome': 'G1', 'query_span_cov': '0.62'}, {'genome': 'G2', 'query_span_cov': '0.10'}]
    d = candidate_signal(prot, tb, 50)
    assert d['other_protein_max_qcov'] == '30'
    assert d['other_protein_n_filtered'] == '1'
    assert float(d['other_protein_min_paralog_evalue']) == 1e-100
    assert d['other_tblastn_max_cov'] == '0.620'
    assert d['other_tblastn_n_genomes'] == '2'
    assert d['other_signal_class'] == 'broad'        # 62 >= 50


def test_not_measured_stays_blank():
    d = candidate_signal(None, None, 50)
    assert all(v == '' for v in d.values())


def test_no_rows_is_none_when_a_sidecar_was_supplied():
    assert candidate_signal([], None, 50)['other_signal_class'] == 'none'


def _run(tmp_path, extra):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    subprocess.run([sys.executable, str(REPO / 'bin' / 'make_novelties.py'),
                    '--matrix', str(tmp_path / 'matrix.tsv'), '--config', str(tmp_path / 'config.csv'),
                    '--ingroup_min', '0.5', '--output_dir', str(tmp_path), *extra],
                   check=True, capture_output=True, text=True)
    with open(tmp_path / 'novelties.Ncra.tsv', newline='') as fh:
        return {r['protein_id']: r for r in csv.DictReader(fh, delimiter='\t')}


def test_make_novelties_writes_signal_columns(tmp_path):
    ev = tmp_path / 'ev.tsv.gz'
    with gzip.open(ev, 'wt') as fh:
        fh.write('protein_id\tsource_proteome\tother_proteome\tevidence\tbest_target_id\tbest_evalue\t'
                 'best_bitscore\tbest_qcov\tmax_qcov\tn_hits\tstatus\tparalog_id\tparalog_evalue\tparalog_delta\n')
        fh.write('n1\tNcra\tSpom\tprotein\tt\t1e-10\t60\t20\t20\t1\tparalog_filtered\tp\t1e-100\t90\n')
    cov = tmp_path / 'cov.tsv.gz'
    with gzip.open(cov, 'wt') as fh:
        fh.write('protein_id\trep_id\tgenome\tn_hsps\tbest_evalue\tquery_span_cov\tspan_start\tspan_end\n')
        fh.write('n2\tn2\tG1\t1\t1e-8\t0.800000\t1\t80\n')
    rows = _run(tmp_path, ['--other_evidence', str(ev), '--tblastn_coverage', str(cov)])
    assert rows['n1']['other_signal_class'] == 'domain_only'
    assert rows['n1']['other_protein_n_filtered'] == '1'
    assert rows['n2']['other_signal_class'] == 'broad'
    assert rows['n3']['other_signal_class'] == 'none'


def test_threshold_option_changes_the_class(tmp_path):
    cov = tmp_path / 'cov.tsv.gz'
    with gzip.open(cov, 'wt') as fh:
        fh.write('protein_id\trep_id\tgenome\tn_hsps\tbest_evalue\tquery_span_cov\tspan_start\tspan_end\n')
        fh.write('n2\tn2\tG1\t1\t1e-8\t0.600000\t1\t60\n')
    assert _run(tmp_path, ['--tblastn_coverage', str(cov)])['n2']['other_signal_class'] == 'broad'
    assert _run(tmp_path, ['--tblastn_coverage', str(cov), '--other_signal_qcov', '70'])['n2']['other_signal_class'] == 'domain_only'


def test_columns_absent_without_sidecars(tmp_path):
    rows = _run(tmp_path, [])
    assert 'other_signal_class' not in rows['n1']
