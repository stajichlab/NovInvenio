"""Issue #208: runtime behaviour of the Outgroup signal controls on novelties.html (jsdom)."""
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DRIVER = REPO / 'tests' / 'js' / 'drive_other_signal.mjs'
sys.path.insert(0, str(REPO / 'tests'))
from test_report_data import CONFIG, FASTA, MATRIX, TBLASTN  # noqa: E402


def _jsdom():
    env = os.environ.get('NOVINVENIO_JSDOM')
    cands = [env] if env else []
    cands += [str(Path.home() / '.cache/novinvenio-jsdom/node_modules/jsdom')]
    for c in cands:
        if c and Path(c).exists():
            return c
    return None


@pytest.mark.skipif(shutil.which('node') is None or _jsdom() is None, reason='node + jsdom not available')
def test_signal_controls_filter_and_reclassify(tmp_path):
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'matrix.tsv').write_text(MATRIX)
    (tmp_path / 'tblastn.tsv').write_text(TBLASTN)
    (tmp_path / 'candidates.fa').write_text(FASTA)
    with gzip.open(tmp_path / 'ev.tsv.gz', 'wt') as fh:
        fh.write('protein_id\tsource_proteome\tother_proteome\tevidence\tbest_target_id\tbest_evalue\t'
                 'best_bitscore\tbest_qcov\tmax_qcov\tn_hits\tstatus\tparalog_id\tparalog_evalue\tparalog_delta\n')
        fh.write('n1\tNcra\tSpom\tprotein\tt1\t1e-10\t60\t20\t20\t1\tparalog_filtered\tp1\t1e-100\t90\n')
    with gzip.open(tmp_path / 'cov.tsv.gz', 'wt') as fh:
        fh.write('protein_id\trep_id\tgenome\tn_hsps\tbest_evalue\tquery_span_cov\tspan_start\tspan_end\n')
        fh.write('n2\tn2\tScer\t1\t1e-8\t0.800000\t5\t85\n')
    subprocess.run([sys.executable, str(REPO / 'bin' / 'make_report.py'),
                    '--config', str(tmp_path / 'config.csv'), '--matrix', str(tmp_path / 'matrix.tsv'),
                    '--tblastn_summary', str(tmp_path / 'tblastn.tsv'),
                    '--candidates_fa', str(tmp_path / 'candidates.fa'),
                    '--other_evidence', str(tmp_path / 'ev.tsv.gz'),
                    '--tblastn_coverage', str(tmp_path / 'cov.tsv.gz'),
                    '--output', str(tmp_path / 'novelties.html')],
                   check=True, capture_output=True, text=True)
    r = subprocess.run(['node', str(DRIVER), str(tmp_path / 'novelties.html'), _jsdom()],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr


@pytest.mark.skipif(shutil.which('node') is None or _jsdom() is None, reason='node + jsdom not available')
def test_losses_signal_controls_filter_and_reclassify(tmp_path):
    from test_report_data import LOSSES_MATRIX
    (tmp_path / 'config.csv').write_text(CONFIG)
    (tmp_path / 'm.tsv').write_text(LOSSES_MATRIX)
    with gzip.open(tmp_path / 'ev.tsv.gz', 'wt') as fh:
        fh.write('protein_id\tsource_proteome\tother_proteome\tevidence\tbest_target_id\tbest_evalue\t'
                 'best_bitscore\tbest_qcov\tmax_qcov\tn_hits\tstatus\tparalog_id\tparalog_evalue\tparalog_delta\n')
        fh.write('loss1\tSpom\tNcra\tprotein\tt1\t1e-10\t60\t30\t30\t1\tparalog_filtered\tp1\t1e-90\t80\n')
    with gzip.open(tmp_path / 'cov.tsv.gz', 'wt') as fh:
        fh.write('protein_id\trep_id\tgenome\tn_hsps\tbest_evalue\tquery_span_cov\tspan_start\tspan_end\n')
        fh.write('loss1b\tloss1b\tAfum\t1\t1e-8\t0.700000\t1\t70\n')
    subprocess.run([sys.executable, str(REPO / 'bin' / 'make_losses_report.py'),
                    '--config', str(tmp_path / 'config.csv'), '--matrix', str(tmp_path / 'm.tsv'),
                    '--other_evidence', str(tmp_path / 'ev.tsv.gz'),
                    '--tblastn_coverage', str(tmp_path / 'cov.tsv.gz'),
                    '--output', str(tmp_path / 'losses.html')],
                   check=True, capture_output=True, text=True)
    r = subprocess.run(['node', str(DRIVER), str(tmp_path / 'losses.html'), _jsdom()],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
