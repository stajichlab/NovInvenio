"""Unit tests for bin/calibrate_family_hmms.py."""
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / 'bin' / 'calibrate_family_hmms.py'


def _domtblout_line(target, query, full_e=1e-10, qlen=100, hmm_from=1, hmm_to=60):
    """One hmmsearch --domtblout data row (23 whitespace fields)."""
    f = ['-'] * 23
    f[0] = target
    f[3] = query
    f[5] = str(qlen)
    f[6] = f'{full_e:g}'
    f[15] = str(hmm_from)
    f[16] = str(hmm_to)
    return ' '.join(f)


def _write_domtblout(path, hits):
    with open(path, 'w') as fh:
        fh.write("# hmmsearch domtblout\n")
        for target, query, *rest in hits:
            evalue = rest[0] if rest else 1e-10
            fh.write(_domtblout_line(target, query, full_e=evalue) + "\n")


def test_negative_control_calibration(tmp_path):
    """Family with a DISCOVERY_OUT hit → threshold = best outgroup E-value."""
    families = tmp_path / 'families.tsv'
    families.write_text("family_index\trepresentative_id\tn_members\n"
                        "fam_000001\tfamA\t2\n"
                        "fam_000002\tfamB\t2\n")

    dom = tmp_path / 'disc_out.domtblout'
    _write_domtblout(dom, [('target1', 'famA', 1e-5)])

    out = tmp_path / 'thresholds.tsv'
    subprocess.run(
        [sys.executable, str(BIN),
         '--domtblout', str(dom),
         '--families', str(families),
         '--default-evalue', '1e-3',
         '--output', str(out)],
        check=True,
    )

    lines = out.read_text().strip().split('\n')
    # famA has a hit at 1e-5 → threshold = 1e-5, source = negative_control
    famA_line = [ln for ln in lines if ln.startswith('famA\t')][0]
    parts = famA_line.split('\t')
    assert parts[2] == 'negative_control'
    assert float(parts[1]) == 1e-5


def test_global_default_for_no_hit(tmp_path):
    """Family with no DISCOVERY_OUT hit → threshold = global default."""
    families = tmp_path / 'families.tsv'
    families.write_text("family_index\trepresentative_id\tn_members\n"
                        "fam_000001\tfamA\t2\n")

    dom = tmp_path / 'disc_out.domtblout'
    _write_domtblout(dom, [])  # no hits

    out = tmp_path / 'thresholds.tsv'
    subprocess.run(
        [sys.executable, str(BIN),
         '--domtblout', str(dom),
         '--families', str(families),
         '--default-evalue', '1e-3',
         '--output', str(out)],
        check=True,
    )

    lines = out.read_text().strip().split('\n')
    famA_line = [ln for ln in lines if ln.startswith('famA\t')][0]
    parts = famA_line.split('\t')
    assert parts[2] == 'global_default'
    assert float(parts[1]) == 1e-3


def test_best_evalue_across_multiple_domtblouts(tmp_path):
    """Multiple DISCOVERY_OUT proteomes → threshold = lowest (best) E-value across all."""
    families = tmp_path / 'families.tsv'
    families.write_text("family_index\trepresentative_id\tn_members\n"
                        "fam_000001\tfamA\t2\n")

    dom1 = tmp_path / 'disc1.domtblout'
    dom2 = tmp_path / 'disc2.domtblout'
    _write_domtblout(dom1, [('t1', 'famA', 1e-5)])
    _write_domtblout(dom2, [('t2', 'famA', 1e-8)])  # better hit

    out = tmp_path / 'thresholds.tsv'
    subprocess.run(
        [sys.executable, str(BIN),
         '--domtblout', str(dom1), str(dom2),
         '--families', str(families),
         '--default-evalue', '1e-3',
         '--output', str(out)],
        check=True,
    )

    lines = out.read_text().strip().split('\n')
    famA_line = [ln for ln in lines if ln.startswith('famA\t')][0]
    parts = famA_line.split('\t')
    assert float(parts[1]) == 1e-8
