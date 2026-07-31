"""Unit tests for bin/extract_family_seqs.py, including the --max-members oversized-family
filter (issue #22 follow-up): a family this large is an ancient multi-copy superfamily that
can't pass the novelty/loss predicate and dominates famsa's runtime, so it's skipped and
reported rather than profiled.
"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'extract_family_seqs.py'

CLUSTER = (
    "small\tsmall\n"                              # 1 member -> below min (2), skipped
    "mid\tmid\nmid\tmid2\n"                        # 2 members -> in range
    "big\tbig\nbig\tb2\nbig\tb3\nbig\tb4\n"        # 4 members -> above max (3), skipped
)
FASTA = ">small\nMA\n>mid\nMB\n>mid2\nMC\n>big\nMD\n>b2\nME\n>b3\nMF\n>b4\nMG\n"


def _setup(tmp_path):
    (tmp_path / 'cluster.tsv').write_text(CLUSTER)
    (tmp_path / 'seed.faa').write_text(FASTA)


def _run(tmp_path, *extra):
    outdir = tmp_path / 'out'
    outdir.mkdir()
    cmd = [
        sys.executable, str(BIN),
        '--cluster-tsv', str(tmp_path / 'cluster.tsv'),
        '--fasta', str(tmp_path / 'seed.faa'),
        '--min-members', '2',
        '--outdir', str(outdir),
        *extra,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r, outdir


def test_max_members_skips_oversized_family_and_reports_it(tmp_path):
    _setup(tmp_path)
    r, outdir = _run(tmp_path, '--max-members', '3',
                     '--oversized-report', str(tmp_path / 'oversized.tsv'))
    assert r.returncode == 0, r.stderr

    families = (outdir / 'families.tsv').read_text().splitlines()
    reps = [line.split('\t')[1] for line in families[1:]]
    assert reps == ['mid']   # small: below min; big: above max — both excluded

    oversized = (tmp_path / 'oversized.tsv').read_text().splitlines()
    assert oversized[1].split('\t') == ['big', '4']
    assert 'oversized' in r.stderr


def test_no_max_members_keeps_large_families(tmp_path):
    _setup(tmp_path)
    r, outdir = _run(tmp_path)  # no --max-members => unlimited, matches prior behaviour
    assert r.returncode == 0, r.stderr
    families = (outdir / 'families.tsv').read_text().splitlines()
    reps = sorted(line.split('\t')[1] for line in families[1:])
    assert reps == ['big', 'mid']


def test_oversized_report_omitted_when_no_family_exceeds_max(tmp_path):
    _setup(tmp_path)
    report = tmp_path / 'oversized.tsv'
    r, outdir = _run(tmp_path, '--max-members', '10', '--oversized-report', str(report))
    assert r.returncode == 0, r.stderr
    assert report.exists()
    assert report.read_text().splitlines() == ['representative_id\tn_members']


def test_fam_fasta_written_only_for_in_range_families(tmp_path):
    _setup(tmp_path)
    r, outdir = _run(tmp_path, '--max-members', '3')
    assert r.returncode == 0, r.stderr
    fam_files = sorted(outdir.glob('fam_*.faa'))
    assert len(fam_files) == 1
    assert fam_files[0].read_text() == '>mid\nMB\n>mid2\nMC\n'
