"""Issue #208: --output-coverage for bin/summarize_tblastn.py."""
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / 'bin' / 'summarize_tblastn.py'


def tb(q, ev, qs, qe):
    # qseqid sseqid evalue bitscore pident length qstart qend sstart send sframe qseq sseq
    return f'{q}\tctg1\t{ev}\t50\t40.0\t30\t{qs}\t{qe}\t1\t90\t1\tAAA\tAAA\n'


def run(tmp_path, files, cluster, fasta, extra=()):
    hits = []
    for name, text in files.items():
        p = tmp_path / f'{name}.tblastn.tsv'
        p.write_text(text)
        hits.append(str(p))
    (tmp_path / 'cl.tsv').write_text(cluster)
    (tmp_path / 'reps.fa').write_text(fasta)
    cmd = [sys.executable, str(SCRIPT), '--hits', *hits, '--cluster_tsv', str(tmp_path / 'cl.tsv'),
           '--output', str(tmp_path / 'sum.tsv'), *extra]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def rows(tmp_path):
    with gzip.open(tmp_path / 'cov.tsv.gz', 'rt') as fh:
        lines = fh.read().splitlines()
    head = lines[0].split('\t')
    return [dict(zip(head, l.split('\t'))) for l in lines[1:]]


def go(tmp_path, files, cluster='repA\trepA\nrepA\tmemB\n', length=100, evalue='1e-5'):
    run(tmp_path, files, cluster, f'>repA\n{"M" * length}\n',
        ['--output-coverage', str(tmp_path / 'cov.tsv.gz'),
         '--query_fasta', str(tmp_path / 'reps.fa'), '--evalue', evalue])
    return rows(tmp_path)


def test_union_of_overlapping_hsps_and_separate_hsp(tmp_path):
    # 10-40 and 30-60 overlap -> 10..60 = 51 residues; 80-90 adds 11 -> 62 of 100.
    text = tb('repA', '1e-10', 10, 40) + tb('repA', '1e-8', 30, 60) + tb('repA', '1e-6', 80, 90)
    r = go(tmp_path, {'G1': text})
    row = [x for x in r if x['protein_id'] == 'repA'][0]
    assert row['genome'] == 'G1' and row['rep_id'] == 'repA'
    assert row['n_hsps'] == '3'
    assert abs(float(row['query_span_cov']) - 0.62) < 1e-9
    assert row['span_start'] == '10' and row['span_end'] == '90'
    assert float(row['best_evalue']) == 1e-10


def test_hits_above_the_cutoff_are_ignored(tmp_path):
    text = tb('repA', '1e-10', 10, 20) + tb('repA', '1e-2', 50, 100)
    r = go(tmp_path, {'G1': text})
    row = r[0]
    assert row['n_hsps'] == '1'
    assert abs(float(row['query_span_cov']) - 0.11) < 1e-9


def test_rep_values_are_copied_to_cluster_members(tmp_path):
    r = go(tmp_path, {'G1': tb('repA', '1e-10', 1, 50)})
    by = {x['protein_id']: x for x in r}
    assert set(by) == {'repA', 'memB'}
    assert by['memB']['rep_id'] == 'repA'
    assert by['memB']['query_span_cov'] == by['repA']['query_span_cov']


def test_one_row_per_genome_with_a_hit(tmp_path):
    r = go(tmp_path, {'G1': tb('repA', '1e-10', 1, 50), 'G2': ''})
    assert {x['genome'] for x in r} == {'G1'}


def test_presence_summary_is_unchanged_by_coverage_output(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir(); b.mkdir()
    files = {'G1': tb('repA', '1e-10', 1, 50)}
    go(a, files)
    run(b, files, 'repA\trepA\nrepA\tmemB\n', '>repA\n' + 'M' * 100 + '\n')
    assert (a / 'sum.tsv').read_bytes() == (b / 'sum.tsv').read_bytes()
