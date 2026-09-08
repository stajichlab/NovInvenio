import gzip
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))

from build_alignment_shards import (  # noqa: E402
    build_genome_shard, genome_short, load_candidate_ids, main, parse_cluster_tsv,
)

# Column order matching modules/tblastn.nf's -outfmt 6 string (issue #70):
# qseqid sseqid evalue bitscore pident length qstart qend sstart send sframe qseq sseq
_HEADER = ('qseqid', 'sseqid', 'evalue', 'bitscore', 'pident', 'length',
           'qstart', 'qend', 'sstart', 'send', 'sframe', 'qseq', 'sseq')


def _row(**kw):
    defaults = dict(
        qseqid='rep1', sseqid='scaffold_1', evalue='1e-40', bitscore='150',
        pident='75.0', length='40', qstart='1', qend='40', sstart='1000',
        send='1120', sframe='1', qseq='MKVL-ACDEF', sseq='MKVLQACDEF',
    )
    defaults.update(kw)
    return '\t'.join(str(defaults[c]) for c in _HEADER)


def test_genome_short_strips_known_suffixes():
    assert genome_short('/x/Afum.tblastn.tsv') == 'Afum'
    assert genome_short('Afum.tblastn.tsv.gz') == 'Afum'
    assert genome_short('Weird.name.tblastn.tsv') == 'Weird.name'


def test_load_candidate_ids_strips_short_prefix(tmp_path):
    p = tmp_path / 'candidates.txt'
    p.write_text('Ncra::protA\nAfum::protB\n')
    assert load_candidate_ids(p) == {'protA', 'protB'}


def test_parse_cluster_tsv_maps_member_to_rep(tmp_path):
    p = tmp_path / 'clusters.tsv'
    p.write_text('rep1\trep1\nrep1\tmemberA\nrep2\trep2\n')
    assert parse_cluster_tsv(p) == {'rep1': 'rep1', 'memberA': 'rep1', 'rep2': 'rep2'}


def test_build_genome_shard_expands_representative_to_candidate_members(tmp_path):
    tsv = tmp_path / 'Afum.tblastn.tsv'
    tsv.write_text(_row() + '\n')
    member_to_rep = {'rep1': 'rep1', 'memberA': 'rep1'}

    shard = build_genome_shard(tsv, 'Afum', 1e-5, member_to_rep, {'memberA'})

    assert set(shard.keys()) == {'memberA'}
    hit = shard['memberA'][0]
    assert hit['genome'] == 'Afum'
    assert hit['aligned_as'] == 'rep1'
    assert hit['evalue'] == 1e-40
    assert hit['qseq'] == 'MKVL-ACDEF' and hit['sseq'] == 'MKVLQACDEF'


def test_build_genome_shard_no_aligned_as_when_member_is_the_representative(tmp_path):
    tsv = tmp_path / 'Afum.tblastn.tsv'
    tsv.write_text(_row() + '\n')
    member_to_rep = {'rep1': 'rep1'}

    shard = build_genome_shard(tsv, 'Afum', 1e-5, member_to_rep, {'rep1'})

    assert 'aligned_as' not in shard['rep1'][0]


def test_build_genome_shard_drops_hits_above_evalue_cutoff(tmp_path):
    tsv = tmp_path / 'Afum.tblastn.tsv'
    tsv.write_text(_row(evalue='1.0') + '\n')
    member_to_rep = {'rep1': 'rep1'}

    shard = build_genome_shard(tsv, 'Afum', 1e-5, member_to_rep, {'rep1'})

    assert shard == {}


def test_build_genome_shard_drops_non_candidate_members(tmp_path):
    tsv = tmp_path / 'Afum.tblastn.tsv'
    tsv.write_text(_row() + '\n')
    member_to_rep = {'rep1': 'rep1', 'memberA': 'rep1'}

    shard = build_genome_shard(tsv, 'Afum', 1e-5, member_to_rep, {'someOtherCandidate'})

    assert shard == {}


def test_build_genome_shard_collects_multiple_hsps_as_a_list(tmp_path):
    tsv = tmp_path / 'Afum.tblastn.tsv'
    tsv.write_text(_row(sseqid='scaffold_1') + '\n' + _row(sseqid='scaffold_2', sstart='5000', send='5120') + '\n')
    member_to_rep = {'rep1': 'rep1'}

    shard = build_genome_shard(tsv, 'Afum', 1e-5, member_to_rep, {'rep1'})

    assert len(shard['rep1']) == 2
    assert {h['sseqid'] for h in shard['rep1']} == {'scaffold_1', 'scaffold_2'}


def test_end_to_end_writes_shard_and_manifest(tmp_path, monkeypatch):
    hits = tmp_path / 'Afum.tblastn.tsv'
    hits.write_text(_row() + '\n')
    candidates = tmp_path / 'candidates.txt'
    candidates.write_text('Ncra::memberA\n')
    cluster_tsv = tmp_path / 'clusters.tsv'
    cluster_tsv.write_text('rep1\trep1\nrep1\tmemberA\n')
    outdir = tmp_path / 'alignments'

    argv = [
        'build_alignment_shards.py',
        '--hits', str(hits),
        '--candidates', str(candidates),
        '--cluster_tsv', str(cluster_tsv),
        '--project', 'demo',
        '--outdir', str(outdir),
    ]
    monkeypatch.setattr(sys, 'argv', argv)
    main()

    manifest = json.loads((outdir / 'manifest.json').read_text())
    assert manifest['project'] == 'demo'
    assert manifest['schema_version'] == 1
    assert manifest['shards'] == ['Afum.json.gz']

    with gzip.open(outdir / 'Afum.json.gz', 'rt') as fh:
        shard = json.load(fh)
    assert shard['memberA'][0]['aligned_as'] == 'rep1'


def test_end_to_end_skips_genomes_with_no_candidate_hits(tmp_path, monkeypatch):
    hits = tmp_path / 'Afum.tblastn.tsv'
    hits.write_text(_row() + '\n')
    candidates = tmp_path / 'candidates.txt'
    candidates.write_text('Ncra::nobodyHome\n')
    cluster_tsv = tmp_path / 'clusters.tsv'
    cluster_tsv.write_text('rep1\trep1\n')
    outdir = tmp_path / 'alignments'

    argv = [
        'build_alignment_shards.py',
        '--hits', str(hits),
        '--candidates', str(candidates),
        '--cluster_tsv', str(cluster_tsv),
        '--project', 'demo',
        '--outdir', str(outdir),
    ]
    monkeypatch.setattr(sys, 'argv', argv)
    main()

    manifest = json.loads((outdir / 'manifest.json').read_text())
    assert manifest['shards'] == []
    assert not (outdir / 'Afum.json.gz').exists()
