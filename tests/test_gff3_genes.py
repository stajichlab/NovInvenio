import gzip

from gff3_genes import (
    gene_id_from_protein_id,
    load_gff3_index,
    lookup_gene_position,
    parse_gff3,
    resolve_gff3_path,
    resolve_gff3_paths,
)

# Ncra1 has a gene with a matching top-level ID (Afu1g01620-T-p1 -> Afu1g01620 style).
# Ncra2's gene ID never separately matches its protein ID, but its mRNA feature keeps
# the -T1 transcript suffix -- the case the user specifically flagged (NCU00499-T1-p1).
GFF3 = """\
##gff-version 3
scaffold_1\tfunannotate\tgene\t100\t500\t.\t+\t.\tID=Ncra1;Name=Ncra1
scaffold_1\tfunannotate\tmRNA\t100\t500\t.\t+\t.\tID=Ncra1-T1;Parent=Ncra1
scaffold_2\tfunannotate\tgene\t9000\t9800\t.\t-\t.\tID=gene_NCU00499_locus
scaffold_2\tfunannotate\tmRNA\t9000\t9800\t.\t-\t.\tID=NCU00499-T1;Parent=gene_NCU00499_locus
scaffold_3\tfunannotate\tgene\t42\t142\t.\t+\t.\tID=Direct1
scaffold_3\tfunannotate\tCDS\t42\t142\t.\t+\t0\tID=Direct1-cds;Parent=Direct1
"""


def test_parse_gff3_indexes_gene_and_mrna_ids(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    index = parse_gff3(p)
    assert index['Ncra1'] == ('scaffold_1', 100)
    assert index['Ncra1-T1'] == ('scaffold_1', 100)
    assert index['NCU00499-T1'] == ('scaffold_2', 9000)
    assert index['gene_NCU00499_locus'] == ('scaffold_2', 9000)
    # CDS features are not indexed at all.
    assert 'Direct1-cds' not in index


def test_parse_gff3_ignores_comments_and_short_lines(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text('# comment\n\nscaffold_1\tsrc\tgene\t1\t2\n' + GFF3)
    index = parse_gff3(p)
    assert index['Ncra1'] == ('scaffold_1', 100)


def test_parse_gff3_gene_wins_on_id_collision(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(
        'scaffold_1\tsrc\tmRNA\t1\t2\t.\t+\t.\tID=Dup\n'
        'scaffold_1\tsrc\tgene\t500\t600\t.\t+\t.\tID=Dup\n'
    )
    index = parse_gff3(p)
    assert index['Dup'] == ('scaffold_1', 500)


def test_parse_gff3_handles_gzip(tmp_path):
    p = tmp_path / 'test.gff3.gz'
    with gzip.open(p, 'wt') as fh:
        fh.write(GFF3)
    index = parse_gff3(p)
    assert index['Ncra1'] == ('scaffold_1', 100)


def test_gene_id_strips_transcript_suffixes():
    assert gene_id_from_protein_id('Afu1g01620-T-p1') == 'Afu1g01620'
    assert gene_id_from_protein_id('NCU00499-t26_1-p1') == 'NCU00499'
    assert gene_id_from_protein_id('AOBS_000123') == 'AOBS_000123'


def test_lookup_gene_position_direct_id_match(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    index = parse_gff3(p)
    assert lookup_gene_position('Direct1', index) == ('scaffold_3', 42)


def test_lookup_gene_position_falls_back_to_gene_id_strip(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    index = parse_gff3(p)
    # Ncra1-T-p1 strips (via gene_id_from_protein_id) to Ncra1, which matches the gene ID.
    assert lookup_gene_position('Ncra1-T-p1', index) == ('scaffold_1', 100)


def test_lookup_gene_position_falls_back_to_mrna_id_when_gene_id_does_not_match(tmp_path):
    # The case the user flagged: NCU00499-T1-p1's gene-id-strip heuristic collapses to
    # "NCU00499", which the GFF3 never gives the gene as a top-level ID -- but stripping
    # only the protein suffix recovers "NCU00499-T1", which matches the mRNA feature.
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    index = parse_gff3(p)
    assert lookup_gene_position('NCU00499-T1-p1', index) == ('scaffold_2', 9000)


def test_lookup_gene_position_returns_none_when_nothing_matches(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    index = parse_gff3(p)
    assert lookup_gene_position('totally_unrelated_id-p1', index) is None


def test_load_gff3_index_caches_by_path(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    cache = {}
    first = load_gff3_index(p, cache)
    p.write_text('scaffold_9\tsrc\tgene\t1\t2\t.\t+\t.\tID=Changed\n')
    second = load_gff3_index(p, cache)
    assert first is second
    assert 'Changed' not in second  # proves the second call was served from cache


def test_load_gff3_index_without_cache_reparses_every_time(tmp_path):
    p = tmp_path / 'test.gff3'
    p.write_text(GFF3)
    assert load_gff3_index(p) == parse_gff3(p)


def test_resolve_gff3_path_flat_layout(tmp_path):
    (tmp_path / 'Ncra.gff3').write_text(GFF3)
    resolved = resolve_gff3_path(tmp_path, 'Ncra.gff3')
    assert resolved == tmp_path / 'Ncra.gff3'


def test_resolve_gff3_path_searches_gff3_subdir(tmp_path):
    (tmp_path / 'gff3').mkdir()
    (tmp_path / 'gff3' / 'Ncra.gff3').write_text(GFF3)
    resolved = resolve_gff3_path(tmp_path, 'Ncra.gff3')
    assert resolved == tmp_path / 'gff3' / 'Ncra.gff3'


def test_resolve_gff3_path_returns_none_when_missing():
    assert resolve_gff3_path('/nonexistent/dir', 'Ncra.gff3') is None


def test_resolve_gff3_path_returns_none_for_empty_inputs(tmp_path):
    assert resolve_gff3_path(tmp_path, '') is None
    assert resolve_gff3_path('', 'Ncra.gff3') is None


class _Sample:
    def __init__(self, short, gff3=''):
        self.short = short
        self.gff3 = gff3


def test_resolve_gff3_paths_skips_samples_with_no_gff3_or_unresolvable_path(tmp_path):
    (tmp_path / 'Ncra.gff3').write_text(GFF3)
    samples = [
        _Sample('Ncra', 'Ncra.gff3'),
        _Sample('Afum', ''),  # no GFF3 column value
        _Sample('Spom', 'Spom.gff3'),  # doesn't exist under tmp_path
    ]
    paths = resolve_gff3_paths(samples, tmp_path)
    assert set(paths) == {'Ncra'}
    assert paths['Ncra'] == str(tmp_path / 'Ncra.gff3')


def test_resolve_gff3_paths_empty_without_data_dir():
    samples = [_Sample('Ncra', 'Ncra.gff3')]
    assert resolve_gff3_paths(samples, None) == {}
