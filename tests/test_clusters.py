from clusters import build_families, read_cluster_tsv

CLUSTER_TSV = """\
n1\tn1
n1\tn2
n1\tn3
lonely_rep\tlonely_rep
"""


def test_read_cluster_tsv_maps_member_to_rep(tmp_path):
    path = tmp_path / 'clusters_cluster.tsv'
    path.write_text(CLUSTER_TSV)
    assert read_cluster_tsv(path) == {
        'n1': 'n1', 'n2': 'n1', 'n3': 'n1', 'lonely_rep': 'lonely_rep',
    }


def test_build_families_drops_singleton_clusters(tmp_path):
    path = tmp_path / 'clusters_cluster.tsv'
    path.write_text(CLUSTER_TSV)
    families = build_families(read_cluster_tsv(path))
    assert families == {'n1': ['n1', 'n2', 'n3']}
