import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'refine_ambiguous_families.py'
sys.path.insert(0, str(REPO / 'lib'))
_spec = importlib.util.spec_from_file_location('refine_ambiguous_families', BIN)
raf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(raf)

from hits import Hit  # noqa: E402


def test_has_species_duplication():
    p2p = {'a': 'In1', 'b': 'In2', 'c': 'In1'}
    assert raf.has_species_duplication(['a', 'b', 'c'], p2p)  # In1 appears twice
    assert not raf.has_species_duplication(['a', 'b'], p2p)


def test_family_fractions():
    presence_by_protein = {'a': {'In1': 1, 'In2': 1, 'Out1': 1, 'Out2': 0}}
    ing, out = raf.family_fractions('a', ['a', 'b', 'c'], presence_by_protein,
                                    {'In1', 'In2'}, {'Out1', 'Out2'})
    assert ing == 1.0   # both In1, In2 present
    assert out == 0.5   # Out1 present, Out2 absent


def test_detect_ambiguous_families_matches_hex1_ada1_split():
    # HEX1-like family: species-duplicated (In1 x2) AND a near-miss novelty candidate
    # (ingroup_frac=1.0 >= 0.75, outgroup_frac=1.0 > 0.0) -> ambiguous.
    # ADA1-like family: no duplication -> never ambiguous regardless of presence.
    fam_members = {
        'hex1_rep': ['hex1', 'eif5a', 'orth_in2'],  # In1, In1, In2
        'ada1_rep': ['ada1_in1', 'ada1_in2'],        # In1, In2
    }
    p2p = {'hex1': 'In1', 'eif5a': 'In1', 'orth_in2': 'In2',
           'ada1_in1': 'In1', 'ada1_in2': 'In2'}
    ingroup_ids, outgroup_ids = {'In1', 'In2'}, {'Out1'}
    # Fake presence: both families present in all ingroup + all outgroup (Out1).
    # (family_fractions is computed from cluster membership + a presence lookup in the
    # real pipeline via presence_matrix.tsv; here we inject it directly via a stub.)
    def fake_fractions(rep, members, presence_by_protein, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, {}, ingroup_ids, outgroup_ids, oversized_reps=set(),
        ingroup_min_frac=0.75, other_max_frac=0.0)
    assert ambiguous == {'hex1_rep'}


def test_detect_ambiguous_families_excludes_oversized():
    fam_members = {'big_rep': ['a', 'b', 'c']}
    p2p = {'a': 'In1', 'b': 'In1', 'c': 'In2'}
    def fake_fractions(rep, members, presence_by_protein, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, {}, {'In1', 'In2'}, {'Out1'}, oversized_reps={'big_rep'},
        ingroup_min_frac=0.75, other_max_frac=0.0)
    assert ambiguous == set()


def test_extract_needed_sequences(tmp_path):
    fa1 = tmp_path / 'In1.pep.fa'
    fa1.write_text(">a desc one\nMKV\n>b desc two\nMKL\n")
    fa2 = tmp_path / 'In2.pep.fa'
    fa2.write_text(">c desc three\nMKA\n>d desc four\nMKD\n")
    seqs = raf.extract_needed_sequences({'In1': fa1, 'In2': fa2}, needed_ids={'a', 'd'})
    assert seqs == {'a': 'MKV', 'd': 'MKD'}


def test_write_fasta(tmp_path):
    out = tmp_path / 'out.fa'
    raf.write_fasta({'a': 'MKV', 'b': 'MKL'}, out)
    text = out.read_text()
    assert text == ">a\nMKV\n>b\nMKL\n"


def test_load_profiled_reps(tmp_path):
    families_tsv = tmp_path / 'families.tsv'
    families_tsv.write_text(
        'family_index\trepresentative_id\tn_members\n'
        'fam_000001\trep_a\t3\n'
        'fam_000002\trep_b\t2\n'
    )
    reps = raf.load_profiled_reps(families_tsv)
    assert reps == {'rep_a', 'rep_b'}


def test_filter_to_profiled():
    fam_members = {'rep_a': ['x', 'y'], 'rep_b': ['z'], 'rep_c': ['w']}
    profiled_reps = {'rep_a', 'rep_b'}
    filtered = raf.filter_to_profiled(fam_members, profiled_reps)
    assert filtered == {'rep_a': ['x', 'y'], 'rep_b': ['z']}


def test_write_refined_families_passthrough_and_split(tmp_path):
    fam_members = {'rep1': ['x', 'y'], 'hex1_rep': ['hex1', 'eif5a', 'orth_in2']}
    ambiguous = {'hex1_rep'}
    subfamilies = {'hex1_rep': [['hex1', 'orth_in2'], ['eif5a']]}
    out_cluster = tmp_path / 'refined_cluster.tsv'
    out_families = tmp_path / 'refined_families.tsv'
    raf.write_refined_families(fam_members, ambiguous, subfamilies, out_cluster, out_families)

    cluster_lines = set(out_cluster.read_text().splitlines())
    # rep1 passes through unchanged (not ambiguous).
    assert 'rep1\tx' in cluster_lines
    assert 'rep1\ty' in cluster_lines
    # hex1_rep is split: each subfamily keyed by its own first member as new rep.
    assert 'hex1\thex1' in cluster_lines
    assert 'hex1\torth_in2' in cluster_lines
    assert 'eif5a\teif5a' in cluster_lines
    assert 'hex1_rep\thex1' not in cluster_lines  # old rep gone for the split family

    families_lines = out_families.read_text().splitlines()
    assert families_lines[0] == 'family_index\trepresentative_id\tn_members'
    body = {line.split('\t')[1]: line.split('\t')[2] for line in families_lines[1:]}
    assert body['rep1'] == '2'
    assert body['hex1'] == '2'
    assert body['eif5a'] == '1'


# --------------------------------------------------------------- graph splitting

def test_load_paralog_map(tmp_path):
    cutoffs = tmp_path / 'Ncra.paralog_cutoffs.tsv'
    cutoffs.write_text(
        'protein_ID\tparalog_protein_ID\tbitscore\tevalue\n'
        'HEX1_NEUCR\tIF5A_NEUCR\t512.3\t1e-140\n'
        'IF5A_NEUCR\tHEX1_NEUCR\t512.3\t1e-140\n'
        'OTHER_PROT\tANOTHER_PROT\t88.0\t1e-20\n'
    )
    paralog_map = raf.load_paralog_map([cutoffs])
    assert paralog_map == {
        'HEX1_NEUCR': 'IF5A_NEUCR',
        'IF5A_NEUCR': 'HEX1_NEUCR',
        'OTHER_PROT': 'ANOTHER_PROT',
    }


def test_load_paralog_map_merges_multiple_files(tmp_path):
    f1 = tmp_path / 'Ncra.paralog_cutoffs.tsv'
    f1.write_text('protein_ID\tparalog_protein_ID\tbitscore\tevalue\nA\tB\t1.0\t1e-5\n')
    f2 = tmp_path / 'Afum.paralog_cutoffs.tsv'
    f2.write_text('protein_ID\tparalog_protein_ID\tbitscore\tevalue\nC\tD\t1.0\t1e-5\n')
    paralog_map = raf.load_paralog_map([f1, f2])
    assert paralog_map == {'A': 'B', 'C': 'D'}


def test_build_within_family_edges_excludes_self_hits_and_high_evalue():
    hits = [
        Hit(query_id='a', target_id='b', evalue=1e-13, bitscore=200.0),
        Hit(query_id='b', target_id='a', evalue=1e-13, bitscore=200.0),  # reverse dup
        Hit(query_id='a', target_id='a', evalue=1e-50, bitscore=999.0),  # self-hit
        Hit(query_id='a', target_id='c', evalue=1e-5, bitscore=50.0),   # AT cutoff -> excluded
        Hit(query_id='a', target_id='d', evalue=1e-4, bitscore=10.0),   # ABOVE cutoff -> excluded
    ]
    edges = raf.build_within_family_edges(hits, evalue_cutoff=1e-5)
    assert edges == {frozenset({'a', 'b'})}


def test_find_paralog_edges_to_cut_direct_pair():
    members = ['hex1', 'eif5a', 'other']
    paralog_map = {'hex1': 'eif5a', 'eif5a': 'hex1'}
    cuts = raf.find_paralog_edges_to_cut(members, paralog_map)
    assert cuts == {frozenset({'hex1', 'eif5a'})}


def test_find_paralog_edges_to_cut_ignores_partner_outside_family():
    # hex1's registered paralog is NOT a member of this family -> no cut edge.
    members = ['hex1', 'other']
    paralog_map = {'hex1': 'some_other_genome_paralog'}
    cuts = raf.find_paralog_edges_to_cut(members, paralog_map)
    assert cuts == set()


def test_connected_components_singleton_with_no_surviving_edges():
    members = ['a', 'b', 'c']
    edges = {frozenset({'a', 'b'})}
    components = raf.connected_components(members, edges)
    groups = {frozenset(c) for c in components}
    assert groups == {frozenset({'a', 'b'}), frozenset({'c'})}


def test_connected_components_stays_connected_via_redundant_path():
    # Real HEX1 case: cutting the direct hex1-eif5a edge still leaves them in one
    # component because both are also independently connected to a shared ortholog
    # ('bridge') -- the second path that made refinement fail to split on real data.
    members = ['hex1', 'eif5a', 'bridge']
    edges = {frozenset({'hex1', 'bridge'}), frozenset({'eif5a', 'bridge'})}  # direct edge already cut
    components = raf.connected_components(members, edges)
    assert len(components) == 1
    assert set(components[0]) == {'hex1', 'eif5a', 'bridge'}


def test_split_family_splits_when_only_direct_edge_connects_pair():
    members = ['hex1', 'eif5a']
    edges = {frozenset({'hex1', 'eif5a'})}
    paralog_map = {'hex1': 'eif5a', 'eif5a': 'hex1'}
    subfamilies = raf.split_family(members, edges, paralog_map)
    groups = {frozenset(g) for g in subfamilies}
    assert groups == {frozenset({'hex1'}), frozenset({'eif5a'})}


def test_split_family_does_not_split_with_redundant_bridging_edge():
    # Matches the real HEX1/IF5A_NEUCR finding (task-5-fix-report.md): the direct
    # registered-paralog edge is cut, but a redundant path through a shared ortholog
    # keeps the family in one component -- refinement does NOT split it.
    members = ['hex1', 'eif5a', 'bridge']
    edges = {
        frozenset({'hex1', 'eif5a'}),   # direct paralog edge -- will be cut
        frozenset({'hex1', 'bridge'}),  # redundant path
        frozenset({'eif5a', 'bridge'}),
    }
    paralog_map = {'hex1': 'eif5a', 'eif5a': 'hex1'}
    subfamilies = raf.split_family(members, edges, paralog_map)
    assert len(subfamilies) == 1
    assert set(subfamilies[0]) == {'hex1', 'eif5a', 'bridge'}
