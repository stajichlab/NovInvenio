import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / 'bin' / 'refine_ambiguous_families.py'
sys.path.insert(0, str(REPO / 'lib'))
_spec = importlib.util.spec_from_file_location('refine_ambiguous_families', BIN)
raf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(raf)


def test_has_species_duplication():
    p2p = {'a': 'In1', 'b': 'In2', 'c': 'In1'}
    assert raf.has_species_duplication(['a', 'b', 'c'], p2p)  # In1 appears twice
    assert not raf.has_species_duplication(['a', 'b'], p2p)


def test_family_fractions():
    p2p = {'a': 'In1', 'b': 'In2', 'c': 'Out1'}
    ing, out = raf.family_fractions('a', ['a', 'b', 'c'], p2p, {'In1', 'In2'}, {'Out1', 'Out2'})
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
    def fake_fractions(rep, members, p2p, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, ingroup_ids, outgroup_ids, oversized_reps=set(),
        ingroup_min_frac=0.75, other_max_frac=0.0)
    assert ambiguous == {'hex1_rep'}


def test_detect_ambiguous_families_excludes_oversized():
    fam_members = {'big_rep': ['a', 'b', 'c']}
    p2p = {'a': 'In1', 'b': 'In1', 'c': 'In2'}
    def fake_fractions(rep, members, p2p, ingroup_ids, outgroup_ids):
        return (1.0, 1.0)
    raf.family_fractions = fake_fractions
    ambiguous = raf.detect_ambiguous_families(
        fam_members, p2p, {'In1', 'In2'}, {'Out1'}, oversized_reps={'big_rep'},
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
