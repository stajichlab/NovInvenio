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
