"""--species_tree integration for bin/pangenome_cooccurrence.py (issue #140).

Drives the real find_cooccurring_pairs() with a hand-built matrix and tree, so the
wiring itself is under test: that carriers are taken from the FULL matrix (not the
ingroup-only strain list every other statistic here uses), that the per-family Dollo
result is memoised across pairs, and that the tree columns are added WITHOUT disturbing
the count-based direction_a.
"""
import importlib.util
import io
import sys
from pathlib import Path

from Bio import Phylo

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'lib'))
from pangenome_matrix import PresenceMatrix  # noqa: E402

spec = importlib.util.spec_from_file_location(
    'pangenome_cooccurrence', REPO / 'bin' / 'pangenome_cooccurrence.py')
coocc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coocc)

# Ingroup A-D, outgroup E-F. Families X and Y co-occur perfectly in the ingroup.
INGROUP = ['A', 'B', 'C', 'D']
ALL_STRAINS = INGROUP + ['E', 'F']
NEWICK = '(((A,B),(C,D)),(E,F));'


def build_matrix(present_by_family):
    m = PresenceMatrix(families=list(present_by_family), strains=list(ALL_STRAINS))
    for fam, carriers in present_by_family.items():
        for s in carriers:
            m.set_call(fam, s, 'present')
    return m


def run(present_by_family, tree=True):
    matrix = build_matrix(present_by_family)
    freq = [{'family': f, 'bin': 'shell'} for f in present_by_family]
    outgroup_presence = {
        f: (sum(1 for s in ('E', 'F') if matrix.is_present(f, s)), 2)
        for f in present_by_family
    }
    species_tree = Phylo.read(io.StringIO(NEWICK), 'newick') if tree else None
    return coocc.find_cooccurring_pairs(
        matrix, freq, {s: 'c1' for s in ALL_STRAINS}, outgroup_presence,
        min_strain_count=2, fdr_alpha=1.0, strains=INGROUP, screen_alpha=1.0,
        species_tree=species_tree,
    )


# X and Y are both ingroup-only and co-occur in A,B -- enough to survive the filters.
INGROUP_ONLY = {'X': {'A', 'B'}, 'Y': {'A', 'B'}}
# X is ancestral (spans the outgroup) and C,D lost it.
ANCESTRAL_LOST = {'X': {'A', 'B', 'E', 'F'}, 'Y': {'A', 'B'}}


def test_no_tree_produces_no_tree_columns():
    rows = run(INGROUP_ONLY, tree=False)
    assert rows, 'fixture produced no pairs'
    assert 'direction_a_tree' not in rows[0]
    assert 'n_loss_events' not in rows[0]


def test_ingroup_only_family_is_called_a_gain_by_the_tree():
    rows = run(INGROUP_ONLY)
    assert rows[0]['direction_a_tree'] == 'gain'
    assert rows[0]['n_loss_events'] == 0
    assert rows[0]['asr_method'] == 'dollo'


def test_carriers_come_from_the_full_matrix_not_the_ingroup_subset():
    # X's carriers include E and F. If the wiring used the ingroup-only strain list,
    # the gain would look confined to the ingroup and be miscalled "gain"; taking
    # carriers from the full matrix makes it an ancestral family that C,D lost.
    rows = run(ANCESTRAL_LOST)
    assert rows[0]['direction_a_tree'] == 'loss'
    assert rows[0]['n_loss_events'] == 1
    assert rows[0]['loss_clades'] == 'C,D'


def test_count_based_direction_is_untouched_by_the_tree():
    # The whole point of issue #140's additive scope: direction_a keeps its old value.
    with_tree = run(ANCESTRAL_LOST)[0]
    without = run(ANCESTRAL_LOST, tree=False)[0]
    assert with_tree['direction_a'] == without['direction_a']
