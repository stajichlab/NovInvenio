"""Dollo-parsimony gain/loss reconstruction (issue #140).

Every test builds its tree from a Newick string by hand, so these are real assertions
about the algorithm with no pipeline or data dependency.

The rule under test: a family is gained ONCE, on the branch leading to the MRCA of all
tips that carry it; below that node, every MAXIMAL all-absent clade is ONE loss event.
The maximal-clade part is the whole point -- two absent sister tips are one loss, not two.
"""
import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from ancestral_states import dollo_polarize  # noqa: E402

from Bio import Phylo  # noqa: E402


def tree(newick):
    return Phylo.read(io.StringIO(newick), 'newick')


BALANCED = '((A,B),(C,D));'
NESTED = '(((A,B),(C,D)),(E,F));'


def test_family_in_every_tip_is_gained_at_the_root_with_no_losses():
    r = dollo_polarize(tree(BALANCED), {'A', 'B', 'C', 'D'})
    assert r.gain_tips == ('A', 'B', 'C', 'D')
    assert r.n_loss_events == 0


def test_family_in_one_clade_is_gained_on_that_clades_stem():
    r = dollo_polarize(tree(BALANCED), {'A', 'B'})
    assert r.gain_tips == ('A', 'B')
    assert r.n_loss_events == 0


def test_single_absent_tip_is_one_loss():
    r = dollo_polarize(tree(BALANCED), {'A', 'B', 'C'})
    assert r.gain_tips == ('A', 'B', 'C', 'D')
    assert r.n_loss_events == 1
    assert r.loss_clades == [('D',)]


def test_two_absent_tips_in_different_clades_are_two_losses():
    r = dollo_polarize(tree(BALANCED), {'A', 'C'})
    assert r.n_loss_events == 2
    assert r.loss_clades == [('B',), ('D',)]


def test_two_absent_sister_tips_are_ONE_loss_not_two():
    # The defining behaviour of the method. C and D are sisters and both absent, so a
    # single loss on their shared stem explains both -- a raw count of absent strains
    # would say two, which is what the current count-based polarize_direction implies.
    r = dollo_polarize(tree(NESTED), {'A', 'B', 'E', 'F'})
    assert r.n_loss_events == 1
    assert r.loss_clades == [('C', 'D')]


def test_a_single_carrier_is_gained_on_its_own_tip_branch():
    r = dollo_polarize(tree(BALANCED), {'A'})
    assert r.gain_tips == ('A',)
    assert r.n_loss_events == 0


def test_absent_everywhere_has_no_gain_and_no_losses():
    r = dollo_polarize(tree(BALANCED), set())
    assert r.gain_tips == ()
    assert r.gain_node is None
    assert r.n_loss_events == 0


def test_nested_loss_below_a_deeper_gain():
    # Gain is at MRCA(A,B,C,D) -- E,F never had it, so they are NOT losses. Within that
    # subtree, D alone is absent: exactly one loss. A method that counted absent tips
    # across the whole tree would wrongly report three.
    r = dollo_polarize(tree(NESTED), {'A', 'B', 'C'})
    assert r.gain_tips == ('A', 'B', 'C', 'D')
    assert r.n_loss_events == 1
    assert r.loss_clades == [('D',)]


def test_tips_not_present_in_the_tree_are_rejected():
    with pytest.raises(ValueError, match='not in the tree'):
        dollo_polarize(tree(BALANCED), {'A', 'ZZZ'})


def test_unrooted_multifurcating_root_is_rejected():
    # Dollo cannot polarise without a root; a trifurcation at the root is the standard
    # signature of an unrooted tree and must fail loudly rather than guess.
    with pytest.raises(ValueError, match='rooted'):
        dollo_polarize(tree('(A,B,C);'), {'A'})


# --- turning a reconstruction into a gain/loss label -------------------------
#
# The count-based polarize_direction() answers from outgroup counts alone. The tree
# version asks where the single gain sits relative to the ingroup.

from ancestral_states import tree_direction  # noqa: E402

INGROUP = {'A', 'B', 'C', 'D'}
OUT = {'E', 'F'}


def test_gain_confined_to_the_ingroup_is_a_gain():
    r = dollo_polarize(tree(NESTED), {'A', 'B'})
    assert tree_direction(r, carriers={'A', 'B'}, ingroup=INGROUP) == 'gain'


def test_ancestral_gain_with_an_ingroup_strain_lacking_it_is_a_loss():
    # Gain spans the outgroup, so the family predates the ingroup; D lacks it.
    r = dollo_polarize(tree(NESTED), {'A', 'B', 'C', 'E', 'F'})
    assert tree_direction(r, carriers={'A', 'B', 'C', 'E', 'F'},
                          ingroup=INGROUP) == 'loss'


def test_ancestral_gain_retained_by_every_ingroup_strain_is_conserved():
    # Neither gained nor lost within the ingroup -- the count-based rule calls this
    # "loss" whenever the outgroup is complete, which is wrong when nothing was lost.
    r = dollo_polarize(tree(NESTED), set('ABCDEF'))
    assert tree_direction(r, carriers=set('ABCDEF'), ingroup=INGROUP) == 'conserved'


def test_family_absent_from_every_strain_is_absent():
    r = dollo_polarize(tree(NESTED), set())
    assert tree_direction(r, carriers=set(), ingroup=INGROUP) == 'absent'


def test_gain_confined_to_the_outgroup_is_not_an_ingroup_event():
    r = dollo_polarize(tree(NESTED), {'E', 'F'})
    assert tree_direction(r, carriers={'E', 'F'}, ingroup=INGROUP) == 'outgroup_only'


# ---- soft polytomies (issue #185) --------------------------------------------------
# A node collapsed for low support is a SOFT polytomy: the true order is unknown. The
# minimum-loss resolution groups all absent children into one clade, so they are ONE
# loss, not one per child. Without this, collapsing poorly supported nodes (the point
# of doing it) would inflate loss counts instead of making them conservative.

def test_absent_children_of_a_polytomy_are_one_loss():
    r = dollo_polarize(tree('((A,B,C,D),E);'), {'A', 'B', 'E'})
    assert r.n_loss_events == 1
    assert r.loss_clades == [('C', 'D')]


def test_one_absent_child_of_a_polytomy_is_one_loss():
    r = dollo_polarize(tree('((A,B,C),D);'), {'A', 'B', 'D'})
    assert r.loss_clades == [('C',)]


def test_polytomy_absent_children_still_count_separately_from_other_losses():
    r = dollo_polarize(tree('(((A,B,C,D),(E,F)),G);'), {'A', 'B', 'E', 'G'})
    assert r.n_loss_events == 2
    assert r.loss_clades == [('C', 'D'), ('F',)]


def test_bifurcating_behaviour_is_unchanged():
    r = dollo_polarize(tree(NESTED), {'A', 'E'})
    assert r.n_loss_events == 3
