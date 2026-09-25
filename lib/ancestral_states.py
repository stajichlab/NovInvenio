"""Dollo-parsimony ancestral state reconstruction for gene-family presence/absence.

Replaces nothing on its own -- bin/pangenome_cooccurrence.py's count-based
`polarize_direction()` keeps its current meaning and columns (issue #140). This module
supplies the phylogeny-aware alternative that is reported ALONGSIDE it, so the two can
be compared on a real run before anything is retired.

THE MODEL, AND ITS ASSUMPTION

Dollo parsimony assumes a gene family is gained EXACTLY ONCE and may then be lost many
times independently. The justification is asymmetry: independently evolving the same
complex gene twice is far less likely than losing it twice.

    1. Place the single gain on the branch leading to the MRCA of every tip that
       carries the family.
    2. Below that node, each MAXIMAL all-absent clade is ONE loss event, on its stem.
    3. A node with more than two children is treated as a SOFT polytomy (unknown
       order, e.g. collapsed for low support): its absent children together are ONE
       loss, the minimum-loss resolution (issue #185).

Step 2 is the point of using a tree at all. Two absent sister tips are ONE loss, not
two, because a single loss on their shared stem explains both. A raw count of absent
strains -- which is what a count-based rule implicitly reports -- says two.

WHERE THIS IS WRONG, AND WHY IT SHIPS ANYWAY

Dollo FORBIDS independent gains. For families acquired by horizontal transfer -- which
in a pangenome is plausible for exactly the accessory content that is usually the point
of the analysis -- repeated independent gain is real, and Dollo must explain it as one
ancient gain plus a pile of losses. This SYSTEMATICALLY INFLATES LOSS COUNTS in the most
interesting families, and no amount of tree quality fixes it; it is the model, not the
data.

It ships as the starting point because it is the standard baseline in the gene-family
gain/loss literature and needs no rate parameters to estimate. The intended successors
are Wagner parsimony with an explicit gain:loss cost ratio, or a likelihood model
(Count, BadiRate). `DolloResult.method` names the assumption so downstream output can
carry it rather than presenting Dollo loss counts as model-free fact.

ROOTING

Dollo cannot polarise an unrooted tree: without a root there is no ancestral direction,
so "the MRCA of the carriers" is undefined. A trifurcating root is the usual signature
of an unrooted tree and is rejected. A bifurcating root that was never actually rooted
cannot be detected here -- root the tree on the outgroup before calling this.
"""
from dataclasses import dataclass, field

__all__ = ['DolloResult', 'dollo_polarize', 'tree_direction']


@dataclass(frozen=True)
class DolloResult:
    """One family's reconstruction.

    gain_node     readable label for the branch carrying the single gain (None if the
                  family is absent from every tip).
    gain_tips     every tip descended from that node -- carriers AND the ones below it
                  that lost it. NOT the same as the set of carriers.
    n_loss_events number of INDEPENDENT losses, not the number of tips lacking it.
    loss_clades   tip sets that lost it, one tuple per event, sorted.
    method        the reconstruction assumption, for the caller to record.
    """
    gain_node: str | None
    gain_tips: tuple[str, ...]
    n_loss_events: int
    loss_clades: list[tuple[str, ...]] = field(default_factory=list)
    method: str = 'dollo'


def _tips(clade):
    return tuple(sorted(t.name for t in clade.get_terminals()))


def _label(clade, tips):
    """A stable, readable identifier. Internal nodes are usually unnamed in Newick."""
    if clade.name:
        return clade.name
    if len(tips) == 1:
        return tips[0]
    return f'MRCA({tips[0]}..{tips[-1]})'


def dollo_polarize(tree, present_tips) -> DolloResult:
    """Reconstruct one family's single gain and its independent losses.

    tree          a ROOTED Bio.Phylo tree.
    present_tips  set of tip names carrying the family.
    """
    if len(tree.root.clades) > 2:
        raise ValueError(
            'tree root has more than two children, which is the usual signature of an '
            'unrooted tree; Dollo parsimony needs a rooted tree to polarise gain vs loss'
        )

    all_tips = set(_tips(tree.root))
    present = set(present_tips)
    unknown = present - all_tips
    if unknown:
        raise ValueError(
            f'{len(unknown)} carrier tip(s) not in the tree: {sorted(unknown)[:5]} -- '
            'tree tip labels must match the presence matrix strain identifiers exactly'
        )

    if not present:
        return DolloResult(gain_node=None, gain_tips=(), n_loss_events=0)

    gain_clade = (tree.common_ancestor(*present) if len(present) > 1
                  else next(t for t in tree.root.get_terminals() if t.name in present))
    gain_tips = _tips(gain_clade)

    losses: list[tuple[str, ...]] = []

    def walk(clade):
        tips = _tips(clade)
        if not (set(tips) & present):
            # Maximal all-absent clade: one loss on this stem. Do NOT descend -- its
            # children are explained by this single event.
            losses.append(tips)
            return
        children = clade.clades
        if len(children) > 2:
            # Soft polytomy (e.g. a node collapsed for low support): the order of
            # its children is unknown. The minimum-loss resolution puts all absent
            # children in one clade, so together they are ONE loss (issue #185).
            absent = [c for c in children if not (set(_tips(c)) & present)]
            if absent:
                losses.append(tuple(sorted(t for c in absent for t in _tips(c))))
            children = [c for c in children if c not in absent]
        for child in children:
            walk(child)

    walk(gain_clade)

    return DolloResult(
        gain_node=_label(gain_clade, gain_tips),
        gain_tips=gain_tips,
        n_loss_events=len(losses),
        loss_clades=sorted(losses),
    )


def tree_direction(result: DolloResult, carriers, ingroup) -> str:
    """Label one family's ingroup history from its reconstruction.

    The count-based polarize_direction() decides from outgroup counts alone, which
    conflates two different histories: a family every ingroup strain still carries and
    one several ingroup strains have lost both look like "loss" to it whenever the
    outgroup is complete. Asking where the single gain sits separates them.

      gain           the gain is inside the ingroup -- the family arose there
      loss           the gain is ancestral (spans the outgroup) and some ingroup
                     strain lacks it
      conserved      ancestral and every ingroup strain still has it: no ingroup event
      outgroup_only  the gain sits entirely outside the ingroup
      absent         no strain carries it
    """
    carriers, ingroup = set(carriers), set(ingroup)
    if not carriers:
        return 'absent'

    gain_tips = set(result.gain_tips)
    if not gain_tips & ingroup:
        return 'outgroup_only'
    if gain_tips <= ingroup:
        return 'gain'
    # Ancestral gain: did anything in the ingroup lose it?
    return 'loss' if (ingroup - carriers) else 'conserved'
