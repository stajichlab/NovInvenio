"""Column labels for the island and locus views (uniform within one island).

A gene family is identified everywhere by its mmseqs tier-1 representative, for example
`Asfu_K148s1|KAH2869786.1`. mmseqs picks that protein from whichever strain it likes, and a
second run picks differently (73% ID overlap between two runs of one study). Used as a
column label it names the same island by genes from several arbitrary strains.

Each island already has one display strain (the island's `example_strain`, or a locus's
`exemplar`), and `lib/island_synteny.column_locations()` resolves each column to that
strain's own gene for the hover popup. The visible label is that same gene's ID, so the
label, the popup and the coordinates all come from one strain. The representative ID stays
in the payload as the stable key for joining to tables, and is shown in the popup.

Where the display strain has no annotated gene for a family (absent, or only a TBLASTN
rescue hit), the label falls back to the representative ID. The page draws those labels in
italics, and the italic ID carries its own `Short|` prefix, so a reader can see it is
from another strain.
"""
from __future__ import annotations

LABEL_REF = "ref"   # the display strain's own gene ID
LABEL_REP = "rep"   # fallback: the family's mmseqs representative ID


def column_labels(families: list[str],
                  family_locations: list | None) -> tuple[list[str], list[str]]:
    """(labels, kinds), one entry per column, in `families` order.

    `family_locations` is column_locations()' result for the same columns, or None
    when no gene coordinates were available (every label is then the representative ID).
    """
    labels: list[str] = []
    kinds: list[str] = []
    for i, fam in enumerate(families):
        loc = family_locations[i] if family_locations and i < len(family_locations) else None
        protein = loc.get("protein") if isinstance(loc, dict) else None
        if protein:
            labels.append(str(protein))
            kinds.append(LABEL_REF)
        else:
            labels.append(fam)
            kinds.append(LABEL_REP)
    return labels, kinds
