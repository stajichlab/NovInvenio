"""Lineage-proximity ranking shared by the master-pool builder and the
targeted config renderer. See docs/superpowers/specs/2026-09-05-config-builder-design.md
("Lineage-proximity ranking") for the full rationale.
"""

RANK_NAMES = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']


def lineage_match(lineage_a: list[str], lineage_b: list[str]) -> str:
    """Deepest rank NAME both lineages still agree on.

    Both lineages must be exactly len(RANK_NAMES) long (one slot per rank,
    '' for a rank that isn't recorded). A rank where either side is '' is
    skipped -- it neither extends nor breaks the match -- so a missing
    SUBPHYLUM/SUBCLASS never causes a same-genus pair to be misjudged as
    merely same-order. Comparison is case-insensitive. Returns '' if no
    rank matches (including if they diverge at PHYLUM itself).
    """
    if len(lineage_a) != len(RANK_NAMES) or len(lineage_b) != len(RANK_NAMES):
        raise ValueError(
            f"lineage must have exactly {len(RANK_NAMES)} fields "
            f"(got {len(lineage_a)} and {len(lineage_b)})"
        )
    deepest = ''
    for name, a, b in zip(RANK_NAMES, lineage_a, lineage_b):
        if not a or not b:
            continue
        if a.lower() != b.lower():
            break
        deepest = name
    return deepest
