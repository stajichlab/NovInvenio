"""Per-strain gene-order synteny helpers for the pangenome-profiling
subworkflow's pair-classification step (physical linkage / Starship
mechanism between co-occurring gene families).

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/synteny_windows.py) per its own migration note: this module's CLI was
always a deliberate stub ("use these functions directly for a real run"),
so it belongs in `lib/` as importable logic, not as its own Nextflow
process (see NEXTFLOW_MIGRATION_NOTES.md section A item 10).

`linkage_fraction` is the one function `bin/pangenome_pair_classification.py`
actually calls today; `parse_gff3_gene_order`/`sliding_windows`/
`accessory_islands` are carried over unused by any process yet (a fuller
bp-distance/accessory-island physical-linkage test is still a documented
TODO -- see pangenome_pair_classification.py's module docstring), but kept
here rather than dropped so a future process can reach them without another
port.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

# Anchored to the start of the attributes string or right after a ';' so a
# decoy attribute like "orig_protein_ID=XP_1;ID=geneA" does not match on
# "ID=" buried inside "orig_protein_ID=" -- unanchored .search() would
# incorrectly return "XP_1" there.
_ID_RE = re.compile(r"(?:^|;)ID=([^;\n]+)")


def parse_gff3_gene_order(gff3_path: str | Path) -> list[tuple[str, str, int, int]]:
    """Return [(gene_id, contig_id, start, end), ...] for every 'gene'
    feature, sorted by contig then start position."""
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            contig, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            m = _ID_RE.search(attrs)
            if not m:
                continue
            genes.append((m.group(1), contig, start, end))
    genes.sort(key=lambda g: (g[1], g[2]))
    return genes


def sliding_windows(gene_order: list[tuple[str, str, int, int]], n: int) -> list[list[tuple]]:
    windows = []
    for i in range(len(gene_order) - n + 1):
        window = gene_order[i:i + n]
        contigs = {g[1] for g in window}
        if len(contigs) == 1:
            windows.append(window)
    return windows


def accessory_islands(
    gene_order: list[tuple[str, str, int, int]], is_core: dict[str, bool]
) -> list[list[tuple]]:
    """Maximal runs of consecutive non-core genes, never crossing a
    contig boundary.

    A gene_id absent from `is_core` is treated as core (`is_core.get(gene_id,
    True)`), i.e. it splits/ends any in-progress island rather than being
    folded into one. This is deliberately conservative: an unrecognized gene
    is more likely a family the clustering step never called than a genuine,
    silent accessory gene, so defaulting it to non-core would risk
    manufacturing spurious islands out of noise. Warns (via `warnings.warn`)
    whenever more than 10% of the genes it sees are absent from `is_core` --
    a real run should treat that warning as a hard stop, since it usually
    means a GFF3-`ID=` vs. presence-matrix protein-ID format mismatch.
    """
    islands: list[list[tuple]] = []
    current: list[tuple] = []
    prev_contig = None
    missing = 0
    for gene in gene_order:
        gene_id, contig = gene[0], gene[1]
        if gene_id not in is_core:
            missing += 1
        core = is_core.get(gene_id, True)
        if contig != prev_contig and current:
            islands.append(current)
            current = []
        if not core:
            current.append(gene)
        elif current:
            islands.append(current)
            current = []
        prev_contig = contig
    if current:
        islands.append(current)
    if gene_order and missing / len(gene_order) > 0.1:
        warnings.warn(
            f"accessory_islands: {missing}/{len(gene_order)} genes "
            f"({missing / len(gene_order):.0%}) were absent from is_core and "
            "defaulted to core -- check for a GFF3 ID vs presence-matrix "
            "protein ID format mismatch before trusting these islands.",
            stacklevel=2,
        )
    return islands


def linkage_fraction(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    k: int = 10,
) -> float:
    """Fraction of strains carrying BOTH families where SOME copy of
    family_a is within k genes of SOME copy of family_b on the same contig.

    `gene_position[strain][family]` is a LIST of (contig, rank) positions,
    not a single position: a multi-copy family (e.g. an ancient paralogous
    domain family present at several copies per strain) would otherwise have
    its physical linkage mislabeled by whichever copy happened to be stored.
    Checking the MINIMUM distance over every copy-pair answers "is ANY copy
    of family_a near ANY copy of family_b", not "is the arbitrary stored
    copy of each near the other".
    """
    both_present = [
        s for s, positions in gene_position.items()
        if positions.get(family_a) and positions.get(family_b)
    ]
    if not both_present:
        return 0.0
    linked = 0
    for s in both_present:
        positions = gene_position[s]
        is_linked = any(
            contig_a == contig_b and abs(pos_a - pos_b) <= k
            for contig_a, pos_a in positions[family_a]
            for contig_b, pos_b in positions[family_b]
        )
        if is_linked:
            linked += 1
    return linked / len(both_present)
