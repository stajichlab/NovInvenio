#!/usr/bin/env python3
"""Classify each FDR-significant co-occurring family pair as physically
linked (captain-gene-explained or not) vs. trans (candidate non-physical
interaction/co-evolution).

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/pair_classification.py). `linkage_fraction` is imported from this
repo's lib/pangenome_synteny.py (a real Nextflow-repo module) instead of
the study's own sibling bin/synteny_windows.py CLI stub.

Labels (a pair gets exactly one):
  - insufficient_data:    fewer than `min_co_carrying` strains carry both
                           families with a resolvable genomic position.
  - starship_explained:   linkage_fraction >= `physical_threshold` AND at
                           least one co-carrying strain has a captain-gene
                           family (--captain_tblout) within `k` genes of
                           either family. Named for the fungal
                           Starship-transposon mechanism this was designed
                           around; the captain-gene evidence itself is
                           OPTIONAL (see --captain_tblout below) and not
                           fungal/Starship-specific in its mechanics -- any
                           "mobile element marker gene" hmmsearch tblout
                           works the same way.
  - unexplained_physical: linkage_fraction >= `physical_threshold`, no
                           captain-gene evidence found -- a real candidate
                           for an unannotated/novel mobile element, not
                           discarded as "expected".
  - ambiguous_linkage:    `trans_threshold` < linkage_fraction <
                           `physical_threshold` -- reported, not discarded.
  - trans:                linkage_fraction <= `trans_threshold` AND the
                           pair ALSO clears `permutation_p < perm_alpha`
                           AND spans >= `min_clades` distinct clade labels.
  - trans_unconfirmed:    linkage_fraction <= `trans_threshold` but fails
                           the permutation/clade gate above -- physically
                           non-linked but not statistically robust enough
                           to call a confident interaction candidate.

Known partial-implementation caveats (carried over from the study,
recorded here rather than hidden):
  - `linkage_fraction` is rank-window only ("within k genes" statistic) --
    an additional same-accessory-island check and a real bp-distance window
    sized to an actual mobile-element footprint is NOT yet implemented. A
    physical block much larger than `k` genes is not yet caught.
  - No `taxon_scheme` restriction: `min_clades` counts distinct entries in
    `clade_composition` as given, which may mix incompatible
    population-clustering schemes if a study's TaxonGroup column does.
    Treat `trans` calls as provisional until a caller restricts to one
    consistent clade scheme.

Usage:
  pangenome_pair_classification.py --cooccurring_pairs cooccurring_pairs.tsv \\
      --family_positions family_positions.tsv \\
      --cluster_tsv tier1_cluster.tsv \\
      --captain_tblout captain_vs_study.tblout \\
      --output pair_classification.tsv

`--captain_tblout` is OPTIONAL (pass an empty file, e.g. from
EMPTY_EVALUES_STUB, when a study has no captain-gene/mobile-element marker
to search for) -- every pair then classifies as `unexplained_physical`
rather than `starship_explained` when linkage_fraction clears
`physical_threshold`, since there is no captain-gene evidence to check.
"""
from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import read_cluster_tsv, iter_tblout_family_hits  # noqa: E402
from compressed_io import open_maybe_compressed  # noqa: E402
from pangenome_synteny import linkage_fraction  # noqa: E402

DEFAULT_K = 10
DEFAULT_PHYSICAL_THRESHOLD = 0.5
DEFAULT_TRANS_THRESHOLD = 0.05
DEFAULT_MIN_CO_CARRYING = 5
DEFAULT_PERM_ALPHA = 0.05
DEFAULT_MIN_CLADES = 2


def load_family_positions(path: str) -> dict[str, dict[str, list[tuple[str, int]]]]:
    gene_position: dict[str, dict[str, list[tuple[str, int]]]] = {}
    with open_maybe_compressed(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            short, family, contig, rank = parts[0], parts[1], parts[2], int(parts[3])
            gene_position.setdefault(short, {}).setdefault(family, []).append((contig, rank))
    return gene_position


def load_captain_families(
    tblout_path: str, member_to_rep: dict[str, str], id_sep: str = "|"
) -> dict[str, set[str]]:
    """{Short: {family_id, ...}} for every strain with at least one
    captain/mobile-element-marker hmmsearch hit, resolved to that protein's
    OWN tier-1 family -- every captain-hit protein was part of the same
    all-strains clustering input, so it always has a family assignment
    already. An empty/missing tblout yields {} (no captain evidence
    anywhere), the documented behavior for a study with no marker gene.
    Parsing itself lives in lib/pangenome_matrix.iter_tblout_family_hits
    (shared with pangenome_build_islands.py's load_hit_families)."""
    captain_families: dict[str, set[str]] = {}
    for short, family in iter_tblout_family_hits(tblout_path, member_to_rep, id_sep=id_sep):
        captain_families.setdefault(short, set()).add(family)
    return captain_families


def has_captain_evidence(
    family_a: str,
    family_b: str,
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    captain_families: dict[str, set[str]],
    k: int,
) -> bool:
    """True if, in ANY strain, a captain-gene family sits within k genes of
    family_a or family_b (reuses linkage_fraction's own proximity test by
    treating each captain family as one more family to check linkage
    against)."""
    for short, captains in captain_families.items():
        positions = gene_position.get(short, {})
        if family_a not in positions and family_b not in positions:
            continue
        for captain_family in captains:
            if captain_family not in positions:
                continue
            single_strain_positions = {short: positions}
            for fam in (family_a, family_b):
                if fam not in positions:
                    continue
                if linkage_fraction(fam, captain_family, single_strain_positions, k=k) > 0:
                    return True
    return False


def classify_pair(
    family_a: str,
    family_b: str,
    permutation_p: float,
    clade_composition: dict[str, int],
    gene_position: dict[str, dict[str, list[tuple[str, int]]]],
    captain_families: dict[str, set[str]],
    k: int = DEFAULT_K,
    physical_threshold: float = DEFAULT_PHYSICAL_THRESHOLD,
    trans_threshold: float = DEFAULT_TRANS_THRESHOLD,
    min_co_carrying: int = DEFAULT_MIN_CO_CARRYING,
    perm_alpha: float = DEFAULT_PERM_ALPHA,
    min_clades: int = DEFAULT_MIN_CLADES,
) -> tuple[str, float]:
    """Returns (classification, linkage_fraction)."""
    n_co_carrying = sum(
        1 for positions in gene_position.values()
        if positions.get(family_a) and positions.get(family_b)
    )
    if n_co_carrying < min_co_carrying:
        return "insufficient_data", 0.0

    frac = linkage_fraction(family_a, family_b, gene_position, k=k)

    if frac >= physical_threshold:
        if has_captain_evidence(family_a, family_b, gene_position, captain_families, k):
            return "starship_explained", frac
        return "unexplained_physical", frac

    if frac <= trans_threshold:
        n_clades = len(clade_composition)
        if permutation_p < perm_alpha and n_clades >= min_clades:
            return "trans", frac
        return "trans_unconfirmed", frac

    return "ambiguous_linkage", frac


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cooccurring_pairs", required=True)
    ap.add_argument("--family_positions", required=True)
    ap.add_argument("--cluster_tsv", required=True)
    ap.add_argument("--captain_tblout", required=True,
                    help="hmmsearch --tblout for an optional captain/mobile-element "
                         "marker gene; pass an empty file if the study has none")
    ap.add_argument("--k", type=int, default=DEFAULT_K)
    ap.add_argument("--physical_threshold", type=float, default=DEFAULT_PHYSICAL_THRESHOLD)
    ap.add_argument("--trans_threshold", type=float, default=DEFAULT_TRANS_THRESHOLD)
    ap.add_argument("--min_co_carrying", type=int, default=DEFAULT_MIN_CO_CARRYING)
    ap.add_argument("--perm_alpha", type=float, default=DEFAULT_PERM_ALPHA)
    ap.add_argument("--min_clades", type=int, default=DEFAULT_MIN_CLADES)
    ap.add_argument("--id_sep", default="|",
                    help="Short-prefix separator in clustering-input FASTA headers "
                         "(default: '|'); must match --pangenome_id_sep used "
                         "everywhere else in this pipeline run.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    print("Loading family positions...", file=sys.stderr)
    gene_position = load_family_positions(args.family_positions)
    print("Loading cluster membership...", file=sys.stderr)
    member_to_rep = read_cluster_tsv(args.cluster_tsv)
    print("Loading captain-gene evidence...", file=sys.stderr)
    captain_families = load_captain_families(args.captain_tblout, member_to_rep, id_sep=args.id_sep)
    print(
        f"{sum(len(v) for v in captain_families.values())} captain-gene family "
        f"assignments across {len(captain_families)} strains", file=sys.stderr,
    )

    with open_maybe_compressed(args.cooccurring_pairs) as fh, open(args.output, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        out.write(
            "family_a\tfamily_b\tclassification\tlinkage_fraction\tjaccard\t"
            "fisher_p\tfdr_q\tpermutation_p\tdirection_a\tclade_composition\n"
        )
        n = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            family_a, family_b = parts[idx["family_a"]], parts[idx["family_b"]]
            permutation_p = float(parts[idx["permutation_p"]])
            clade_composition = ast.literal_eval(parts[idx["clade_composition"]])
            classification, frac = classify_pair(
                family_a, family_b, permutation_p, clade_composition,
                gene_position, captain_families,
                k=args.k, physical_threshold=args.physical_threshold,
                trans_threshold=args.trans_threshold,
                min_co_carrying=args.min_co_carrying,
                perm_alpha=args.perm_alpha, min_clades=args.min_clades,
            )
            out.write(
                f"{family_a}\t{family_b}\t{classification}\t{frac:.4f}\t"
                f"{parts[idx['jaccard']]}\t{parts[idx['fisher_p']]}\t"
                f"{parts[idx['fdr_q']]}\t{parts[idx['permutation_p']]}\t"
                f"{parts[idx['direction_a']]}\t{parts[idx['clade_composition']]}\n"
            )
            n += 1
            if n % 100_000 == 0:
                print(f"pair_classification: classified {n} pairs", file=sys.stderr)

    print(f"pair_classification: {n} pairs classified", file=sys.stderr)


if __name__ == "__main__":
    main()
