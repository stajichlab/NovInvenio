#!/usr/bin/env python3
"""Build the family x strain presence matrix from a tier-1 cluster TSV --
the producer of the `presence_matrix.tsv` that pangenome_rescue_pass.py,
pangenome_frequency_bins.py and pangenome_cooccurrence.py all consume.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/build_presence_matrix.py); unlike that copy, `config_parser` is
imported directly from THIS repo's own lib/ (no cross-repo NOVINVENIO_ROOT
lookup needed -- this script already lives next to lib/config_parser.py).

Usage:
  pangenome_build_presence_matrix.py --cluster_tsv tier1_cluster.tsv --config config.csv \\
      --output presence_matrix.tsv

## ID conventions established here (several other pangenome_* scripts depend on them)

1. **Protein IDs must be Short-prefixed BEFORE clustering** (see
   pangenome_prefix_fasta.py) -- the cluster TSV carries only sequence IDs,
   so the only way to recover which strain a family member came from is to
   put the strain there: ``><Short><id_sep><original_protein_id>``.

2. **Family ID = the tier-1 cluster representative ID, as-is.** No
   rewriting, no renumbering -- pangenome_rescue_pass.py's `--matrix` family
   column and pangenome_cooccurrence.py's family column both assume the
   family ID is literally the representative sequence ID mmseqs/diamond
   reported (Short prefix included).

3. **Columns = `config.csv`'s ingroup/outgroup strains by default**
   (`--groups IN,OUT`, or whatever `--ingroup_label`/`--outgroup_label`
   values a study's samplesheet uses). pangenome_cooccurrence.py's
   gain/loss polarization reads outgroup columns out of this same matrix --
   an ingroup-only matrix (`--groups IN`) makes every family's polarization
   "ambiguous" rather than a mis-reported confident gain/loss.

4. **`copy_number` = how many of that strain's proteins fall in that
   family.** Persisted alongside the matrix in `<output>.copy_number.tsv`
   (see lib/pangenome_matrix.py).

A member whose prefix is not a known strain is counted and reported; if
EVERY member is unrecognized the run hard-errors, because that means the
FASTA was not Short-prefixed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import (  # noqa: E402
    PRESENT,
    PresenceMatrix,
    build_families,
    read_cluster_tsv,
)

DEFAULT_ID_SEP = "|"


def split_member_id(member_id: str, id_sep: str = DEFAULT_ID_SEP) -> tuple[str, str]:
    """Split a Short-prefixed member ID into (strain_short, protein_id).

    Splits on the FIRST separator only, so a protein ID that itself contains
    the separator survives intact. An ID with no separator yields
    ``("", member_id)`` -- an unrecognized strain, handled by the caller.
    """
    strain, sep, protein = member_id.partition(id_sep)
    if not sep:
        return "", member_id
    return strain, protein


def count_family_members(
    families: dict[str, list[str]],
    strains: list[str],
    id_sep: str = DEFAULT_ID_SEP,
) -> tuple[dict[tuple[str, str], int], int]:
    """Count each strain's members per family.

    Returns ``(counts, n_unrecognized)`` where `counts` maps
    (family, strain) -> number of that strain's proteins in that family, and
    `n_unrecognized` counts members whose prefix is not in `strains`.
    """
    known = set(strains)
    counts: dict[tuple[str, str], int] = {}
    unrecognized = 0
    for family, members in families.items():
        for member in members:
            strain, _protein = split_member_id(member, id_sep)
            if strain not in known:
                unrecognized += 1
                continue
            key = (family, strain)
            counts[key] = counts.get(key, 0) + 1
    return counts, unrecognized


def build_matrix(
    families: dict[str, list[str]],
    strains: list[str],
    id_sep: str = DEFAULT_ID_SEP,
) -> tuple[PresenceMatrix, int]:
    """Build the PresenceMatrix for `families` over `strains`.

    Every family is a row (singletons included -- they are a real frequency
    bin here), every strain is a column. A (family, strain) cell is PRESENT
    with `copies` = that strain's member count in the family; strains with no
    member are left at the default ABSENT (the genome-level rescue pass, not
    this script, is what may later upgrade one to GENOME_ONLY).

    Returns ``(matrix, n_unrecognized_members)``.
    """
    counts, unrecognized = count_family_members(families, strains, id_sep)
    matrix = PresenceMatrix(families=sorted(families), strains=list(strains))
    for (family, strain), copies in counts.items():
        matrix.set_call(family, strain, PRESENT, copies=copies)
    return matrix, unrecognized


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cluster_tsv", required=True,
                    help="tier-1 cluster TSV (rep<TAB>member) from pangenome_cluster_backend.py")
    ap.add_argument("--config", required=True, help="samplesheet CSV (GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup)")
    ap.add_argument("--groups", default="IN,OUT",
                    help="comma-separated GROUP values to use as matrix columns "
                         "(default: IN,OUT -- pangenome_cooccurrence.py's gain/loss "
                         "polarization needs outgroup columns present; pass 'IN' to "
                         "build an ingroup-only matrix for a study with no outgroup)")
    ap.add_argument("--id_sep", default=DEFAULT_ID_SEP,
                    help=f"Short-prefix separator in member IDs (default: {DEFAULT_ID_SEP!r})")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    from config_parser import parse_config  # noqa: E402

    wanted = {g.strip() for g in args.groups.split(",") if g.strip()}
    samples = parse_config(args.config)
    strains = [s.short for s in samples if s.group in wanted]
    if not strains:
        print(f"ERROR: no strains in {args.config} with GROUP in {sorted(wanted)}",
              file=sys.stderr)
        sys.exit(1)

    families = build_families(read_cluster_tsv(args.cluster_tsv))
    matrix, unrecognized = build_matrix(families, strains, args.id_sep)

    total_members = sum(len(m) for m in families.values())
    print(
        f"Presence matrix: {len(matrix.families)} families x {len(strains)} strains "
        f"from {total_members} clustered proteins "
        f"({unrecognized} members skipped: prefix not a known strain)",
        file=sys.stderr,
    )
    if total_members and unrecognized == total_members:
        print(
            f"ERROR: all {total_members} cluster members had an unrecognized "
            f"'<Short>{args.id_sep}<protein_id>' prefix -- the clustering input FASTA "
            "was probably not Short-prefixed (see this script's docstring)",
            file=sys.stderr,
        )
        sys.exit(1)

    matrix.to_tsv(args.output)


if __name__ == "__main__":
    main()
