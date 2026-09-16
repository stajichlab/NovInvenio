#!/usr/bin/env python3
"""Core/soft-core/shell/cloud/singleton frequency binning.

Bins are computed over the INGROUP strains only (samplesheet `GROUP ==
--ingroup_label`); counting outgroup reference strains in the denominator
would shift every frequency and therefore every bin boundary. Pass
`--inventory` (pangenome_dereplicate_strains.py's `strain_inventory.tsv`) to
further restrict the denominator to one representative per mash dedup
group.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/frequency_bins.py); `config_parser` is imported directly from this
repo's own lib/ (no cross-repo lookup needed), and `--ingroup_label`
replaces a hardcoded literal `"IN"` (NEXTFLOW_MIGRATION_NOTES.md section E).

Usage:
  pangenome_frequency_bins.py --matrix presence_matrix.rescued.tsv --config config.csv \\
      [--inventory strain_inventory.tsv] --output frequency_table.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import PresenceMatrix  # noqa: E402
from pangenome_strain_inventory import read_representative_shorts  # noqa: E402


def assign_bin(
    freq: float,
    strain_count: int,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
) -> str:
    if strain_count <= 1:
        return "singleton"
    if freq >= core_cutoff:
        return "core"
    if freq >= softcore_cutoff:
        return "soft_core"
    if freq >= shell_cutoff:
        return "shell"
    return "cloud"


def compute_frequency_table(
    matrix: PresenceMatrix,
    core_cutoff: float = 0.95,
    softcore_cutoff: float = 0.90,
    shell_cutoff: float = 0.15,
    strains: list[str] | None = None,
) -> list[dict]:
    """`strains` restricts both the numerator and the denominator of every
    frequency to that subset (ingroup-only, optionally dereplicated) instead
    of always using `matrix.strains` -- the same pattern as
    pangenome_cooccurrence.py's `find_cooccurring_pairs`. Defaults to
    `matrix.strains` for backward compatibility with callers whose matrix is
    already exactly the right set."""
    strains = matrix.strains if strains is None else strains
    n = len(strains)
    rows = []
    for fam in matrix.families:
        count = sum(1 for s in strains if matrix.is_present(fam, s))
        freq = count / n if n else 0.0
        rows.append({
            "family": fam,
            "frequency": freq,
            "strain_count": count,
            "bin": assign_bin(freq, count, core_cutoff, softcore_cutoff, shell_cutoff),
        })
    return rows


def resolve_strains(
    matrix: PresenceMatrix, config_path: str, inventory_path: str | None, ingroup_label: str,
) -> list[str]:
    """Ingroup strains from the samplesheet, intersected with the matrix's
    own columns and (when given) with the inventory's representative
    strains."""
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(config_path)
    strains = [s.short for s in samples if s.group == ingroup_label]
    if inventory_path:
        reps = set(read_representative_shorts(inventory_path))
        strains = [s for s in strains if s in reps]
    in_matrix = set(matrix.strains)
    missing = [s for s in strains if s not in in_matrix]
    if missing:
        print(
            f"WARNING: {len(missing)} ingroup strains are absent from the matrix "
            f"columns and are ignored (first few: {missing[:5]})",
            file=sys.stderr,
        )
    return [s for s in strains if s in in_matrix]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--config", required=True,
                    help="samplesheet CSV -- its ingroup strains are the frequency denominator")
    ap.add_argument("--ingroup_label", default="IN",
                    help="samplesheet GROUP value identifying ingroup strains (default: IN)")
    ap.add_argument("--inventory",
                    help="strain_inventory.tsv from pangenome_dereplicate_strains.py; when "
                         "given, only is_representative == 1 strains are counted")
    ap.add_argument("--core_cutoff", type=float, default=0.95)
    ap.add_argument("--softcore_cutoff", type=float, default=0.90)
    ap.add_argument("--shell_cutoff", type=float, default=0.15)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    strains = resolve_strains(matrix, args.config, args.inventory, args.ingroup_label)
    if not strains:
        print("ERROR: no ingroup strains left to bin over (check --config/--inventory "
              "against the matrix's columns)", file=sys.stderr)
        sys.exit(1)
    print(
        f"Binning over {len(strains)} ingroup strains"
        f"{' (dereplicated)' if args.inventory else ''}",
        file=sys.stderr,
    )

    table = compute_frequency_table(
        matrix, args.core_cutoff, args.softcore_cutoff, args.shell_cutoff,
        strains=strains,
    )
    with open(args.output, "w") as fh:
        fh.write("family\tfrequency\tstrain_count\tbin\n")
        for row in table:
            fh.write(f"{row['family']}\t{row['frequency']:.4f}\t{row['strain_count']}\t{row['bin']}\n")


if __name__ == "__main__":
    main()
