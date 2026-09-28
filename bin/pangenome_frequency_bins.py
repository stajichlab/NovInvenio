#!/usr/bin/env python3
"""Core/soft-core/shell/cloud/singleton frequency binning.

Bins are computed per group (issue #212): over the INGROUP strains
(samplesheet `GROUP == --ingroup_label`) for `bin`, and over the OUTGROUP
strains for `bin_out`; one group's strains never enter the other group's
denominator. A family counted in no representative of a group gets
`nonrep_only` (present only in that group's non-representative strains),
`outgroup_only` / `ingroup_only` (present only in the other group) or
`absent`, instead of `singleton`. Pass
`--inventory` (pangenome_dereplicate_strains.py's `strain_inventory.tsv`) to
further restrict the denominator to one representative per mash dedup
group.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/frequency_bins.py); `config_parser` is imported directly from this
repo's own lib/ (no cross-repo lookup needed), and `--ingroup_label`
replaces a hardcoded literal `"IN"` (NEXTFLOW_MIGRATION_NOTES.md section E).

Usage:
  pangenome_frequency_bins.py --matrix presence_matrix.rescued.tsv --config config.csv \\
      [--inventory strain_inventory.tsv] [--outgroup_label OUT --outgroup_min_bin_strains 3] \\
      --output frequency_table.tsv
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
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


def group_class(
    count: int, n_reps: int, in_nonrep: bool, in_other_group: bool, is_ingroup: bool,
    core_cutoff: float = 0.95, softcore_cutoff: float = 0.90, shell_cutoff: float = 0.15,
) -> str:
    """Class of one family within one group (spec section 1), rules in order:
    counted in >= 1 representative -> the frequency bin; else present in a
    non-representative of this group -> nonrep_only; else present in the
    other group -> outgroup_only / ingroup_only; else absent."""
    if count >= 1:
        return assign_bin(count / n_reps, count, core_cutoff, softcore_cutoff, shell_cutoff)
    if in_nonrep:
        return "nonrep_only"
    if in_other_group:
        return "outgroup_only" if is_ingroup else "ingroup_only"
    return "absent"


@dataclass
class GroupStrains:
    in_reps: list[str]
    in_nonreps: list[str] = field(default_factory=list)
    out_reps: list[str] = field(default_factory=list)
    out_nonreps: list[str] = field(default_factory=list)


def compute_group_frequency_table(
    matrix: PresenceMatrix, groups: GroupStrains, outgroup_min_bin_strains: int = 3,
    core_cutoff: float = 0.95, softcore_cutoff: float = 0.90, shell_cutoff: float = 0.15,
) -> list[dict]:
    """One row per family: ingroup frequency/count/class (`frequency`,
    `strain_count`, `bin`) and, when the outgroup has at least
    `outgroup_min_bin_strains` representatives, the outgroup ones
    (`frequency_out`, `strain_count_out`, `bin_out`; None/None/"-" otherwise)."""
    cut = (core_cutoff, softcore_cutoff, shell_cutoff)
    n_in, n_out = len(groups.in_reps), len(groups.out_reps)
    bin_out_on = n_out >= outgroup_min_bin_strains
    rows = []
    for fam in matrix.families:
        def present(strains):
            return [s for s in strains if matrix.is_present(fam, s)]
        c_in = len(present(groups.in_reps))
        c_out = len(present(groups.out_reps))
        in_nonrep = bool(present(groups.in_nonreps))
        out_nonrep = bool(present(groups.out_nonreps))
        in_any = c_in > 0 or in_nonrep
        out_any = c_out > 0 or out_nonrep
        row = {
            "family": fam,
            "frequency": c_in / n_in if n_in else 0.0,
            "strain_count": c_in,
            "bin": group_class(c_in, n_in, in_nonrep, out_any, True, *cut),
            "frequency_out": None, "strain_count_out": None, "bin_out": "-",
        }
        if bin_out_on:
            row["frequency_out"] = c_out / n_out
            row["strain_count_out"] = c_out
            row["bin_out"] = group_class(c_out, n_out, out_nonrep, in_any, False, *cut)
        rows.append(row)
    return rows


def _representatives(inventory_path: str | None) -> set[str] | None:
    """Representative Shorts, or None when dereplication is off (no path, or
    the 0-byte EMPTY_INVENTORY_STUB file)."""
    if not inventory_path or Path(inventory_path).stat().st_size == 0:
        return None
    return set(read_representative_shorts(inventory_path))


def resolve_group_strains(
    matrix: PresenceMatrix, config_path: str, inventory_path: str | None,
    ingroup_label: str, outgroup_label: str,
) -> GroupStrains:
    """Ingroup/outgroup strains from the samplesheet that are matrix columns,
    split into representatives and non-representatives (spec section 1)."""
    from config_parser import parse_config  # noqa: E402

    samples = parse_config(config_path)
    in_matrix = set(matrix.strains)
    reps = _representatives(inventory_path)
    missing = [s.short for s in samples
               if s.group in (ingroup_label, outgroup_label) and s.short not in in_matrix]
    if missing:
        print(f"WARNING: {len(missing)} ingroup/outgroup strains are absent from the matrix "
              f"columns and are ignored (first few: {missing[:5]})", file=sys.stderr)

    def split(label):
        members = [s.short for s in samples if s.group == label and s.short in in_matrix]
        if reps is None:
            return members, []
        return [s for s in members if s in reps], [s for s in members if s not in reps]

    in_reps, in_nonreps = split(ingroup_label)
    out_reps, out_nonreps = split(outgroup_label)
    return GroupStrains(in_reps, in_nonreps, out_reps, out_nonreps)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--config", required=True,
                    help="samplesheet CSV -- its ingroup strains are the frequency denominator")
    ap.add_argument("--ingroup_label", default="IN",
                    help="samplesheet GROUP value identifying ingroup strains (default: IN)")
    ap.add_argument("--outgroup_label", default="OUT",
                    help="samplesheet GROUP value identifying outgroup strains (default: OUT)")
    ap.add_argument("--outgroup_min_bin_strains", type=int, default=3,
                    help="bin the outgroup only with at least this many representatives (default: 3)")
    ap.add_argument("--inventory",
                    help="strain_inventory.tsv from pangenome_dereplicate_strains.py; when "
                         "given, only is_representative == 1 strains are counted")
    ap.add_argument("--core_cutoff", type=float, default=0.95)
    ap.add_argument("--softcore_cutoff", type=float, default=0.90)
    ap.add_argument("--shell_cutoff", type=float, default=0.15)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    groups = resolve_group_strains(matrix, args.config, args.inventory,
                                   args.ingroup_label, args.outgroup_label)
    if not groups.in_reps:
        print("ERROR: no ingroup strains left to bin over (check --config/--inventory "
              "against the matrix's columns)", file=sys.stderr)
        sys.exit(1)
    binned_out = len(groups.out_reps) >= args.outgroup_min_bin_strains
    print(f"Binning over {len(groups.in_reps)} ingroup representatives "
          f"({len(groups.in_nonreps)} non-representatives); outgroup: "
          f"{len(groups.out_reps)} representatives, "
          f"{'binned' if binned_out else 'not binned'}", file=sys.stderr)

    table = compute_group_frequency_table(
        matrix, groups, args.outgroup_min_bin_strains,
        args.core_cutoff, args.softcore_cutoff, args.shell_cutoff,
    )
    with open(args.output, "w") as fh:
        fh.write("family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out\n")
        for row in table:
            f_out = "-" if row["frequency_out"] is None else f"{row['frequency_out']:.4f}"
            c_out = "-" if row["strain_count_out"] is None else str(row["strain_count_out"])
            fh.write(f"{row['family']}\t{row['frequency']:.4f}\t{row['strain_count']}\t{row['bin']}\t"
                     f"{f_out}\t{c_out}\t{row['bin_out']}\n")

if __name__ == "__main__":
    main()
