#!/usr/bin/env python3
"""Per-strain query extraction for the genome-level tblastn rescue pass.

For each strain, writes a FASTA of only the tier-1 family representative
sequences currently ABSENT in that strain (per `--matrix`) -- rescue only
ever needs to check cells the protein-level clustering already called
absent, so this also cuts total tblastn work relative to searching the
full family set against every strain.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/extract_absent_family_queries.py), which itself replaced an earlier
combined-genome-database rescue design after a real correctness bug: a
shared ``-max_target_seqs`` cap was found to truncate hits ACROSS strains
per query, not per strain (89.4% of hit queries in that run landed at
exactly the cap). The per-strain scatter this feeds (one BLAST database per
strain, restricted to that strain's own absent-family query set) is the
structural fix -- see workflows/pangenome_profile.nf's RESCUE_PASS section.

Usage:
  pangenome_extract_absent_family_queries.py --matrix presence_matrix.tsv \\
      --rep_fasta tier1_rep_seq.fasta --out_dir per_strain_queries/
  # writes per_strain_queries/<Short>.absent.fa for every strain with at
  # least one absent family, plus per_strain_queries/manifest.tsv
  # (Short<TAB>n_absent_families<TAB>query_fasta_path)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import ABSENT, PresenceMatrix  # noqa: E402


def read_fasta_records(path: str) -> dict[str, list[str]]:
    """Parse a multi-FASTA file into {record_id: [header_line, seq_lines...]},
    keyed by the first whitespace-delimited token after '>' -- matches the
    "family ID = tier-1 cluster representative ID, verbatim" convention, so
    a presence-matrix family name looks up its representative sequence
    directly."""
    records: dict[str, list[str]] = {}
    current_id: str | None = None
    current_lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = current_lines
                current_id = line[1:].split(None, 1)[0].rstrip("\n")
                current_lines = [line]
            else:
                current_lines.append(line)
        if current_id is not None:
            records[current_id] = current_lines
    return records


def absent_families_by_strain(matrix: PresenceMatrix) -> dict[str, list[str]]:
    """{strain: [family, ...]} for every (family, strain) cell called
    ABSENT -- the per-strain rescue query set."""
    by_strain: dict[str, list[str]] = {s: [] for s in matrix.strains}
    for family in matrix.families:
        for strain in matrix.strains:
            if matrix.call(family, strain) == ABSENT:
                by_strain[strain].append(family)
    return by_strain


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--rep_fasta", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    matrix = PresenceMatrix.from_tsv(args.matrix)
    records = read_fasta_records(args.rep_fasta)
    by_strain = absent_families_by_strain(matrix)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    missing_reps: set[str] = set()
    for strain, families in sorted(by_strain.items()):
        if not families:
            continue
        query_path = out_dir / f"{strain}.absent.fa"
        with open(query_path, "w") as fh:
            for family in families:
                lines = records.get(family)
                if lines is None:
                    missing_reps.add(family)
                    continue
                fh.writelines(lines)
        manifest_rows.append((strain, len(families), str(query_path)))

    manifest_path = out_dir / "manifest.tsv"
    with open(manifest_path, "w") as fh:
        fh.write("strain\tn_absent_families\tquery_fasta\n")
        for strain, n, path in manifest_rows:
            fh.write(f"{strain}\t{n}\t{path}\n")

    if missing_reps:
        print(
            f"WARNING: {len(missing_reps)} family IDs in the matrix had no "
            f"matching representative sequence in {args.rep_fasta} -- "
            "skipped for every strain that needed them (mismatched "
            "matrix/rep_fasta inputs?)", file=sys.stderr,
        )
    print(
        f"Wrote {len(manifest_rows)} per-strain query FASTAs "
        f"({sum(n for _, n, _ in manifest_rows)} total absent-family "
        f"entries) to {out_dir}, manifest at {manifest_path}", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
