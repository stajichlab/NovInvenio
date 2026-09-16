#!/usr/bin/env python3
"""Short-prefix every header of a FASTA file for the pangenome-profiling
subworkflow's tier-1 clustering step.

Ported from a manual awk one-liner documented in NovInvenio_Investigations'
Afumigatus_pangenome study (see NEXTFLOW_MIGRATION_NOTES.md section A item 1
and pangenome_build_presence_matrix.py's "ID conventions" docstring section)
-- the study ran this by hand per-strain before concatenating; this gives
the pangenome-profiling subworkflow a real, tested PREFIX_PROTEOMES /
PREFIX_GENOME process body instead of a shell one-liner with hardcoded
column indices.

Rewrites every header ``>original_id [description...]`` to
``><Short><id_sep><original_id> [description...]`` -- only the first
whitespace-delimited token is touched, so any trailing description text is
preserved unchanged. `--id_sep` must match whatever every other
`pangenome_*` script's `--id_sep` is set to (default ``|``); a family ID
downstream is only recoverable back to its strain if this prefix is applied
consistently to every strain's proteome AND to every strain's genome FASTA
before the rescue-pass BLAST database is built from it.

Usage:
  pangenome_prefix_fasta.py --short AF293 --input AF293.pep.fa \\
      --output AF293.prefixed.pep.fa
"""
from __future__ import annotations

import argparse
import sys


def prefix_fasta(input_path: str, output_path: str, short: str, id_sep: str = "|") -> int:
    """Rewrite every FASTA header's first token as ``><short><id_sep><id>``.
    Returns the number of records rewritten."""
    n = 0
    with open(input_path) as fh, open(output_path, "w") as out:
        for line in fh:
            if line.startswith(">"):
                rest = line[1:].rstrip("\n")
                first, sep, tail = rest.partition(" ")
                out.write(f">{short}{id_sep}{first}")
                if sep:
                    out.write(f" {tail}")
                out.write("\n")
                n += 1
            else:
                out.write(line)
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--short", required=True, help="Short strain ID (config.csv 'Short' column)")
    ap.add_argument("--input", required=True, help="FASTA to prefix (protein or genome/DNA)")
    ap.add_argument("--id_sep", default="|", help="prefix separator (default: '|')")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    n = prefix_fasta(args.input, args.output, args.short, args.id_sep)
    print(f"pangenome_prefix_fasta: {n} records prefixed with '{args.short}{args.id_sep}'", file=sys.stderr)


if __name__ == "__main__":
    main()
