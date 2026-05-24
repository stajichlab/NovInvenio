#!/usr/bin/env python3
"""Split a multi-FASTA file into one file per sequence."""

import argparse
import re
from pathlib import Path


def safe_filename(seq_id):
    return re.sub(r"[^\w.-]", "_", seq_id)


def split_fasta(input_path, prefix, suffix):
    current_id = None
    lines = []

    def flush(seq_id, seq_lines):
        fname = f"{prefix}{safe_filename(seq_id)}{suffix}"
        Path(fname).write_text("".join(seq_lines))

    with open(input_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    flush(current_id, lines)
                current_id = line[1:].split()[0]
                lines = [line]
            else:
                lines.append(line)
    if current_id is not None:
        flush(current_id, lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input multi-FASTA file")
    parser.add_argument("--prefix", default="", help="Output filename prefix")
    parser.add_argument("--suffix", default=".fa", help="Output filename suffix")
    args = parser.parse_args()

    split_fasta(args.input, args.prefix, args.suffix)


if __name__ == "__main__":
    main()
