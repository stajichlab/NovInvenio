#!/usr/bin/env python3
"""Split a multi-FASTA file into fixed-size chunks (N sequences per chunk file).

Unlike split_fasta.py (one file per sequence -- used for per-family FASTAs where a
downstream step needs to address a single sequence by name), this is for scatter-gather
parallelism over a large flat sequence set: fixed-size batches so a fan-out step has a
bounded, predictable number of tasks regardless of the input's sequence count.
"""
import argparse
from pathlib import Path


def split_fasta_chunks(input_path, chunk_size, outdir, prefix):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    idx = 0
    n_in_chunk = 0
    out_fh = None

    def open_chunk(i):
        return open(outdir / f"{prefix}{i:06d}.fa", 'w')

    with open(input_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if n_in_chunk >= chunk_size:
                    out_fh.close()
                    idx += 1
                    n_in_chunk = 0
                    out_fh = None
                if out_fh is None:
                    out_fh = open_chunk(idx)
                n_in_chunk += 1
            out_fh.write(line)

    if out_fh is not None:
        out_fh.close()
    elif idx == 0:
        # Empty input -- still emit one (empty) chunk so downstream fan-out has something
        # to consume, matching SPLIT_HMM_DB's empty-input handling.
        open_chunk(0).close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True, help='Input multi-FASTA file')
    ap.add_argument('--chunk-size', type=int, required=True, dest='chunk_size',
                    help='Sequences per chunk file')
    ap.add_argument('--outdir', default='.', help='Output directory')
    ap.add_argument('--prefix', default='chunk_', help='Output filename prefix')
    args = ap.parse_args()

    split_fasta_chunks(args.input, args.chunk_size, args.outdir, args.prefix)


if __name__ == '__main__':
    main()
