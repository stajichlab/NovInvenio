#!/usr/bin/env python3
"""Extract singleton protein sequences from an mmseqs cluster TSV.

Singletons are clusters with exactly one member (rep == member).  In the
novelty_discovery workflow, multi-member families are profiled with famsa
+ hmmbuild, while singletons are handled via pairwise search against the
DISCOVERY_OUT proteomes.  This script pulls the singleton sequences from the
concatenated target proteome FASTA and writes them as a single FASTA file.

Output: one FASTA file containing all singleton protein sequences.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import read_fasta  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep_id<TAB>member_id)')
    ap.add_argument('--fasta', required=True,
                    help='Concatenated FASTA the clustering was run on')
    ap.add_argument('--output', required=True,
                    help='Output FASTA file of singleton sequences')
    args = ap.parse_args()

    records = read_fasta(args.fasta)

    singletons = []
    with open(args.cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep == member:
                if member not in records:
                    sys.exit(f"ERROR: singleton '{member}' not found in {args.fasta}")
                singletons.append(member)

    with open(args.output, 'w') as out:
        for sid in singletons:
            rec = records[sid]
            out.write(f'>{sid}\n{rec.seq}\n')

    print(f"Wrote {len(singletons)} singleton sequences to {args.output}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
