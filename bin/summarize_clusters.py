#!/usr/bin/env python3
"""Summarise mmseqs2 cluster output into a per-cluster table."""
import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cluster-tsv', required=True,
                    help='mmseqs easy-cluster *_cluster.tsv output')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    clusters: dict[str, list] = defaultdict(list)
    with open(args.cluster_tsv) as fh:
        for line in fh:
            rep, member = line.strip().split('\t')
            clusters[rep].append(member)

    with open(args.output, 'w') as out:
        out.write('representative\tcluster_size\tmembers\n')
        for rep, members in sorted(clusters.items(), key=lambda x: -len(x[1])):
            out.write(f'{rep}\t{len(members)}\t{",".join(members)}\n')


if __name__ == '__main__':
    main()
