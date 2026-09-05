#!/usr/bin/env python3
"""Extend singletons.fa with each singleton's own within-genome paralog (issue #52).

The singleton pairwise search (workflows/novelty_discovery.nf's SINGLETON_*_SEARCH)
applies the same significance + paralog-competition filters bin/build_presence_matrix.py
uses (bin/novelty_presence_matrix.py). The paralog-competition check needs the
*paralog's* own hits against the same targets to compare against -- which only exist if
the paralog was searched too. This script just
extends the query FASTA; the actual filtering happens in
bin/novelty_presence_matrix.py, which only ever emits presence for true singletons
(from the mmseqs cluster_tsv), never for a paralog added here purely to support the
comparison.

Paralogs are looked up from paralog_cutoffs.tsv (protein_ID -> paralog_protein_ID), the
same self-vs-self search output bin/build_presence_matrix.py reads. A singleton with no
detected paralog contributes nothing extra. Paralog sequences are pulled from
--seed-fasta (the concatenated DISCOVERY_TARGET proteomes); a paralog already present in
singletons.fa (i.e. it's also a singleton) is not duplicated.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import read_fasta  # noqa: E402


def load_paralog_of(cutoff_files):
    paralog_of = {}
    for path in cutoff_files:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        paralog_of.update(dict(zip(df['protein_ID'], df['paralog_protein_ID'])))
    return paralog_of


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--singletons-fa', required=True, dest='singletons_fa',
                    help='singletons.fa from bin/extract_singletons.py')
    ap.add_argument('--seed-fasta', required=True, dest='seed_fasta',
                    help='Concatenated DISCOVERY_TARGET proteome FASTA (seed_all.faa), '
                         'to pull paralog sequences not already in --singletons-fa')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[], dest='paralog_cutoffs',
                    help='paralog_cutoffs.tsv files from the DISCOVERY_TARGET self-vs-self search')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    singletons = read_fasta(args.singletons_fa)
    paralog_of = load_paralog_of(args.paralog_cutoffs)

    wanted_paralogs = {
        paralog_of[sid] for sid in singletons
        if sid in paralog_of and pd.notna(paralog_of[sid]) and paralog_of[sid] not in singletons
    }

    extra = {}
    if wanted_paralogs:
        seed = read_fasta(args.seed_fasta)
        extra = {pid: seed[pid] for pid in wanted_paralogs if pid in seed}

    with open(args.output, 'w') as out:
        for sid, rec in singletons.items():
            out.write(f'>{sid}\n{rec.seq}\n')
        for pid, rec in extra.items():
            out.write(f'>{pid}\n{rec.seq}\n')

    print(f"Extended singleton query: {len(singletons)} singletons + "
          f"{len(extra)} paralogs (of {len(wanted_paralogs)} wanted)", file=sys.stderr)


if __name__ == '__main__':
    main()
