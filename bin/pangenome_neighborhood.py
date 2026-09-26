#!/usr/bin/env python3
"""View B1 (issue #182): per-strain neighbourhood statistics for Leiden trans
modules -- how close each module's member genes sit inside each strain's own
assembly, against a within-strain permutation null. The statistics, the
exclusion classes and the null pools are documented in
lib/pangenome_neighborhood.py.

Inputs are pipeline outputs of the same run (family IDs are not stable across
runs, so never mix files from two runs):
  --gene_positions   gene_positions.tsv[.zst] (Short, protein_id, contig, start, end)
  --cluster_tsv      tier1_cluster.tsv (family representative, member "Short|protein_id")
  --family_modules   family_modules.tsv (family, module_id, module_size)
  --frequency_table  frequency_table.tsv (family, ..., bin); needed for --null_pool accessory

Only annotated genes (with a protein in gene_positions) are placed; rescued
genome-only calls have no gene model and are not scored.

Usage:
  pangenome_neighborhood.py --gene_positions gene_positions.tsv.zst \\
      --cluster_tsv tier1_cluster.tsv --family_modules family_modules.tsv \\
      --frequency_table frequency_table.tsv --max_kb 100 --min_gene_gap 11 \\
      --n_perm 200 --seed 0 --null_pool accessory --output module_neighborhood.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from compressed_io import open_maybe_compressed  # noqa: E402
from pangenome_neighborhood import NULL_POOLS, score_modules  # noqa: E402

CORE_BINS = {'core', 'soft_core'}
COLUMNS = ['module_id', 'module_size', 'strains_scored', 'member_genes',
           'pairs_total', 'pairs_cross_contig', 'pairs_adjacent', 'pairs_informative',
           'pairs_colocalized', 'cross_contig_frac',
           'obs_frac', 'null_mean_frac', 'null_sd_frac', 'effect_ratio', 'p_empirical',
           'obs_frac_total', 'null_mean_frac_total', 'null_sd_frac_total',
           'effect_ratio_total', 'p_empirical_total',
           'null_pool', 'max_kb', 'min_gene_gap', 'n_perm', 'seed']


def read_tsv(path, **kw) -> pd.DataFrame:
    with open_maybe_compressed(path) as fh:
        return pd.read_csv(fh, sep='\t', **kw)


def load_genes(gene_positions, cluster_tsv) -> pd.DataFrame:
    genes = read_tsv(gene_positions, dtype={'Short': str, 'protein_id': str, 'contig': str})
    clusters = read_tsv(cluster_tsv, header=None, names=['family', 'member'], dtype=str)
    fam_of = dict(zip(clusters['member'], clusters['family']))
    genes['family'] = (genes['Short'] + '|' + genes['protein_id']).map(fam_of)
    n_missing = int(genes['family'].isna().sum())
    if n_missing:
        print(f'WARNING: {n_missing} of {len(genes)} genes have no tier-1 family and '
              'are ignored.', file=sys.stderr)
    return genes.dropna(subset=['family'])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gene_positions', required=True)
    ap.add_argument('--cluster_tsv', required=True)
    ap.add_argument('--family_modules', required=True)
    ap.add_argument('--frequency_table', default=None)
    ap.add_argument('--max_kb', type=float, default=100.0)
    ap.add_argument('--min_gene_gap', type=int, default=11)
    ap.add_argument('--n_perm', type=int, default=200)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--null_pool', choices=NULL_POOLS, default='accessory')
    ap.add_argument('--min_module_size', type=int, default=2,
                    help='Skip modules with fewer member families (default 2).')
    ap.add_argument('--output', required=True)
    args = ap.parse_args(argv)

    modules = read_tsv(args.family_modules, dtype={'family': str, 'module_id': str})
    sizes = modules.groupby('module_id')['family'].nunique()
    keep = set(sizes[sizes >= args.min_module_size].index)
    modules = modules[modules['module_id'].isin(keep)]
    module_of = dict(zip(modules['family'], modules['module_id']))

    rows = []
    if module_of:
        accessory = None
        if args.null_pool == 'accessory':
            if not args.frequency_table:
                sys.exit('ERROR: --null_pool accessory needs --frequency_table.')
            freq = read_tsv(args.frequency_table, dtype={'family': str})
            accessory = set(freq.loc[~freq['bin'].isin(CORE_BINS), 'family'])
        genes = load_genes(args.gene_positions, args.cluster_tsv)
        rows = score_modules(genes, module_of, max_kb=args.max_kb,
                             min_gene_gap=args.min_gene_gap, n_perm=args.n_perm,
                             seed=args.seed, null_pool=args.null_pool,
                             accessory_families=accessory)
        for r in rows:
            r['module_size'] = int(sizes[r['module_id']])
    out = pd.DataFrame(rows, columns=COLUMNS)
    out.to_csv(args.output, sep='\t', index=False, float_format='%.6g')
    print(f'Scored {len(out)} module(s); scale max_kb={args.max_kb:g}, '
          f'min_gene_gap={args.min_gene_gap}, null_pool={args.null_pool}, '
          f'n_perm={args.n_perm}, seed={args.seed}', file=sys.stderr)


if __name__ == '__main__':
    main()
