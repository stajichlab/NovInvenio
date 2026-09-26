"""View B1 (issue #182): per-strain neighbourhood statistics for Leiden trans
modules. Design: docs/superpowers/specs/2026-09-19-pangenome-gainloss-
visualization-design.md, "View B -- B1".

For each module and each strain carrying >= 2 member genes, count pairs of
member genes (from different families) in that strain's own assembly:

  total        every member-gene pair from different families
  cross_contig the two genes are on different contigs. Excluded as
               UNINFORMATIVE, never counted as "distant": with a median of 283
               scaffolds per assembly, two genes on different contigs may be
               neighbours in the real chromosome (spec constraint 2).
  adjacent     same contig, gene-rank gap < min_gene_gap. Excluded: the `trans`
               class is defined by low linkage inside a k=10 gene window, so a
               statistic at that scale would be circular (spec constraint 1).
  informative  same contig, gene-rank gap >= min_gene_gap.
  colocalized  informative, and gene starts within max_bp of each other.

Positions are never compared across strains. Two module statistics, each with
its own null:

  frac_informative = sum(colocalized) / sum(informative)   (the spec's statistic)
  frac_total       = sum(colocalized) / sum(total)

frac_informative conditions on "same contig", so it ignores the evidence that
two genes share a contig at all. On a fragmented assembly, where many contigs
are about as long as max_kb, almost every same-contig pair is "within max_kb"
and the statistic loses power. frac_total keeps that evidence. It is not a
fragmentation measure, because the null draws from the same strain and so
carries the same contig breaks; a contig break costs it power, never a false
signal. Which one to lead with is decided by the empirical pass (issue #182). The permutation null
redraws, in every scored strain, the same number of genes from that strain's
own pool and recomputes the same fraction; p = (1 + #null >= obs) / (1 + n_perm).

The null pool matters. `all` draws from every gene in the strain. `accessory`
draws only from genes whose family is accessory (not core/soft-core). Accessory
genes cluster in general (islands, subtelomeres), so an all-genes null can call
any accessory module "clustered"; the accessory null asks whether the module
clusters MORE than accessory genes already do. Every output row records the
pool, the scale (max_kb, min_gene_gap), n_perm and the seed.

Distances are gene START to gene START; ranks are each gene's 0-based order by
start among ALL genes on its contig in that strain (the same rank
bin/pangenome_build_family_positions.py produces), not among module members.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

NULL_POOLS = ('accessory', 'all')


def rank_genes(genes: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of `genes` (Short, contig, start, ...) with a 0-based
    `rank`: order by start among all genes on the same contig of the same strain."""
    out = genes.sort_values(['Short', 'contig', 'start'], kind='mergesort').copy()
    out['rank'] = out.groupby(['Short', 'contig'], sort=False).cumcount()
    return out


def _pairs(counts: np.ndarray) -> int:
    counts = counts[counts > 1].astype(np.int64)
    return int((counts * (counts - 1) // 2).sum())


def _window_pairs(group: np.ndarray, rank: np.ndarray, start: np.ndarray,
                  min_gap: int, max_bp: int, rank_span: int, bp_span: int):
    """Pairs (i < j) in the same group with rank gap >= min_gap, and those also
    within max_bp. Inputs must be sorted by (group, rank). The group sits in the
    high part of each search key, so a rank or bp window cannot reach into the
    next group, and the group's end index bounds the rank window."""
    key_rank = group * rank_span + rank
    key_bp = group * bp_span + start
    lo = np.searchsorted(key_rank, key_rank + min_gap, side='left')
    hi = np.searchsorted(key_bp, key_bp + max_bp, side='right')
    end = np.searchsorted(key_rank, (group + 1) * rank_span, side='left')
    informative = int(np.maximum(0, end - lo).sum())
    colocalized = int(np.maximum(0, np.minimum(hi, end) - lo).sum())
    return informative, colocalized


def _counts(contig: np.ndarray, rank: np.ndarray, start: np.ndarray, fam: np.ndarray,
            min_gap: int, max_bp: int) -> dict:
    """Pair counts for one strain's gene set. All inputs are int arrays."""
    n = len(rank)
    if n < 2:
        return {'total': 0, 'cross_contig': 0, 'adjacent': 0, 'informative': 0,
                'colocalized': 0}
    rank_span = int(rank.max()) + min_gap + 1
    bp_span = int(start.max()) + max_bp + 1
    n_contig = int(contig.max()) + 1

    order = np.lexsort((rank, contig))
    inf_all, col_all = _window_pairs(contig[order], rank[order], start[order],
                                     min_gap, max_bp, rank_span, bp_span)

    # Same-family pairs are paralogs, not a module relationship: subtract them.
    fc = fam.astype(np.int64) * n_contig + contig
    order_f = np.lexsort((rank, fc))
    inf_same, col_same = _window_pairs(fc[order_f], rank[order_f], start[order_f],
                                       min_gap, max_bp, rank_span, bp_span)

    _, fam_counts = np.unique(fam, return_counts=True)
    _, contig_counts = np.unique(contig, return_counts=True)
    _, fc_counts = np.unique(fc, return_counts=True)
    total = n * (n - 1) // 2 - _pairs(fam_counts)
    same_contig = _pairs(contig_counts) - _pairs(fc_counts)
    informative = inf_all - inf_same
    return {
        'total': int(total),
        'cross_contig': int(total - same_contig),
        'adjacent': int(same_contig - informative),
        'informative': int(informative),
        'colocalized': int(col_all - col_same),
    }


def pair_counts(members: pd.DataFrame, *, max_bp: int, min_gene_gap: int) -> dict:
    """Pair counts for one strain's module members (a rank_genes() slice)."""
    contig = pd.factorize(members['contig'])[0].astype(np.int64)
    fam = pd.factorize(members['family'])[0].astype(np.int64)
    return _counts(contig, members['rank'].to_numpy(np.int64),
                   members['start'].to_numpy(np.int64), fam, min_gene_gap, max_bp)


def _statistic(obs_num: int, obs_den: int, null_num: np.ndarray, null_den: np.ndarray) -> dict:
    """Observed fraction, null mean/sd, effect ratio and empirical p for one
    num/den statistic. Permutations with a zero denominator are dropped."""
    obs_frac = obs_num / obs_den if obs_den else math.nan
    keep = null_den > 0
    null = null_num[keep] / null_den[keep]
    null_mean = float(null.mean()) if len(null) else math.nan
    null_sd = float(null.std(ddof=1)) if len(null) > 1 else math.nan
    if math.isnan(obs_frac) or not len(null):
        p_emp = math.nan
    else:
        p_emp = (1 + int((null >= obs_frac).sum())) / (1 + len(null))
    effect = obs_frac / null_mean if null_mean and not math.isnan(obs_frac) else math.nan
    return {'obs_frac': obs_frac, 'null_mean_frac': null_mean, 'null_sd_frac': null_sd,
            'effect_ratio': effect, 'p_empirical': p_emp}


def score_modules(genes: pd.DataFrame, module_of: dict, *, max_kb: float = 100,
                  min_gene_gap: int = 11, n_perm: int = 200, seed: int = 0,
                  null_pool: str = 'accessory', accessory_families=None) -> list[dict]:
    """One result dict per module with >= 1 scored strain.

    genes: one row per gene with Short, contig, start, family (every gene of every
    strain, not only module members -- ranks and the null need the full complement).
    module_of: family -> module_id. accessory_families: required for
    null_pool='accessory'."""
    if null_pool not in NULL_POOLS:
        raise ValueError(f'null_pool must be one of {NULL_POOLS}, not {null_pool!r}')
    if null_pool == 'accessory' and accessory_families is None:
        raise ValueError("null_pool='accessory' needs accessory_families")
    if not module_of:
        return []
    max_bp = int(round(max_kb * 1000))
    ranked = rank_genes(genes)
    ranked['_fam'] = pd.factorize(ranked['family'])[0].astype(np.int64)
    ranked['_mod'] = ranked['family'].map(module_of)
    if null_pool == 'accessory':
        ranked['_pool'] = ranked['family'].isin(set(accessory_families))
    else:
        ranked['_pool'] = True

    modules = sorted(ranked['_mod'].dropna().unique(), key=str)
    rng = np.random.default_rng(seed)
    strains = []
    for short, g in ranked.groupby('Short', sort=True):
        contig = pd.factorize(g['contig'])[0].astype(np.int64)
        strains.append({
            'contig': contig, 'rank': g['rank'].to_numpy(np.int64),
            'start': g['start'].to_numpy(np.int64), 'fam': g['_fam'].to_numpy(),
            'mod': g['_mod'].to_numpy(), 'pool': np.flatnonzero(g['_pool'].to_numpy()),
        })

    rows = []
    for m in modules:
        obs = dict.fromkeys(('total', 'cross_contig', 'adjacent', 'informative',
                             'colocalized'), 0)
        null_inf = np.zeros(n_perm, dtype=np.int64)
        null_col = np.zeros(n_perm, dtype=np.int64)
        null_tot = np.zeros(n_perm, dtype=np.int64)
        n_strains = n_genes = 0
        for s in strains:
            idx = np.flatnonzero(s['mod'] == m)
            if len(idx) < 2:
                continue
            n_strains += 1
            n_genes += len(idx)
            c = _counts(s['contig'][idx], s['rank'][idx], s['start'][idx], s['fam'][idx],
                        min_gene_gap, max_bp)
            for k in obs:
                obs[k] += c[k]
            # Members always belong to their own strain's pool, so a small pool
            # can still supply a draw of the same size.
            pool = np.union1d(s['pool'], idx)
            for p in range(n_perm):
                draw = rng.choice(pool, size=len(idx), replace=False)
                nc = _counts(s['contig'][draw], s['rank'][draw], s['start'][draw],
                             s['fam'][draw], min_gene_gap, max_bp)
                null_inf[p] += nc['informative']
                null_col[p] += nc['colocalized']
                null_tot[p] += nc['total']
        if not n_strains:
            continue
        inf_stat = _statistic(obs['colocalized'], obs['informative'], null_col, null_inf)
        tot_stat = _statistic(obs['colocalized'], obs['total'], null_col, null_tot)
        rows.append({
            'module_id': m,
            'strains_scored': n_strains,
            'member_genes': n_genes,
            'pairs_total': obs['total'],
            'pairs_cross_contig': obs['cross_contig'],
            'pairs_adjacent': obs['adjacent'],
            'pairs_informative': obs['informative'],
            'pairs_colocalized': obs['colocalized'],
            'cross_contig_frac': obs['cross_contig'] / obs['total'] if obs['total'] else math.nan,
            **{f'{k}': v for k, v in inf_stat.items()},
            **{f'{k}_total': v for k, v in tot_stat.items()},
            'null_pool': null_pool,
            'max_kb': max_kb,
            'min_gene_gap': min_gene_gap,
            'n_perm': n_perm,
            'seed': seed,
        })
    return rows
