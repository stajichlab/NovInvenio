"""View B1 (issue #182): per-strain neighbourhood statistics for trans modules.

Synthetic genomes with known structure -- the spec's validation requirement
(docs/superpowers/specs/2026-09-19-pangenome-gainloss-visualization-design.md,
"Phase 2 -- View B").
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
from pangenome_neighborhood import pair_counts, rank_genes, score_modules  # noqa: E402

KB = 1000


def _genome(strain, contigs, spacing=KB):
    """One gene every `spacing` bp; each gene is its own background family
    unless relabelled. contigs: {contig: n_genes}."""
    rows = []
    for contig, n in contigs.items():
        for i in range(n):
            rows.append({'Short': strain, 'contig': contig, 'start': i * spacing + 1,
                         'end': i * spacing + 500, 'family': f'bg_{strain}_{contig}_{i}'})
    return pd.DataFrame(rows)


def _set_family(genes, strain, contig, index, family):
    mask = ((genes['Short'] == strain) & (genes['contig'] == contig)
            & (genes['start'] == index * KB + 1))
    assert mask.sum() == 1
    genes.loc[mask, 'family'] = family


def _members(genes, families):
    g = rank_genes(genes)
    return g[g['family'].isin(families)]


# ---- pair_counts: the per-strain pair classification ---------------------------

def test_pairs_within_distance_and_outside_window_are_colocalized():
    g = _genome('S', {'c1': 100})
    for i, f in [(0, 'A'), (20, 'B'), (40, 'C')]:
        _set_family(g, 'S', 'c1', i, f)
    c = pair_counts(_members(g, {'A', 'B', 'C'}), max_bp=100 * KB, min_gene_gap=11)
    assert c == {'total': 3, 'cross_contig': 0, 'adjacent': 0, 'informative': 3,
                 'colocalized': 3}


def test_pairs_inside_the_trans_window_are_adjacent_not_informative():
    # Rank gap 5 < 11: the trans definition already excludes these, so they must
    # not count either way (constraint 1: no circularity).
    g = _genome('S', {'c1': 100})
    _set_family(g, 'S', 'c1', 0, 'A')
    _set_family(g, 'S', 'c1', 5, 'B')
    c = pair_counts(_members(g, {'A', 'B'}), max_bp=100 * KB, min_gene_gap=11)
    assert (c['adjacent'], c['informative'], c['colocalized']) == (1, 0, 0)


def test_gap_exactly_min_gene_gap_is_informative():
    g = _genome('S', {'c1': 100})
    _set_family(g, 'S', 'c1', 0, 'A')
    _set_family(g, 'S', 'c1', 11, 'B')
    c = pair_counts(_members(g, {'A', 'B'}), max_bp=100 * KB, min_gene_gap=11)
    assert (c['adjacent'], c['informative'], c['colocalized']) == (0, 1, 1)


def test_cross_contig_pairs_are_excluded_not_distant():
    # Constraint 2: different contigs may be neighbours in the real chromosome.
    g = _genome('S', {'c1': 50, 'c2': 50})
    _set_family(g, 'S', 'c1', 0, 'A')
    _set_family(g, 'S', 'c2', 0, 'B')
    c = pair_counts(_members(g, {'A', 'B'}), max_bp=100 * KB, min_gene_gap=11)
    assert (c['total'], c['cross_contig'], c['informative']) == (1, 1, 0)


def test_same_contig_beyond_max_distance_is_informative_but_not_colocalized():
    g = _genome('S', {'c1': 300})
    _set_family(g, 'S', 'c1', 0, 'A')
    _set_family(g, 'S', 'c1', 200, 'B')        # 200 kb apart
    c = pair_counts(_members(g, {'A', 'B'}), max_bp=100 * KB, min_gene_gap=11)
    assert (c['informative'], c['colocalized']) == (1, 0)


def test_pairs_of_the_same_family_are_ignored():
    # Two copies of one family are paralogs, not a module relationship.
    g = _genome('S', {'c1': 100})
    _set_family(g, 'S', 'c1', 0, 'A')
    _set_family(g, 'S', 'c1', 20, 'A')
    _set_family(g, 'S', 'c1', 40, 'B')
    c = pair_counts(_members(g, {'A', 'B'}), max_bp=100 * KB, min_gene_gap=11)
    assert c['total'] == 2
    assert c['colocalized'] == 2


def test_rank_counts_every_gene_on_the_contig_not_only_members():
    # Ranks come from the whole gene complement, so a gap of 20 genes is 20 even
    # though only two genes are module members.
    g = _genome('S', {'c1': 30})
    ranked = rank_genes(g)
    assert list(ranked.sort_values('start')['rank'])[:3] == [0, 1, 2]
    assert ranked['rank'].max() == 29


# ---- score_modules: aggregation and the permutation null -------------------------

def _clustered_and_random(seed=1, n_strains=6):
    """Module 'clus': 6 families placed 12 genes apart inside one 70 kb block of
    each strain. Module 'rand': 6 families at random positions. 8 contigs of 750
    genes (750 kb) per strain, so 100 kb is local, not contig-wide."""
    rng = np.random.default_rng(seed)
    frames = []
    for s in range(n_strains):
        strain = f'S{s}'
        g = _genome(strain, {f'c{k}': 750 for k in range(8)})
        block_contig = f'c{rng.integers(8)}'
        start = int(rng.integers(0, 750 - 72))
        for j in range(6):
            _set_family(g, strain, block_contig, start + 12 * j, f'clus{j}')
        free = g.index[~g['family'].str.startswith('clus')]
        for j, idx in enumerate(rng.choice(free, size=6, replace=False)):
            g.loc[idx, 'family'] = f'rand{j}'
        frames.append(g)
    genes = pd.concat(frames, ignore_index=True)
    module_of = {**{f'clus{j}': 'clus' for j in range(6)},
                 **{f'rand{j}': 'rand' for j in range(6)}}
    return genes, module_of


def test_null_flags_a_clustered_module_and_not_a_random_one():
    genes, module_of = _clustered_and_random()
    res = {r['module_id']: r for r in score_modules(genes, module_of, n_perm=200, seed=0,
                                                    null_pool='all')}
    clus, rand = res['clus'], res['rand']
    assert clus['obs_frac'] == pytest.approx(1.0)
    assert clus['p_empirical'] <= 0.01
    # frac_informative cannot exceed 1, so its effect ratio is capped at
    # 1/null (~4.7 on 750 kb contigs); frac_total has no such ceiling here.
    assert clus['effect_ratio'] > 3
    assert clus['effect_ratio_total'] > 10
    assert rand['p_empirical'] > 0.05
    # The all-pairs statistic sees the same signal.
    assert clus['p_empirical_total'] <= 0.01
    assert rand['p_empirical_total'] > 0.05


def test_score_reports_scale_pool_and_exclusion_accounting():
    genes, module_of = _clustered_and_random()
    r = {x['module_id']: x for x in score_modules(genes, module_of, n_perm=50, seed=0,
                                                  null_pool='all', max_kb=100,
                                                  min_gene_gap=11)}['rand']
    assert (r['max_kb'], r['min_gene_gap'], r['null_pool'], r['n_perm'], r['seed']) == \
        (100, 11, 'all', 50, 0)
    assert r['strains_scored'] == 6
    assert r['pairs_total'] == 6 * 15
    assert r['pairs_total'] == r['pairs_cross_contig'] + r['pairs_adjacent'] + r['pairs_informative']
    assert r['cross_contig_frac'] == pytest.approx(r['pairs_cross_contig'] / r['pairs_total'])


def test_score_is_deterministic_for_a_seed():
    genes, module_of = _clustered_and_random()
    a = score_modules(genes, module_of, n_perm=30, seed=7, null_pool='all')
    b = score_modules(genes, module_of, n_perm=30, seed=7, null_pool='all')
    assert a == b


def test_strain_with_fewer_than_two_member_genes_is_not_scored():
    genes, module_of = _clustered_and_random(n_strains=3)
    genes = genes[~((genes['Short'] == 'S0') & genes['family'].str.startswith('rand')
                    & (genes['family'] != 'rand0'))]
    r = {x['module_id']: x for x in score_modules(genes, module_of, n_perm=10, seed=0,
                                                  null_pool='all')}['rand']
    assert r['strains_scored'] == 2


def _accessory_region(n_strains=5):
    """All accessory genes (module members included) sit on one 150 kb contig
    'acc'; the 8 core contigs are 750 kb. Members are 12 genes apart on 'acc'."""
    frames, accessory = [], set()
    for s in range(n_strains):
        strain = f'S{s}'
        g = _genome(strain, {**{f'c{k}': 750 for k in range(8)}, 'acc': 150})
        for j in range(6):
            _set_family(g, strain, 'acc', 10 + 12 * j, f'm{j}')
        accessory |= set(g.loc[g['contig'] == 'acc', 'family'])
        frames.append(g)
    return pd.concat(frames, ignore_index=True), {f'm{j}': 'M' for j in range(6)}, accessory


def test_accessory_pool_samples_only_accessory_genes():
    # An all-genes null mistakes "accessory genes share a region" for module
    # clustering; the accessory-only null samples from that region too.
    genes, module_of, accessory = _accessory_region()
    all_null = score_modules(genes, module_of, n_perm=200, seed=0, null_pool='all')[0]
    acc_null = score_modules(genes, module_of, n_perm=200, seed=0, null_pool='accessory',
                             accessory_families=accessory)[0]
    assert all_null['p_empirical'] <= 0.01
    assert acc_null['null_mean_frac'] > 3 * all_null['null_mean_frac']
    assert acc_null['effect_ratio'] < all_null['effect_ratio']
    assert acc_null['null_pool'] == 'accessory'


def test_accessory_pool_requires_accessory_families():
    genes, module_of, _ = _accessory_region(n_strains=2)
    with pytest.raises(ValueError):
        score_modules(genes, module_of, n_perm=5, seed=0, null_pool='accessory')


def test_empty_modules_give_no_rows():
    genes, _ = _clustered_and_random(n_strains=2)
    assert score_modules(genes, {}, n_perm=5, seed=0, null_pool='all') == []
