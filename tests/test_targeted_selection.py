import pytest
from master_pool import MasterSample
from targeted_selection import (
    candidate_pool,
    exclude_species,
    rank_candidates,
    select_nearest,
    select_trait,
)


def _sample(species, phylum, klass, order, family, genus):
    return MasterSample(
        species=species, strain='', protein_path='', dna_path='',
        lineage=[phylum, '', klass, '', order, family, genus], ncbi_taxid='',
    )


FOCAL = _sample('Mucor circinelloides', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
POOL = [
    FOCAL,
    _sample('Rhizopus arrhizus', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Rhizopodaceae', 'Rhizopus'),
    _sample('Lichtheimia corymbifera', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Lichtheimiaceae', 'Lichtheimia'),
    _sample('Basidiobolus meristosporus', 'Basidiobolomycota', 'Basidiobolomycetes', 'Basidiobolales', 'Basidiobolaceae', 'Basidiobolus'),
]


def test_candidate_pool_excludes_focal_and_out_of_scope_species():
    candidates = candidate_pool(FOCAL.species, FOCAL.lineage, POOL, scope_rank='ORDER')
    assert {c.species for c in candidates} == {'Rhizopus arrhizus', 'Lichtheimia corymbifera'}


def test_rank_candidates_orders_deepest_match_first_then_alphabetically():
    same_family = _sample('Mucor mucedo', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
    candidates = candidate_pool(FOCAL.species, FOCAL.lineage, POOL + [same_family], scope_rank='ORDER')
    ranked = rank_candidates(candidates)
    assert [c.species for c in ranked] == ['Mucor mucedo', 'Lichtheimia corymbifera', 'Rhizopus arrhizus']


def test_select_nearest_takes_top_n():
    picked = select_nearest(FOCAL.species, FOCAL.lineage, POOL, n=1, scope_rank='ORDER')
    assert len(picked) == 1
    # deterministic tiebreak: both tie at rank ORDER, 'Lichtheimia' sorts before 'Rhizopus'
    assert picked[0].species == 'Lichtheimia corymbifera'


def test_select_nearest_backfills_past_an_excluded_top_candidate():
    # Mucor mucedo shares GENUS with the focal (closer than Rhizopus/Lichtheimia,
    # which only share ORDER) -- excluding it must NOT just shrink the result to
    # fewer than n; it must promote the next-best candidate instead.
    same_genus = _sample('Mucor mucedo', 'Mucoromycota', 'Mucoromycetes', 'Mucorales', 'Mucoraceae', 'Mucor')
    picked = select_nearest(
        FOCAL.species, FOCAL.lineage, POOL + [same_genus], n=1,
        scope_rank='ORDER', excluded={'Mucor mucedo'},
    )
    assert [c.species for c in picked] == ['Lichtheimia corymbifera']


def test_select_trait_filters_then_ranks():
    traits_by_species = {
        'Rhizopus arrhizus': [_TraitRow('Rhizopus arrhizus', 'thermotolerance', 'high')],
    }
    picked = select_trait(
        FOCAL.species, FOCAL.lineage, POOL, 'thermotolerance', 'high', n=2,
        traits_by_species=traits_by_species, scope_rank='ORDER',
    )
    assert [c.species for c in picked] == ['Rhizopus arrhizus']


def test_select_trait_errors_on_empty_filtered_pool():
    with pytest.raises(SystemExit, match='no candidate'):
        select_trait(
            FOCAL.species, FOCAL.lineage, POOL, 'thermotolerance', 'high', n=2,
            traits_by_species={}, scope_rank='ORDER',
        )


def test_exclude_species_removes_outgroup_pool_members():
    candidates = rank_candidates(candidate_pool(FOCAL.species, FOCAL.lineage, POOL, scope_rank='ORDER'))
    remaining = exclude_species(candidates, {'Rhizopus arrhizus'})
    assert [c.species for c in remaining] == ['Lichtheimia corymbifera']


def _TraitRow(species, trait, value):
    from trait_data import TraitRow
    return TraitRow(species=species, trait=trait, value=value, source='', notes='')
