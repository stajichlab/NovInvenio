"""mode: nearest / trait / explicit candidate selection, and the
ingroup/outgroup disjointness rule. See
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Trait mode",
"Ingroup/outgroup disjointness").
"""
from dataclasses import dataclass

from lineage import RANK_NAMES, lineage_match
from trait_data import has_trait

DEFAULT_SCOPE_RANK = 'ORDER'


@dataclass
class Candidate:
    species: str
    rank_name: str  # deepest shared rank name with the focal, e.g. 'FAMILY'


def _depth(rank_name: str) -> int:
    """Sort key: deeper rank = smaller index = closer. '' (no match) sorts last."""
    return RANK_NAMES.index(rank_name) if rank_name else -1


def candidate_pool(focal_species: str, focal_lineage: list[str], pool: list, scope_rank: str = DEFAULT_SCOPE_RANK) -> list[Candidate]:
    """Every species in `pool` (excluding the focal itself) sharing at
    least rank `scope_rank` with the focal."""
    if scope_rank not in RANK_NAMES:
        raise ValueError(f"scope_rank {scope_rank!r} must be one of {RANK_NAMES}")
    scope_depth = RANK_NAMES.index(scope_rank)
    out = []
    for s in pool:
        if s.species == focal_species:
            continue
        rank_name = lineage_match(focal_lineage, s.lineage)
        if _depth(rank_name) >= scope_depth:
            out.append(Candidate(species=s.species, rank_name=rank_name))
    return out


def rank_candidates(candidates: list[Candidate]) -> list[Candidate]:
    """Deepest match first; ties broken alphabetically by species name."""
    return sorted(candidates, key=lambda c: (-_depth(c.rank_name), c.species))


def exclude_species(candidates: list[Candidate], excluded: set[str]) -> list[Candidate]:
    return [c for c in candidates if c.species not in excluded]


def select_nearest(focal_species, focal_lineage, pool, n, scope_rank: str = DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]:
    # excluded is applied BEFORE ranking/truncation -- see this task's Interfaces
    # note: excluding after truncation would silently under-fill a study instead
    # of backfilling from the next-best candidate.
    candidates = exclude_species(candidate_pool(focal_species, focal_lineage, pool, scope_rank), excluded)
    return rank_candidates(candidates)[:n]


def select_trait(focal_species, focal_lineage, pool, trait, value, n, traits_by_species, scope_rank: str = DEFAULT_SCOPE_RANK, excluded: set[str] = frozenset()) -> list[Candidate]:
    candidates = exclude_species(candidate_pool(focal_species, focal_lineage, pool, scope_rank), excluded)
    ranked = rank_candidates(candidates)
    filtered = [c for c in ranked if has_trait(traits_by_species, c.species, trait, value)]
    if not filtered:
        raise SystemExit(
            f"mode: trait -- no candidate for focal {focal_species!r} in the "
            f"{scope_rank}-scoped pool has {trait}={value!r} (after excluding the outgroup pool)"
        )
    return filtered[:n]
