"""Species trait loader: config_support/traits/trait_definitions.yaml (the
controlled vocabulary) + config_support/traits/traits.csv (the actual
species-to-trait rows). See config_support/traits/README.md and
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Trait data").
"""
import csv
from dataclasses import dataclass

import yaml


@dataclass
class TraitRow:
    species: str
    trait: str
    value: str
    source: str
    notes: str


def load_trait_definitions(path) -> dict[str, set[str]]:
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    return {trait: set(spec['values'].keys()) for trait, spec in doc['traits'].items()}


def load_traits(path, definitions: dict[str, set[str]]) -> dict[str, list[TraitRow]]:
    by_species: dict[str, list[TraitRow]] = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            trait = row['trait'].strip()
            value = row['value'].strip()
            species = row['Species'].strip()
            if trait not in definitions:
                raise SystemExit(f"{path}: undeclared trait {trait!r} (species: {species!r})")
            if value not in definitions[trait]:
                raise SystemExit(
                    f"{path}: undeclared value {value!r} for trait {trait!r} (species: {species!r})"
                )
            by_species.setdefault(species, []).append(TraitRow(
                species=species, trait=trait, value=value,
                source=(row.get('source') or '').strip(),
                notes=(row.get('notes') or '').strip(),
            ))

    for species, rows in by_species.items():
        values_by_trait: dict[str, set[str]] = {}
        for r in rows:
            values_by_trait.setdefault(r.trait, set()).add(r.value)
        for trait, values in values_by_trait.items():
            if 'none' in values and len(values) > 1:
                raise SystemExit(
                    f"{path}: {species!r} trait {trait!r} has 'none' coexisting with "
                    f"{sorted(values - {'none'})} -- 'none' must never coexist with another value"
                )

    return by_species


def has_trait(traits_by_species: dict[str, list[TraitRow]], species: str, trait: str, value: str) -> bool:
    return any(r.trait == trait and r.value == value for r in traits_by_species.get(species, []))
