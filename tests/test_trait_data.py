import pytest
from trait_data import has_trait, load_trait_definitions, load_traits

DEFS_YAML = """\
traits:
  spore_motility:
    description: motility
    values:
      motile:
        description: flagellated
        ontology_term:
      nonmotile:
        description: not flagellated
        ontology_term:
  animal_association:
    description: host association
    values:
      none:
        description: none
        ontology_term:
      pathogen:
        description: pathogen
        ontology_term:
"""

TRAITS_CSV = """\
Species,trait,value,source,notes
Mucor circinelloides,spore_motility,nonmotile,,sporangiospores
Mucor circinelloides,animal_association,pathogen,,opportunistic
"""


def test_load_trait_definitions():
    defs = load_trait_definitions_from_text(DEFS_YAML)
    assert defs == {
        'spore_motility': {'motile', 'nonmotile'},
        'animal_association': {'none', 'pathogen'},
    }


def load_trait_definitions_from_text(text, tmp_path=None):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / 'defs.yaml'
        p.write_text(text)
        return load_trait_definitions(p)


def test_load_traits_and_has_trait(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text(TRAITS_CSV)

    defs = load_trait_definitions(defs_path)
    by_species = load_traits(traits_path, defs)

    assert has_trait(by_species, 'Mucor circinelloides', 'animal_association', 'pathogen')
    assert not has_trait(by_species, 'Mucor circinelloides', 'animal_association', 'none')
    assert not has_trait(by_species, 'Phycomyces blakesleeanus', 'spore_motility', 'motile')


def test_load_traits_rejects_undeclared_trait(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text("Species,trait,value,source,notes\nMucor circinelloides,made_up_trait,x,,\n")

    with pytest.raises(SystemExit, match='undeclared trait'):
        load_traits(traits_path, load_trait_definitions(defs_path))


def test_load_traits_rejects_undeclared_value(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text("Species,trait,value,source,notes\nMucor circinelloides,spore_motility,flying,,\n")

    with pytest.raises(SystemExit, match='undeclared value'):
        load_traits(traits_path, load_trait_definitions(defs_path))


def test_load_traits_rejects_none_coexisting_with_another_value(tmp_path):
    defs_path = tmp_path / 'defs.yaml'
    defs_path.write_text(DEFS_YAML)
    traits_path = tmp_path / 'traits.csv'
    traits_path.write_text(
        "Species,trait,value,source,notes\n"
        "Mucor circinelloides,animal_association,none,,\n"
        "Mucor circinelloides,animal_association,pathogen,,\n"
    )

    with pytest.raises(SystemExit, match="'none' must never coexist"):
        load_traits(traits_path, load_trait_definitions(defs_path))
