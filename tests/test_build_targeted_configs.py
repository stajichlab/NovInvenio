import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))
from build_targeted_configs import render_batch  # noqa: E402

MASTER_POOL_CSV = """\
Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID
Mucor circinelloides,1006PhL,/data/Mucci.pep.fa,/data/Mucci.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor,36698
Rhizopus arrhizus,,/data/Rhiar.pep.fa,/data/Rhiar.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Rhizopodaceae;Rhizopus,64495
Lichtheimia corymbifera,,/data/Licor.pep.fa,/data/Licor.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Lichtheimiaceae;Lichtheimia,64752
Basidiobolus meristosporus,,/data/Bamer.pep.fa,/data/Bamer.dna.fa,Basidiobolomycota;;Basidiobolomycetes;;Basidiobolales;Basidiobolaceae;Basidiobolus,,
Neurospora crassa,OR74A,/data/Ncra.pep.fa,/data/Ncra.dna.fa,Ascomycota;Pezizomycotina;Sordariomycetes;;Sordariales;Sordariaceae;Neurospora,367110
Aspergillus nidulans,,/data/Anid.pep.fa,/data/Anid.dna.fa,Ascomycota;Pezizomycotina;Eurotiomycetes;;Eurotiales;Aspergillaceae;Aspergillus,227321
"""

TRAIT_DEFINITIONS_YAML = """\
traits:
  thermotolerance:
    description: x
    values:
      high: {description: x, ontology_term:}
      low: {description: x, ontology_term:}
"""

TRAITS_CSV = """\
Species,trait,value,source,notes
Rhizopus arrhizus,thermotolerance,high,,
"""

BATCH_YAML = """\
batch: test_batch_v1
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 2}
    outgroup_pool: dikarya_v1
"""

BATCH_YAML_TRAIT_MODE = """\
batch: test_batch_v2
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: trait, trait: thermotolerance, value: high, n: 1}
    outgroup_pool: dikarya_v1
"""


def _write_fixtures(tmp_path, batch_yaml=BATCH_YAML):
    pool = tmp_path / 'pool.csv'
    pool.write_text(MASTER_POOL_CSV)
    defs = tmp_path / 'trait_definitions.yaml'
    defs.write_text(TRAIT_DEFINITIONS_YAML)
    traits = tmp_path / 'traits.csv'
    traits.write_text(TRAITS_CSV)
    batch = tmp_path / 'batch.yaml'
    batch.write_text(batch_yaml)
    outdir = tmp_path / 'out'
    outdir.mkdir()
    return pool, defs, traits, batch, outdir


def test_render_batch_nearest_mode_produces_config_and_map(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    summaries = render_batch(pool, defs, traits, batch, outdir)

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary['focal'] == 'Mucor circinelloides'
    assert {c['species'] for c in summary['companions']} == {'Rhizopus arrhizus', 'Lichtheimia corymbifera'}

    config_path = Path(summary['config_path'])
    assert config_path.exists()
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    groups = {r['Species']: r['GROUP'] for r in rows}
    assert groups == {
        'Mucor circinelloides': 'IN',
        'Rhizopus arrhizus': 'IN',
        'Lichtheimia corymbifera': 'IN',
        'Neurospora crassa': 'OUT',
        'Aspergillus nidulans': 'OUT',
    }

    # TaxonGroup for a nearest/trait companion is the actual taxon NAME at the
    # matched rank (e.g. "Mucorales"), not the rank category label ("ORDER").
    taxon_groups = {r['Species']: r['TaxonGroup'] for r in rows}
    assert taxon_groups['Rhizopus arrhizus'] == 'Mucorales'
    assert taxon_groups['Lichtheimia corymbifera'] == 'Mucorales'

    map_path = Path(summary['map_path'])
    assert map_path.exists()
    map_text = map_path.read_text()
    assert 'Mucor circinelloides\t' in map_text


SOURCE_DB_CSV = """\
Species,SourceDB,notes
Neurospora crassa,fungidb,FungiDB release-68
"""


def test_render_batch_populates_source_db_for_listed_species_only(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    source_db = tmp_path / 'source_db.csv'
    source_db.write_text(SOURCE_DB_CSV)

    summaries = render_batch(pool, defs, traits, batch, outdir, source_db_path=source_db)

    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    source_db_by_species = {r['Species']: r['SourceDB'] for r in rows}

    assert source_db_by_species['Neurospora crassa'] == 'fungidb'
    # Not in the seed file -> empty, not a KeyError or a fabricated value.
    assert source_db_by_species['Mucor circinelloides'] == ''
    assert source_db_by_species['Aspergillus nidulans'] == ''


def test_render_batch_without_source_db_arg_leaves_column_empty(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    summaries = render_batch(pool, defs, traits, batch, outdir)

    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    assert all(r['SourceDB'] == '' for r in rows)


ANIMAL_POOL_CSV = """\
Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID
Drosophila melanogaster,ISO1,/data/Dmel.pep.fa,,Arthropoda;Hexapoda;Insecta;Pterygota;Diptera;Drosophilidae;Drosophila,7227
"""


def test_render_batch_merges_extra_pool_for_explicit_mode(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    extra_pool = tmp_path / 'animal_pool.csv'
    extra_pool.write_text(ANIMAL_POOL_CSV)

    batch_yaml = """\
batch: animal_v1
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: explicit, members: ["Drosophila melanogaster"], reason: "cross-kingdom trait comparison"}
    outgroup_pool: dikarya_v1
"""
    batch.write_text(batch_yaml)

    summaries = render_batch(pool, defs, traits, batch, outdir, extra_pool_paths=[extra_pool])

    assert summaries[0]['companions'] == [
        {'species': 'Drosophila melanogaster', 'taxon_group': 'Drosophila', 'reason': 'cross-kingdom trait comparison'}
    ]
    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    groups = {r['Species']: r['GROUP'] for r in rows}
    assert groups['Drosophila melanogaster'] == 'IN'


def test_render_batch_errors_on_species_collision_between_pools(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    extra_pool = tmp_path / 'colliding_pool.csv'
    extra_pool.write_text(
        "Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID\n"
        "Mucor circinelloides,,,,;;;;;;,\n"
    )
    with pytest.raises(SystemExit, match='Mucor circinelloides'):
        render_batch(pool, defs, traits, batch, outdir, extra_pool_paths=[extra_pool])


def test_render_batch_trait_mode_filters_to_matching_candidate(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=BATCH_YAML_TRAIT_MODE)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    assert [c['species'] for c in summaries[0]['companions']] == ['Rhizopus arrhizus']


def test_render_batch_errors_when_focal_is_inside_its_own_outgroup_pool(tmp_path):
    bad_batch = """\
batch: bad
outgroup_pools:
  bad_pool:
    members: ["Mucor circinelloides"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 1}
    outgroup_pool: bad_pool
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=bad_batch)
    with pytest.raises(SystemExit, match='inside its own outgroup pool'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_explicit_mode_renders_named_members(tmp_path):
    batch_yaml = """\
batch: explicit_v1
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: explicit, members: ["Basidiobolus meristosporus"], reason: "known outlier"}
    outgroup_pool: dikarya_v1
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    summaries = render_batch(pool, defs, traits, batch, outdir)

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary['companions'] == [
        {'species': 'Basidiobolus meristosporus', 'taxon_group': 'Basidiobolus', 'reason': 'known outlier'}
    ]

    config_path = Path(summary['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    groups = {r['Species']: r['GROUP'] for r in rows}
    assert groups == {
        'Mucor circinelloides': 'IN',
        'Basidiobolus meristosporus': 'IN',
        'Neurospora crassa': 'OUT',
        'Aspergillus nidulans': 'OUT',
    }


def test_render_batch_explicit_mode_errors_when_member_overlaps_outgroup_pool(tmp_path):
    batch_yaml = """\
batch: explicit_overlap
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: explicit, members: ["Neurospora crassa"]}
    outgroup_pool: dikarya_v1
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    with pytest.raises(SystemExit, match='also in outgroup pool'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_errors_when_outgroup_pool_name_undefined(tmp_path):
    batch_yaml = """\
batch: undefined_pool
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 1}
    outgroup_pool: does_not_exist
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    with pytest.raises(SystemExit, match='not defined'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_errors_when_outgroup_pool_member_not_in_master_pool(tmp_path):
    batch_yaml = """\
batch: unknown_pool_member
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Ghostus fungus"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 1}
    outgroup_pool: dikarya_v1
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    with pytest.raises(SystemExit, match='species not found in master pool'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_errors_when_explicit_member_not_in_master_pool(tmp_path):
    batch_yaml = """\
batch: unknown_explicit_member
outgroup_pools:
  dikarya_v1:
    members: ["Neurospora crassa", "Aspergillus nidulans"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: explicit, members: ["Ghostus fungus"]}
    outgroup_pool: dikarya_v1
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    with pytest.raises(SystemExit, match='species not found in master pool'):
        render_batch(pool, defs, traits, batch, outdir)


def test_render_batch_excludes_outgroup_members_from_ingroup_candidates(tmp_path):
    # Neurospora crassa can never be picked as an ORDER-scoped companion of
    # Mucor here anyway (different phylum), so this test uses a pool member
    # that WOULD otherwise rank -- Rhizopus -- placed into the outgroup pool.
    batch_yaml = """\
batch: disjoint
outgroup_pools:
  contains_rhizopus:
    members: ["Rhizopus arrhizus", "Neurospora crassa"]
studies:
  - focal: Mucor circinelloides
    ingroup_extra: {mode: nearest, n: 2}
    outgroup_pool: contains_rhizopus
"""
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path, batch_yaml=batch_yaml)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    assert [c['species'] for c in summaries[0]['companions']] == ['Lichtheimia corymbifera']


def test_render_batch_writes_ncbi_taxid_column(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    assert 'NCBI_TaxID' in rows[0].keys()
    taxids = {r['Species']: r['NCBI_TaxID'] for r in rows}
    assert taxids['Mucor circinelloides'] == '36698'
    assert taxids['Rhizopus arrhizus'] == '64495'


def test_render_batch_without_link_dir_writes_absolute_paths(tmp_path):
    pool, defs, traits, batch, outdir = _write_fixtures(tmp_path)
    summaries = render_batch(pool, defs, traits, batch, outdir)
    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    mucor_row = next(r for r in rows if r['Species'] == 'Mucor circinelloides')
    assert mucor_row['Protein'] == '/data/Mucci.pep.fa'
    assert Path(mucor_row['Protein']).is_absolute()


def test_render_batch_with_link_dir_writes_basenames_and_real_symlinks(tmp_path):
    # Use a master pool whose ProteinPath/DNAPath point at real files, so we
    # can verify the symlinks --link-dir creates actually resolve.
    data_dir = tmp_path / 'real_data'
    data_dir.mkdir()
    real_files = {}
    for short, species in [
        ('Mucci', 'Mucor circinelloides'), ('Rhiar', 'Rhizopus arrhizus'),
        ('Licor', 'Lichtheimia corymbifera'), ('Bamer', 'Basidiobolus meristosporus'),
        ('Ncra', 'Neurospora crassa'), ('Anid', 'Aspergillus nidulans'),
    ]:
        pep = data_dir / f'{short}.pep.fa'
        dna = data_dir / f'{short}.dna.fa'
        pep.write_text('>x\nMKV\n')
        dna.write_text('>x\nACGT\n')
        real_files[species] = (pep, dna)

    pool_csv = "Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID\n"
    for species, lineage in [
        ('Mucor circinelloides', 'Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor'),
        ('Rhizopus arrhizus', 'Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Rhizopodaceae;Rhizopus'),
        ('Lichtheimia corymbifera', 'Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Lichtheimiaceae;Lichtheimia'),
        ('Basidiobolus meristosporus', 'Basidiobolomycota;;Basidiobolomycetes;;Basidiobolales;Basidiobolaceae;Basidiobolus'),
        ('Neurospora crassa', 'Ascomycota;Pezizomycotina;Sordariomycetes;;Sordariales;Sordariaceae;Neurospora'),
        ('Aspergillus nidulans', 'Ascomycota;Pezizomycotina;Eurotiomycetes;;Eurotiales;Aspergillaceae;Aspergillus'),
    ]:
        pep, dna = real_files[species]
        pool_csv += f"{species},,{pep},{dna},{lineage},\n"

    pool = tmp_path / 'pool.csv'
    pool.write_text(pool_csv)
    defs = tmp_path / 'trait_definitions.yaml'
    defs.write_text(TRAIT_DEFINITIONS_YAML)
    traits = tmp_path / 'traits.csv'
    traits.write_text(TRAITS_CSV)
    batch = tmp_path / 'batch.yaml'
    batch.write_text(BATCH_YAML)
    outdir = tmp_path / 'out'
    outdir.mkdir()
    link_dir = tmp_path / 'links'

    summaries = render_batch(pool, defs, traits, batch, outdir, link_dir=link_dir)

    config_path = Path(summaries[0]['config_path'])
    with open(config_path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    mucor_row = next(r for r in rows if r['Species'] == 'Mucor circinelloides')
    assert mucor_row['Protein'] == 'Mucci.pep.fa'
    assert mucor_row['DNA'] == 'Mucci.dna.fa'
    assert not Path(mucor_row['Protein']).is_absolute()

    pep_link = link_dir / 'pep' / 'Mucci.pep.fa'
    dna_link = link_dir / 'dna' / 'Mucci.dna.fa'
    assert pep_link.is_symlink()
    assert dna_link.is_symlink()
    assert pep_link.resolve() == real_files['Mucor circinelloides'][0].resolve()
    assert dna_link.resolve() == real_files['Mucor circinelloides'][1].resolve()
    # every IN/OUT row's file got linked, not just the focal
    for r in rows:
        assert (link_dir / 'pep' / r['Protein']).is_symlink()
        assert (link_dir / 'dna' / r['DNA']).is_symlink()
