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

    map_path = Path(summary['map_path'])
    assert map_path.exists()
    map_text = map_path.read_text()
    assert 'Mucor circinelloides\t' in map_text


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
