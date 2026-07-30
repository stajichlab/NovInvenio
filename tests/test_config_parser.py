from config_parser import (
    GROUPS,
    get_broad_outgroup,
    get_discovery_out,
    get_discovery_target,
    get_group,
    get_ingroup,
    get_near_ingroup,
    get_outgroup,
    parse_config,
    short_to_group,
)


CSV = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
OUT,Neolecta irregularis,,Nirr.pep.fa,Nirr.dna.fa,Nirr,Taphrinomycotina
IN,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
IN,Pyronema omphalodes,CBS144459,Pomp.pep.fa,,Pomp,Pezizomycotina
"""

NOVELTY_CSV = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
DISCOVERY_TARGET,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
DISCOVERY_TARGET,Pyronema omphalodes,CBS144459,Pomp.pep.fa,,Pomp,Pezizomycotina
DISCOVERY_OUT,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
DISCOVERY_OUT,Saccharomyces cerevisiae,S288C,Scer.pep.fa,Scer.dna.fa,Scer,Pezizomycotina
NEAR_INGROUP,Coccidioides immitis,WA211,Cocci_WA211.pep.fa,Cocci_WA211.dna.fa,Cimm,Pezizomycotina
BROAD_OUTGROUP,Schizosaccharomyces pombe,,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
"""

# The original novelty_discovery/novelty_screen labels (issues #24-#29), still accepted as
# aliases for backward compatibility (todo/rename-novelty-discovery-group-labels.md).
NOVELTY_CSV_OLD_ALIASES = """\
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
TARGET,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
DISC_OUT,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
NEAR_IN,Coccidioides immitis,WA211,Cocci_WA211.pep.fa,Cocci_WA211.dna.fa,Cimm,Pezizomycotina
BROAD_OUT,Schizosaccharomyces pombe,,Spom.pep.fa,Spom.dna.fa,Spom,Taphrinomycotina
"""


def test_parse_config_roundtrip(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(CSV)
    samples = parse_config(p)
    assert [s.short for s in samples] == ['Nirr', 'Ncra', 'Pomp']
    assert samples[0].group == 'OUT'
    assert samples[0].strain == ''
    assert samples[2].dna == ''  # blank DNA cell


def test_group_splits(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(CSV)
    samples = parse_config(p)
    assert [s.short for s in get_ingroup(samples)] == ['Ncra', 'Pomp']
    assert [s.short for s in get_outgroup(samples)] == ['Nirr']


def test_short_to_group(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(CSV)
    samples = parse_config(p)
    assert short_to_group(samples) == {'Nirr': 'OUT', 'Ncra': 'IN', 'Pomp': 'IN'}


def test_groups_constant():
    assert GROUPS == {
        'IN', 'OUT', 'DISCOVERY_TARGET', 'DISCOVERY_OUT', 'NEAR_INGROUP', 'BROAD_OUTGROUP',
    }


def test_novelty_group_helpers(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert [s.short for s in get_discovery_target(samples)] == ['Ncra', 'Pomp']
    assert [s.short for s in get_discovery_out(samples)] == ['Afum', 'Scer']
    assert [s.short for s in get_near_ingroup(samples)] == ['Cimm']
    assert [s.short for s in get_broad_outgroup(samples)] == ['Spom']


def test_get_group_generic(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert [s.short for s in get_group(samples, 'DISCOVERY_TARGET')] == ['Ncra', 'Pomp']
    assert [s.short for s in get_group(samples, 'DISCOVERY_OUT')] == ['Afum', 'Scer']


def test_short_to_group_with_novelty_roles(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert short_to_group(samples) == {
        'Ncra': 'DISCOVERY_TARGET',
        'Pomp': 'DISCOVERY_TARGET',
        'Afum': 'DISCOVERY_OUT',
        'Scer': 'DISCOVERY_OUT',
        'Cimm': 'NEAR_INGROUP',
        'Spom': 'BROAD_OUTGROUP',
    }


def test_old_group_aliases_still_parse_and_normalize(tmp_path):
    # TARGET/DISC_OUT/NEAR_IN/BROAD_OUT (issues #24-#29) are renamed for clarity
    # (todo/rename-novelty-discovery-group-labels.md), but existing config CSVs using the
    # old labels must keep working unchanged -- parse_config() normalizes them.
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV_OLD_ALIASES)
    samples = parse_config(p)
    assert short_to_group(samples) == {
        'Ncra': 'DISCOVERY_TARGET',
        'Afum': 'DISCOVERY_OUT',
        'Cimm': 'NEAR_INGROUP',
        'Spom': 'BROAD_OUTGROUP',
    }
    assert [s.short for s in get_discovery_target(samples)] == ['Ncra']
