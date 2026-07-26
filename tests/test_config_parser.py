from config_parser import (
    GROUPS,
    get_broad_out,
    get_disc_out,
    get_group,
    get_ingroup,
    get_near_in,
    get_outgroup,
    get_target,
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
TARGET,Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Pezizomycotina
TARGET,Pyronema omphalodes,CBS144459,Pomp.pep.fa,,Pomp,Pezizomycotina
DISC_OUT,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
DISC_OUT,Saccharomyces cerevisiae,S288C,Scer.pep.fa,Scer.dna.fa,Scer,Saccharomycotina
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
    assert GROUPS == {'IN', 'OUT', 'TARGET', 'DISC_OUT', 'NEAR_IN', 'BROAD_OUT'}


def test_novelty_group_helpers(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert [s.short for s in get_target(samples)] == ['Ncra', 'Pomp']
    assert [s.short for s in get_disc_out(samples)] == ['Afum', 'Scer']
    assert [s.short for s in get_near_in(samples)] == ['Cimm']
    assert [s.short for s in get_broad_out(samples)] == ['Spom']


def test_get_group_generic(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert [s.short for s in get_group(samples, 'TARGET')] == ['Ncra', 'Pomp']
    assert [s.short for s in get_group(samples, 'DISC_OUT')] == ['Afum', 'Scer']


def test_short_to_group_with_novelty_roles(tmp_path):
    p = tmp_path / 'cfg.csv'
    p.write_text(NOVELTY_CSV)
    samples = parse_config(p)
    assert short_to_group(samples) == {
        'Ncra': 'TARGET',
        'Pomp': 'TARGET',
        'Afum': 'DISC_OUT',
        'Scer': 'DISC_OUT',
        'Cimm': 'NEAR_IN',
        'Spom': 'BROAD_OUT',
    }
