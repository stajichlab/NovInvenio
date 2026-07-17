from config_parser import (
    get_ingroup,
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
