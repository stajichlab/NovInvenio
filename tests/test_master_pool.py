import pytest
from master_pool import (
    MASTER_POOL_FIELDS,
    MasterSample,
    assign_shorts,
    load_master_pool,
    make_short,
)

POOL_CSV = """\
Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID
Mucor circinelloides,1006PhL,/data/Mucor_circinelloides.pep.fa,/data/Mucor_circinelloides.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Mucoraceae;Mucor,36698
Rhizopus arrhizus,,/data/Rhizopus_arrhizus.pep.fa,/data/Rhizopus_arrhizus.dna.fa,Mucoromycota;Mucoromycotina;Mucoromycetes;;Mucorales;Rhizopodaceae;Rhizopus,64495
"""


def test_load_master_pool_parses_fixed_width_lineage(tmp_path):
    p = tmp_path / 'pool.csv'
    p.write_text(POOL_CSV)
    samples = load_master_pool(p)
    assert len(samples) == 2
    mucor = samples[0]
    assert mucor.species == 'Mucor circinelloides'
    assert mucor.strain == '1006PhL'
    assert mucor.protein_path == '/data/Mucor_circinelloides.pep.fa'
    assert mucor.lineage == [
        'Mucoromycota', 'Mucoromycotina', 'Mucoromycetes', '', 'Mucorales', 'Mucoraceae', 'Mucor',
    ]
    assert mucor.ncbi_taxid == '36698'


def test_load_master_pool_rejects_wrong_lineage_width(tmp_path):
    p = tmp_path / 'bad_pool.csv'
    p.write_text(
        "Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID\n"
        "Mucor circinelloides,,x,y,Mucoromycota;Mucorales;Mucor,36698\n"
    )
    with pytest.raises(SystemExit, match='expected 7'):
        load_master_pool(p)


def test_load_master_pool_rejects_duplicate_species(tmp_path):
    p = tmp_path / 'dup_pool.csv'
    p.write_text(POOL_CSV + POOL_CSV.splitlines()[1] + "\n")
    with pytest.raises(SystemExit, match='Duplicate Species'):
        load_master_pool(p)


def test_make_short_disambiguates_collisions():
    used: set[str] = set()
    # genus[:4] upper-first + lower + epithet[:4], concatenated and truncated to 8:
    # 'Muco' + 'circ' -> 'Mucocirc' (already 8 chars, no truncation)
    assert make_short('Mucor circinelloides', used) == 'Mucocirc'
    assert make_short('Mucor mucedo', used) == 'Mucomuce'
    # third species collides on the same 8-char base and gets a numeric suffix
    assert make_short('Mucor circinans', used) == 'Mucocir2'


def test_assign_shorts_is_deterministic_regardless_of_input_order(tmp_path):
    p = tmp_path / 'pool.csv'
    p.write_text(POOL_CSV)
    samples = load_master_pool(p)
    forward = assign_shorts(samples)
    backward = assign_shorts(list(reversed(samples)))
    assert forward == backward == {
        'Mucor circinelloides': 'Mucocirc',
        'Rhizopus arrhizus': 'Rhizarrh',
    }
