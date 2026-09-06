import pytest
from master_pool import (
    assign_shorts,
    load_master_pool,
    load_representative_picks,
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


REPR_TSV = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Mucor_circinelloides\tMucor circinelloides\tFalse\tMucor_circinelloides_1006PhL\t99.46\tTrue\n"
    "Mucor_circinelloides_1006PhL\tMucor circinelloides\tTrue\tMucor_circinelloides_1006PhL\t100.0\tTrue\n"
    "Rhizopus_arrhizus\tRhizopus arrhizus\tTrue\tRhizopus_arrhizus\t100.0\tFalse\n"
)

REPR_TSV_ZERO_TRUE = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Phycomyces_blakesleeanus\tPhycomyces blakesleeanus\tFalse\tPhycomyces_blakesleeanus_NRRL1555\t98.0\tTrue\n"
)

REPR_TSV_MULTIPLE_TRUE = (
    "out\tspecies\tis_representative\trepresentative_out\tani_to_representative\treuse_eligible\n"
    "Basidiobolus_A\tBasidiobolus meristosporus\tTrue\tBasidiobolus_A\t100.0\tFalse\n"
    "Basidiobolus_B\tBasidiobolus meristosporus\tTrue\tBasidiobolus_A\t99.9\tFalse\n"
)


def test_load_representative_picks(tmp_path):
    p = tmp_path / 'repr.tsv'
    p.write_text(REPR_TSV)
    picks = load_representative_picks(p)
    assert picks == {
        'Mucor circinelloides': 'Mucor_circinelloides_1006PhL',
        'Rhizopus arrhizus': 'Rhizopus_arrhizus',
    }


def test_load_representative_picks_rejects_zero_true_rows(tmp_path):
    p = tmp_path / 'repr_zero.tsv'
    p.write_text(REPR_TSV_ZERO_TRUE)
    with pytest.raises(SystemExit, match='Phycomyces blakesleeanus'):
        load_representative_picks(p)


def test_load_representative_picks_rejects_multiple_true_rows(tmp_path):
    p = tmp_path / 'repr_multi.tsv'
    p.write_text(REPR_TSV_MULTIPLE_TRUE)
    with pytest.raises(SystemExit, match='Basidiobolus meristosporus'):
        load_representative_picks(p)
