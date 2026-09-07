import pytest
from source_db import load_source_db

VALID_CSV = """\
Species,SourceDB,notes
Neurospora crassa,fungidb,FungiDB release-68
Mucor circinelloides,mycocosm:Mucci1,JGI MycoCosm portal
Coprinopsis cinerea,ensemblfungi:coprinopsis_cinerea,EnsemblFungi
Cryptosporidium parvum,veupathdb:CryptoDB,VEuPathDB project
Drosophila melanogaster,ncbipep,plain NCBI RefSeq proteome with no FungiDB/MycoCosm record
Custom species,https://example.org/gene/{gene},raw URL template
"""

INVALID_JAVASCRIPT_CSV = """\
Species,SourceDB,notes
Evil species,javascript:alert(1),malicious
"""

INVALID_UNKNOWN_PREFIX_CSV = """\
Species,SourceDB,notes
Some species,notarealdb,typo
"""

INVALID_EMPTY_PORTAL_CSV = """\
Species,SourceDB,notes
Some species,mycocosm:,missing portal argument
"""


def test_load_source_db_accepts_all_known_forms(tmp_path):
    p = tmp_path / 'source_db.csv'
    p.write_text(VALID_CSV)
    mapping = load_source_db(p)
    assert mapping == {
        'Neurospora crassa': 'fungidb',
        'Mucor circinelloides': 'mycocosm:Mucci1',
        'Coprinopsis cinerea': 'ensemblfungi:coprinopsis_cinerea',
        'Cryptosporidium parvum': 'veupathdb:CryptoDB',
        'Drosophila melanogaster': 'ncbipep',
        'Custom species': 'https://example.org/gene/{gene}',
    }


def test_load_source_db_rejects_javascript_scheme(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text(INVALID_JAVASCRIPT_CSV)
    with pytest.raises(SystemExit, match='Evil species'):
        load_source_db(p)


def test_load_source_db_rejects_unknown_prefix(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text(INVALID_UNKNOWN_PREFIX_CSV)
    with pytest.raises(SystemExit, match='Some species'):
        load_source_db(p)


def test_load_source_db_rejects_empty_portal_argument(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text(INVALID_EMPTY_PORTAL_CSV)
    with pytest.raises(SystemExit, match='Some species'):
        load_source_db(p)
