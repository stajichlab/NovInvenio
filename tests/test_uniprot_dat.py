import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
from uniprot_dat import parse_dat_gz  # noqa: E402

FIXTURE_DAT = """\
ID   A7UWL5_NEUCR            Unreviewed;       301 AA.
AC   A7UWL5;
GN   ORFNames=NCU10683;
OX   NCBI_TaxID=367110;
DE   SubName: Full=NRS/ER;
DR   VEuPathDB; FungiDB:NCU10683; -.
DR   GeneID; 5847462; -.
DR   RefSeq; XP_001728253.1; XM_001728201.2.
DR   RefSeq; XP_999999999.1; XM_999999999.1.
DR   KEGG; ncr:NCU10683; -.
DR   EnsemblFungi; EAA34304; EAA34304; NCU06768.
DR   PDB; 6YWS; EM; 2.74 A; K=1-249.
DR   STRING; 367110.A7UWL5; -.
DR   Pfam; PF04321; RmlD_sub_bind; 1.
DR   GO; GO:0016020; C:membrane; IEA:UniProtKB-KW.
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
ID   B0000000_TEST           Unreviewed;       10 AA.
AC   B0000000;
GN   ORFNames=NCU99999;
OX   NCBI_TaxID=367110;
DE   SubName: Full=Test with Scer-style Ensembl field order;
DR   EnsemblFungi; YML051W_mRNA; YML051W; YML051W.
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
"""


def _write_fixture(tmp_path) -> Path:
    p = tmp_path / "fixture.dat.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(FIXTURE_DAT)
    return p


def test_xrefs_packs_allowlisted_db_types_first_colon_only(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"].split("|")
    assert "VEuPathDB:FungiDB:NCU10683" in xrefs
    assert "GeneID:5847462" in xrefs
    assert "KEGG:ncr:NCU10683" in xrefs


def test_xrefs_keeps_only_first_refseq_line(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"].split("|")
    refseq_entries = [x for x in xrefs if x.startswith("RefSeq:")]
    assert refseq_entries == ["RefSeq:XP_001728253.1"]


def test_xrefs_ensemblfungi_takes_last_field_regardless_of_layout(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    a7 = records["A7UWL5"]["xrefs"].split("|")
    b0 = records["B0000000"]["xrefs"].split("|")
    assert "EnsemblFungi:NCU06768" in a7
    assert "EnsemblFungi:YML051W" in b0


def test_xrefs_drops_pdb_and_non_allowlisted_dbs(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    xrefs = records["A7UWL5"]["xrefs"]
    assert "PDB" not in xrefs
    assert "STRING" not in xrefs


def test_pfam_go_and_description_fields(tmp_path):
    path = _write_fixture(tmp_path)
    records = {r["accession"]: r for r in parse_dat_gz(path)}
    a7 = records["A7UWL5"]
    assert a7["pfam_ids"] == "PF04321"
    assert a7["pfam_names"] == "RmlD_sub_bind"
    assert a7["go_ids"] == "GO:0016020:IEA"
    assert a7["description"] == "NRS/ER"
    assert a7["gene_name"] == "NCU10683"
    assert a7["taxon_id"] == "367110"
