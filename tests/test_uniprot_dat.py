import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "lib"))
import hashlib  # noqa: E402
from uniprot_dat import RECORD_COLUMNS, normalize_seq, parse_dat_gz, seq_md5  # noqa: E402

FIXTURE_DAT = """\
ID   A7UWL5_NEUCR            Unreviewed;       6 AA.
AC   A7UWL5;
DE   SubName: Full=NRS/ER;
GN   ORFNames=NCU10683;
OX   NCBI_TaxID=367110;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RX   PubMed=12712197; DOI=10.1038/nature01554;
RA   Galagan J.E.;
RT   "The genome sequence of the filamentous fungus Neurospora crassa.";
RL   Nature 422:859-868(2003).
RN   [2]
RP   FUNCTION; SUBCELLULAR LOCATION.
RX   PubMed=99999999;
RT   "A title; with a semicolon | and a pipe
RT   that wraps.";
RL   J. Test 1:1-2(2020).
DR   VEuPathDB; FungiDB:NCU10683; -.
DR   GeneID; 5847462; -.
DR   RefSeq; XP_001728253.1; XM_001728201.2.
DR   RefSeq; XP_999999999.1; XM_999999999.1.
DR   KEGG; ncr:NCU10683; -.
DR   EnsemblFungi; EAA34304; EAA34304; NCU06768.
DR   PDB; 6YWS; EM; 2.74 A; K=1-249.
DR   STRING; 367110.A7UWL5; -.
DR   PANTHER; PTHR10000; PHOSPHOSERINE PHOSPHATASE; 1.
DR   UnknownDB; zzz; -.
DR   Pfam; PF04321; RmlD_sub_bind; 1.
DR   GO; GO:0016020; C:membrane; IEA:UniProtKB-KW.
PE   3: Inferred from homology;
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTE I
//
ID   B0000000_TEST           Reviewed;         10 AA.
AC   B0000000;
GN   ORFNames=NCU99999;
OX   NCBI_TaxID=367110;
DE   RecName: Full=Test with Scer-style Ensembl field order;
DR   EnsemblFungi; YML051W_mRNA; YML051W; YML051W.
PE   1: Evidence at protein level;
SQ   SEQUENCE   6 AA;  1 MW;  0000000000000000 CRC64;
     MATTEI
//
ID   C0000000_TEST           Unreviewed;       6 AA.
AC   C0000000;
OX   NCBI_TaxID=367110;
DE   SubName: Full=No sequence block;
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


def test_record_has_every_column_in_order(tmp_path):
    rec = next(iter(parse_dat_gz(_write_fixture(tmp_path))))
    assert tuple(rec.keys()) == RECORD_COLUMNS


def test_id_line_reviewed_entry_name_and_pe(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    assert recs["A7UWL5"]["entry_name"] == "A7UWL5_NEUCR"
    assert recs["A7UWL5"]["reviewed"] == "0"
    assert recs["B0000000"]["reviewed"] == "1"
    assert recs["A7UWL5"]["protein_existence"] == "3"
    assert recs["B0000000"]["protein_existence"] == "1"


def test_sequence_md5_and_length(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    assert recs["A7UWL5"]["seq_md5"] == hashlib.md5(b"MATTEI").hexdigest()
    assert recs["A7UWL5"]["seq_len"] == "6"


def test_normalize_seq_strips_case_whitespace_and_stop():
    assert normalize_seq(" matt\nei* ") == "MATTEI"
    assert seq_md5("mattei*") == hashlib.md5(b"MATTEI").hexdigest()


def test_record_without_sequence_is_skipped_and_counted(tmp_path):
    stats = {}
    recs = {r["accession"] for r in parse_dat_gz(_write_fixture(tmp_path), stats=stats)}
    assert "C0000000" not in recs
    assert stats == {"n_records": 2, "n_skipped_no_ac": 0, "n_skipped_no_seq": 1}


def test_pubs_scope_and_title_with_separators(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    pubs = recs["A7UWL5"]["pubs"].split("|")
    assert pubs[0] == ("12712197;10.1038/nature01554;proteome;"
                       "The genome sequence of the filamentous fungus Neurospora crassa.")
    pmid, doi, scope, title = pubs[1].split(";", 3)
    assert (pmid, doi, scope) == ("99999999", "", "protein")
    assert title == "A title; with a semicolon / and a pipe that wraps."
    assert recs["B0000000"]["pubs"] == ""


def test_xrefs_keep_new_allowlist_and_drop_unknown(tmp_path):
    recs = {r["accession"]: r for r in parse_dat_gz(_write_fixture(tmp_path))}
    x = recs["A7UWL5"]["xrefs"].split("|")
    assert "PDB:6YWS" in x
    assert "STRING:367110.A7UWL5" in x
    assert "PANTHER:PTHR10000" in x
    assert not any(e.startswith("UnknownDB") for e in x)
