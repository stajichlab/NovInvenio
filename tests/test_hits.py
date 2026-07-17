import gzip
import io
import tempfile
from pathlib import Path

from hits import (
    Hit,
    PARSERS,
    filter_hits,
    open_input,
    parse_phmmer_tblout,
    parse_tabular,
)


PHMMER_TBLOUT = """\
#                                                                            --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ----------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
geneA                -          query1               -           1.2e-50   200.0   0.0   2.4e-50  199.0   0.0   1.0   1   0   0   1   1   1   1 desc
geneB                -          query2               -           5.0e-03    10.0   0.0   1.0e-02    9.0   0.0   1.0   1   0   0   1   1   1   1 desc
"""

TABULAR = """\
q1\tt1\t1e-20\t150
q2\tt2\t1e-2\t12
# comment
q3\tt3\t3e-10\t75
"""


def test_parse_phmmer():
    hits = list(parse_phmmer_tblout(io.StringIO(PHMMER_TBLOUT)))
    assert len(hits) == 2
    assert hits[0] == Hit('query1', 'geneA', 1.2e-50, 200.0)
    assert hits[1].evalue == 5.0e-03


def test_parse_tabular_skips_comments():
    hits = list(parse_tabular(io.StringIO(TABULAR)))
    assert [h.query_id for h in hits] == ['q1', 'q2', 'q3']
    assert hits[0].bitscore == 150.0


def test_filter_hits_evalue():
    hits = list(parse_tabular(io.StringIO(TABULAR)))
    kept = list(filter_hits(iter(hits), evalue_cutoff=1e-5))
    assert [h.query_id for h in kept] == ['q1', 'q3']


def test_filter_hits_bitscore():
    hits = list(parse_tabular(io.StringIO(TABULAR)))
    kept = list(filter_hits(iter(hits), evalue_cutoff=1.0, bitscore_cutoff=100))
    assert [h.query_id for h in kept] == ['q1']


def test_parser_registry_keys():
    assert set(PARSERS) == {'phmmer', 'diamond', 'blast'}


def test_open_input_plain(tmp_path):
    f = tmp_path / "hits.tsv"
    f.write_text("q1\tt1\t1e-10\t100\n")
    with open_input(f) as fh:
        assert fh.read().startswith("q1")


def test_open_input_gz(tmp_path):
    f = tmp_path / "hits.tsv.gz"
    with gzip.open(f, 'wt') as fh:
        fh.write("q1\tt1\t1e-10\t100\n")
    with open_input(f) as fh:
        assert fh.read().startswith("q1")
