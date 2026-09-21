import gzip
import io

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


# --- Wider tabular format: length, pident, qcov, scov, qlen, slen (issue #129) ---
#
# diamond and blastp can both emit these fields (confirmed against real diamond
# 2.1.x / blastp 2.16.x binaries: `--outfmt 6 ... length pident qcovhsp scovhsp qlen
# slen` / `-outfmt "6 ... length pident qcovhsp qlen slen"`). parse_tabular must read
# them WHEN PRESENT and stay silently backward-compatible with the 4-column format
# already sitting in every project's search_cache/ -- storeDir never invalidates an
# existing cached file just because the code that would produce a NEW one changed, so
# old and new files coexist indefinitely (module docstring covers this).

TABULAR_WIDE_DIAMOND = (
    "q1\tt1\t1e-20\t150\t42\t35.5\t80\t90\t100\t45\n"      # full 10-col diamond layout
    "q2\tt2\t1e-2\t12\t10\t20.0\t50\t50\t50\t50\n"
)

# blastp's qcovhsp gives only QUERY coverage -- no scovhsp equivalent exists, so the
# blastp layout is one column narrower than diamond's.
TABULAR_WIDE_BLAST = "q1\tt1\t1e-20\t150\t42\t35.5\t80\t100\t45\n"


def test_parse_tabular_reads_wide_diamond_columns():
    hits = list(parse_tabular(io.StringIO(TABULAR_WIDE_DIAMOND)))
    h = hits[0]
    assert h.length == 42
    assert h.pident == 35.5
    assert h.qcov == 80.0
    assert h.scov == 90.0
    assert h.qlen == 100
    assert h.slen == 45


def test_parse_tabular_reads_wide_blast_columns_with_no_scov():
    hits = list(parse_tabular(io.StringIO(TABULAR_WIDE_BLAST)))
    h = hits[0]
    assert h.length == 42
    assert h.pident == 35.5
    assert h.qcov == 80.0
    assert h.scov is None      # blastp's outfmt has no subject-coverage keyword
    assert h.qlen == 100
    assert h.slen == 45


def test_parse_tabular_old_four_column_cache_files_still_load():
    # A file written before this change (or by a still-unwidened self-search variant)
    # must load exactly as before, with the new fields defaulting to None -- never 0,
    # which would look like a real "zero coverage" measurement.
    hits = list(parse_tabular(io.StringIO(TABULAR)))
    h = hits[0]
    assert h.length is None
    assert h.pident is None
    assert h.qcov is None
    assert h.scov is None
    assert h.qlen is None
    assert h.slen is None


def test_wide_and_narrow_hits_compare_equal_on_the_original_four_fields():
    # A caller that only ever looked at query_id/target_id/evalue/bitscore (every
    # caller today) must see identical behaviour whether the cache file is old or new.
    old = Hit('q1', 't1', 1e-20, 150.0)
    new = list(parse_tabular(io.StringIO(TABULAR_WIDE_DIAMOND)))[0]
    assert (old.query_id, old.target_id, old.evalue, old.bitscore) == \
           (new.query_id, new.target_id, new.evalue, new.bitscore)


def test_parse_phmmer_leaves_alignment_fields_none():
    # --tblout (per-sequence hits) carries no alignment geometry at all -- no length,
    # no coordinates, no identity. --domtblout would be needed for even coverage, and
    # HMMER reports no percent-identity in either format. Tracked as a follow-up
    # (issue #150), not attempted here: it needs a different output format entirely
    # (per-domain rows, not per-sequence), not just a wider row of the same format.
    h = list(parse_phmmer_tblout(io.StringIO(PHMMER_TBLOUT)))[0]
    assert h.length is None
    assert h.pident is None
    assert h.qlen is None
    assert h.slen is None
