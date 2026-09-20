"""MERGE_PFAM_DOMTBLOUT's concatenation contract (issue #112).

FAMILY_PFAM_SCAN is scattered into one hmmscan per query chunk, and the
chunk outputs are reassembled with a plain `cat`. That is only correct if
parsing the concatenation gives exactly what parsing the pieces gives --
hmmscan writes a `#`-prefixed header AND footer into every domtblout, so a
naive concatenation interleaves those blocks in the middle of the merged
file. These tests pin that the parser is indifferent to it.

The Nextflow wiring itself (splitFasta fan-out, one SLURM job per chunk,
the shared hmmpress'd database) is not unit-testable here and is covered by
a real chunked pangenome smoke run instead.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_domain_enrichment import parse_domtblout  # noqa: E402

# Real hmmscan --domtblout shape: 22 whitespace-separated columns, with the
# target (Pfam) name first and the query (family rep) name in column 4.
HEADER = (
    "#                                                                       "
    "--- full sequence ---- --- best 1 domain ----\n"
    "# target name        accession   tlen query name  accession   qlen   "
    "E-value  score  bias\n"
)
FOOTER = (
    "#\n"
    "# Program:         hmmscan\n"
    "# Query file:      chunk.fa\n"
    "# [ok]\n"
)


def domtblout_line(domain: str, accession: str, family: str, ievalue: str) -> str:
    # columns: target(1) accession(2) tlen(3) query(4) accession(5) qlen(6)
    #          E-value(7) score(8) bias(9) #(10) of(11) c-Evalue(12)
    #          i-Evalue(13) ...
    cols = [domain, accession, "100", family, "-", "250",
            "1e-30", "100.0", "0.0", "1", "1", "1e-30", ievalue]
    cols += ["1", "100", "1", "100", "1", "100", "0.95", "-"]
    return " ".join(cols) + "\n"


CHUNK_A = HEADER + domtblout_line("NACHT", "PF05729.1", "fam1", "1e-25") + FOOTER
CHUNK_B = HEADER + domtblout_line("HET", "PF06985.1", "fam2", "1e-20") + FOOTER


def test_parsing_concatenated_chunks_equals_parsing_them_separately(tmp_path):
    a = tmp_path / "chunk_a.domtblout"
    b = tmp_path / "chunk_b.domtblout"
    a.write_text(CHUNK_A)
    b.write_text(CHUNK_B)
    merged = tmp_path / "pfam.domtblout"
    merged.write_text(CHUNK_A + CHUNK_B)

    assert parse_domtblout([str(merged)]) == parse_domtblout([str(a), str(b)])


def test_merged_domtblout_keeps_every_chunks_families(tmp_path):
    merged = tmp_path / "pfam.domtblout"
    merged.write_text(CHUNK_A + CHUNK_B)
    assert parse_domtblout([str(merged)]) == {"fam1": {"NACHT"}, "fam2": {"HET"}}


def test_interleaved_chunk_headers_are_not_parsed_as_hits(tmp_path):
    # The failure this guards against: a merged file's second header block
    # sits mid-file, where a parser that only skipped leading `#` lines
    # would read it as data and invent a family.
    merged = tmp_path / "pfam.domtblout"
    merged.write_text(CHUNK_A + CHUNK_B)
    families = parse_domtblout([str(merged)])
    assert not any(f.startswith("#") for f in families)
    assert len(families) == 2


def test_empty_chunk_contributes_nothing(tmp_path):
    # A chunk whose sequences all miss Pfam still produces a valid
    # header+footer-only domtblout; it must merge cleanly.
    merged = tmp_path / "pfam.domtblout"
    merged.write_text(CHUNK_A + HEADER + FOOTER + CHUNK_B)
    assert parse_domtblout([str(merged)]) == {"fam1": {"NACHT"}, "fam2": {"HET"}}
