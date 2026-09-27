"""lib/genbank_slice.py: rank rebuild, GFF3 CDS parsing, GenBank records."""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from Bio import SeqIO  # noqa: E402
from genbank_slice import build_record, cds_exons, cds_ids, rank_entries, safe_name  # noqa: E402
from pangenome_build_family_positions import build_family_positions  # noqa: E402
from pangenome_build_gene_positions import _resolve_ids  # noqa: E402

GENES = [
    ("S1", "p3", "c2", 50, 90), ("S1", "p1", "c1", 10, 40), ("S1", "p2", "c1", 60, 99),
    ("S2", "q1", "c1", 5, 30), ("S1", "p4", "c1", 200, 260), ("S1", "p5", "c10", 1, 9),
]
RESCUES = [("S1", "famR", "c1", 60), ("S1", "famZ", "c2", 10)]
MEMBER = {"S1|p1": "famA", "S1|p2": "famB", "S1|p3": "famC", "S1|p4": "famD", "S1|p5": "famE",
          "S2|q1": "famA"}


def test_rank_entries_match_build_family_positions():
    fp = build_family_positions(GENES, MEMBER, RESCUES)["S1"]
    rank_of = {(c, r): fam for fam, copies in fp.items() for c, r in copies}
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1"), ("S1", "c2"), ("S1", "c10")})
    for (strain, contig), entries in idx.items():
        for rank, pid, fam, _s, _e in entries:
            expected = rank_of[(contig, rank)]
            assert (MEMBER[f"S1|{pid}"] if pid else fam) == expected, (contig, rank)


def test_rank_entries_order_ties_gene_before_rescue():
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1")})
    assert [(e[1], e[2]) for e in idx[("S1", "c1")]] == [
        ("p1", None), ("p2", None), (None, "famR"), ("p4", None)]


def test_rank_offset_counts_contigs_that_sort_first():
    idx = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c2")})
    # c1 (4 entries) and c10 (1 entry) sort before c2.
    assert [e[0] for e in idx[("S1", "c2")]] == [5, 6]


def test_cds_ids_matches_gene_positions_rule():
    for attrs in ["ID=x;Parent=m1", "protein_id=P1;Parent=m1", "Parent=a,b", "ID=x"]:
        assert cds_ids(attrs) == _resolve_ids(attrs)


GFF = """##gff-version 3
c1\tsrc\tgene\t10\t40\t.\t+\t.\tID=g1
c1\tsrc\tmRNA\t10\t40\t.\t+\t.\tID=p1;Parent=g1
c1\tsrc\tCDS\t10\t20\t.\t+\t0\tID=cds1;Parent=p1
c1\tsrc\tCDS\t30\t40\t.\t+\t0\tID=cds1;Parent=p1
c1\tsrc\tCDS\t60\t99\t.\t-\t0\tID=cds2;Parent=p2
c2\tsrc\tCDS\t50\t90\t.\t+\t0\tID=cds3;Parent=p3
"""


def test_cds_exons_keeps_wanted_ids_on_wanted_contigs():
    ex = cds_exons(io.StringIO(GFF), {"c1"}, {"p1", "p2", "p3"})
    assert ex == {"p1": (1, [(10, 20), (30, 40)]), "p2": (-1, [(60, 99)])}


def test_build_record_writes_gene_and_cds_features(tmp_path):
    seq = "A" * 300
    entries = rank_entries(iter(GENES), iter(RESCUES), {("S1", "c1")})[("S1", "c1")][:3]
    exons = cds_exons(io.StringIO(GFF), {"c1"}, {"p1", "p2"})
    rec, summ = build_record("S1", "c1", seq, entries, exons, {"p1": "MK*", "p2": "MR"},
                             {"p1": "famA", "p2": "famB"}, lambda pid: pid)
    assert (summ["bp_start"], summ["bp_end"]) == (10, 99)
    assert (summ["n_genes"], summ["n_rescue"], summ["n_missing"]) == (2, 1, 0)
    assert summ["labels"] == {"p1": "famA", "p2": "famB"}
    out = tmp_path / "S1.gbk"
    SeqIO.write(rec, str(out), "genbank")
    back = SeqIO.read(str(out), "genbank")
    cds = [f for f in back.features if f.type == "CDS"]
    assert [f.qualifiers["locus_tag"][0] for f in cds] == ["p1", "p2"]
    assert cds[0].qualifiers["translation"][0] == "MK"
    assert cds[0].qualifiers["note"][0] == "family=famA"
    assert len(cds[0].location.parts) == 2 and cds[1].location.strand == -1
    assert len([f for f in back.features if f.type == "gene"]) == 2


def test_gene_without_a_cds_row_is_counted_missing():
    entries = [(0, "pX", None, 10, 40)]
    _rec, summ = build_record("S1", "c1", "A" * 50, entries, {}, {}, {}, lambda p: p)
    assert (summ["n_genes"], summ["n_missing"]) == (0, 1)


def test_safe_name():
    assert safe_name("B0858-Guatemala") == "B0858-Guatemala"
    assert safe_name("a/b c") == "a_b_c"
