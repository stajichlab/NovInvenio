"""CLI tests for bin/pangenome_island_gbk_slice.py on a one-strain fixture."""
import csv
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from Bio import SeqIO  # noqa: E402
from pangenome_island_gbk_slice import main  # noqa: E402

GFF = """##gff-version 3
c1\tsrc\tCDS\t10\t40\t.\t+\t0\tParent=p1
c1\tsrc\tCDS\t60\t99\t.\t-\t0\tParent=p2
c1\tsrc\tCDS\t120\t150\t.\t+\t0\tParent=p3
"""
REGION_HEADER = "locus_key\tstrain\treason\trow_class\tspecies\tcontig\trank_lo\trank_hi\tanchors\n"


def fixture(d: Path, anchors="0:famA;1:famB") -> list[str]:
    for sub in ("dna", "pep", "gff3"):
        (d / "data" / sub).mkdir(parents=True, exist_ok=True)
    (d / "data" / "dna" / "S1.dna.fa").write_text(">c1\n" + "ACGT" * 50 + "\n>c2\nAAAA\n")
    (d / "data" / "pep" / "S1.pep.fa").write_text(">p1\nMKV*\n>p2\nMRR\n>p3\nMQQ\n")
    (d / "data" / "gff3" / "S1.gff3").write_text(GFF)
    (d / "config.csv").write_text("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
                                  "IN,Sp,S1,S1.pep.fa,S1.dna.fa,S1.gff3,S1,t\n")
    (d / "gene_positions.tsv").write_text("Short\tprotein_id\tcontig\tstart\tend\n"
                                          "S1\tp1\tc1\t10\t40\nS1\tp2\tc1\t60\t99\n"
                                          "S1\tp3\tc1\t120\t150\n")
    (d / "rescue.tsv").write_text("Short\tfamily\tcontig\tstart\nS1\tfamR\tc1\t110\n")
    (d / "cluster.tsv").write_text("famA\tS1|p1\nfamB\tS1|p2\nfamC\tS1|p3\n")
    (d / "regions.tsv").write_text(REGION_HEADER +
                                   f"L001\tS1\texemplar\tfull\tSp\tc1\t0\t3\t{anchors}\n")
    return ["--regions", str(d / "regions.tsv"), "--config", str(d / "config.csv"),
            "--data_dir", str(d / "data"), "--gff3_dir", str(d / "data" / "gff3"),
            "--gene_positions", str(d / "gene_positions.tsv"),
            "--rescue_positions", str(d / "rescue.tsv"),
            "--cluster_tsv", str(d / "cluster.tsv"), "--out_dir", str(d / "gbk")]


def test_writes_one_genbank_per_region_with_families(tmp_path):
    assert main(fixture(tmp_path)) == 0
    rec = SeqIO.read(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank")
    cds = [f for f in rec.features if f.type == "CDS"]
    assert [f.qualifiers["note"][0] for f in cds] == ["family=famA", "family=famB", "family=famC"]
    assert cds[0].qualifiers["translation"][0] == "MKV"
    assert len(rec.seq) == 150 - 10 + 1


def test_groups_csv_maps_locus_tags_to_families(tmp_path):
    main(fixture(tmp_path))
    rows = list(csv.reader(open(tmp_path / "gbk" / "L001" / "groups.csv")))
    assert rows == [["p1", "famA"], ["p2", "famB"], ["p3", "famC"]]


def test_slices_table_reports_bp_and_counts(tmp_path):
    main(fixture(tmp_path))
    (row,) = list(csv.DictReader(open(tmp_path / "gbk" / "island_slices.tsv"), delimiter="\t"))
    assert (row["bp_start"], row["bp_end"], row["n_genes"], row["n_rescue"]) == ("10", "150", "3", "1")
    # C1: below the (default 20000) max_gap, one block, no gap skipped.
    assert (row["n_blocks"], row["gap_bp"], row["max_gap_bp"]) == ("1", "0", "0")
    assert row["drawn_bp"] == str(150 - 10 + 1)


def test_max_gap_splits_the_genbank_record_into_blocks(tmp_path):
    # Entries in this fixture (rank_lo=0, rank_hi=3) are p1(10-40), p2(60-99),
    # the rescue hit famR(c1:110, no gene model) and p3(120-150). --max_gap
    # 10 splits at the p1->p2 gap (20 bp) and the p2->famR gap (11 bp), but
    # not famR->p3 (10 bp, not > 10): three blocks, one .gbk record each.
    args = fixture(tmp_path) + ["--max_gap", "10"]
    assert main(args) == 0
    recs = list(SeqIO.parse(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank"))
    assert [r.id for r in recs] == ["S1_b0", "S1_b1", "S1_b2"]
    assert [len([f for f in r.features if f.type == "CDS"]) for r in recs] == [1, 1, 1]
    (row,) = list(csv.DictReader(open(tmp_path / "gbk" / "island_slices.tsv"), delimiter="\t"))
    assert row["n_blocks"] == "3"
    assert row["gap_bp"] == str((60 - 40) + (110 - 99))
    assert row["max_gap_bp"] == str(60 - 40)
    drawn = (40 - 10 + 1) + (99 - 60 + 1) + (150 - 110 + 1)
    assert row["drawn_bp"] == str(drawn)
    # groups.csv still maps every locus_tag across all blocks to its family.
    rows = list(csv.reader(open(tmp_path / "gbk" / "L001" / "groups.csv")))
    assert rows == [["p1", "famA"], ["p2", "famB"], ["p3", "famC"]]


def test_max_gap_zero_never_splits(tmp_path):
    args = fixture(tmp_path) + ["--max_gap", "0"]
    assert main(args) == 0
    recs = list(SeqIO.parse(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank"))
    assert len(recs) == 1 and recs[0].id == "S1"


def test_anchor_mismatch_fails(tmp_path):
    assert main(fixture(tmp_path, anchors="0:famB")) == 1


def test_strain_without_files_is_skipped(tmp_path):
    args = fixture(tmp_path)
    (tmp_path / "data" / "dna" / "S1.dna.fa").unlink()
    assert main(args) == 0
    lines = (tmp_path / "gbk" / "island_slices.tsv").read_text().splitlines()
    assert len(lines) == 1


def test_no_regions_writes_a_header_only_table(tmp_path):
    args = fixture(tmp_path)
    (tmp_path / "regions.tsv").write_text(REGION_HEADER)
    assert main(args) == 0
    assert (tmp_path / "gbk" / "island_slices.tsv").read_text().startswith("locus_key\t")


def test_protein_id_clash_across_strains_prefixes_labels(tmp_path):
    args = fixture(tmp_path)
    d = tmp_path
    for sub, ext in {"dna": "dna.fa", "pep": "pep.fa", "gff3": "gff3"}.items():
        text = (d / "data" / sub / f"S1.{ext}").read_text()
        (d / "data" / sub / f"S2.{ext}").write_text(text)
    with open(d / "config.csv", "a") as fh:
        fh.write("IN,Sp,S2,S2.pep.fa,S2.dna.fa,S2.gff3,S2,t\n")
    with open(d / "gene_positions.tsv", "a") as fh:
        fh.write("S2\tp1\tc1\t10\t40\nS2\tp2\tc1\t60\t99\nS2\tp3\tc1\t120\t150\n")
    with open(d / "cluster.tsv", "a") as fh:
        fh.write("famA\tS2|p1\nfamB\tS2|p2\nfamC\tS2|p3\n")
    with open(d / "regions.tsv", "a") as fh:
        fh.write("L001\tS2\tfill\tfull\tSp\tc1\t0\t2\t0:famA\n")
    assert main(args) == 0
    labels = [r[0] for r in csv.reader(open(d / "gbk" / "L001" / "groups.csv"))]
    assert "S1|p1" in labels and "S2|p1" in labels


def test_gzipped_genome_and_proteins_are_read(tmp_path):
    # Review Focus 4: study FASTAs may be .gz; lib/compressed_io handles them.
    import gzip
    args = fixture(tmp_path)
    for sub, name in (("dna", "S1.dna.fa"), ("pep", "S1.pep.fa")):
        src = tmp_path / "data" / sub / name
        with gzip.open(str(src) + ".gz", "wt") as fh:
            fh.write(src.read_text())
        src.unlink()
    cfg = tmp_path / "config.csv"
    cfg.write_text(cfg.read_text().replace("S1.pep.fa,S1.dna.fa", "S1.pep.fa.gz,S1.dna.fa.gz"))
    assert main(args) == 0
    rec = SeqIO.read(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank")
    assert len([f for f in rec.features if f.type == "CDS"]) == 3
