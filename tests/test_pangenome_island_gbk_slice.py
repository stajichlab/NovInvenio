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
    # the rescue hit famR(c1:110, no gene model) and p3(120-150). M2: split
    # points come from GENE entries only, so --max_gap 10 splits at the
    # p1->p2 gene gap (20 bp) and the p2->p3 gene gap (21 bp); famR never
    # participates in the split decision, and is attached to the nearest
    # resulting block (p3's, 10 bp away, vs. 11 bp to p2's) instead of
    # forming a block of its own -- still three blocks, one .gbk record
    # each, same grouping as before M2 (only the reported gap stats change).
    args = fixture(tmp_path) + ["--max_gap", "10"]
    assert main(args) == 0
    recs = list(SeqIO.parse(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank"))
    # Review fix round 1, item 1: block records are named "b1", "b2", ...
    # (unique within the file), never "<strain>_bN" -- that truncated to the
    # same 16 chars for every block of a strain name >= 16 chars.
    assert [r.id for r in recs] == ["b1", "b2", "b3"]
    assert [len([f for f in r.features if f.type == "CDS"]) for r in recs] == [1, 1, 1]
    (row,) = list(csv.DictReader(open(tmp_path / "gbk" / "island_slices.tsv"), delimiter="\t"))
    assert row["n_blocks"] == "3"
    assert row["gap_bp"] == str((60 - 40) + (120 - 99))
    assert row["max_gap_bp"] == str(120 - 99)
    drawn = (40 - 10 + 1) + (99 - 60 + 1) + (150 - 110 + 1)
    assert row["drawn_bp"] == str(drawn)
    # groups.csv still maps every locus_tag across all blocks to its family.
    rows = list(csv.reader(open(tmp_path / "gbk" / "L001" / "groups.csv")))
    assert rows == [["p1", "famA"], ["p2", "famB"], ["p3", "famC"]]


def test_max_gap_zero_never_splits(tmp_path):
    args = fixture(tmp_path) + ["--max_gap", "0"]
    assert main(args) == 0
    recs = list(SeqIO.parse(str(tmp_path / "gbk" / "L001" / "S1.gbk"), "genbank"))
    # One rule applied consistently (review fix round 1, item 1): even a
    # single, unsplit block is named "b1", not the strain name.
    assert len(recs) == 1 and recs[0].id == "b1"


def test_long_strain_names_produce_unique_block_record_names(tmp_path):
    # Review fix round 1, item 1, real evidence: strain names like
    # "578-1_L_NEW_CPA0049" (20 chars) used to truncate every block's
    # record name to the same 16 characters, since the old rule appended
    # "_bN" to the (possibly long) strain name and then truncated the
    # whole thing. Two such strains, each split into 3 blocks by a tight
    # --max_gap, must each get 3 distinct record names.
    d = tmp_path
    for sub in ("dna", "pep", "gff3"):
        (d / "data" / sub).mkdir(parents=True, exist_ok=True)
    long_names = ["578-1_L_NEW_CPA0049", "574-0_S_OLD_CPA0039"]
    cfg_lines = ["GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"]
    gene_lines = ["Short\tprotein_id\tcontig\tstart\tend\n"]
    cluster_lines = []
    region_lines = [REGION_HEADER]
    for strain in long_names:
        assert len(strain) >= 16
        (d / "data" / "dna" / f"{strain}.dna.fa").write_text(f">c1\n{'ACGT' * 100}\n")
        (d / "data" / "pep" / f"{strain}.pep.fa").write_text(">p1\nMKV*\n>p2\nMRR\n>p3\nMQQ\n")
        (d / "data" / "gff3" / f"{strain}.gff3").write_text(
            "##gff-version 3\nc1\tsrc\tCDS\t10\t40\t.\t+\t0\tParent=p1\n"
            "c1\tsrc\tCDS\t400\t440\t.\t-\t0\tParent=p2\nc1\tsrc\tCDS\t900\t950\t.\t+\t0\tParent=p3\n")
        cfg_lines.append(f"IN,Sp,{strain},{strain}.pep.fa,{strain}.dna.fa,{strain}.gff3,{strain},t\n")
        gene_lines.append(f"{strain}\tp1\tc1\t10\t40\n{strain}\tp2\tc1\t400\t440\n"
                          f"{strain}\tp3\tc1\t900\t950\n")
        cluster_lines.append(f"famA\t{strain}|p1\nfamB\t{strain}|p2\nfamC\t{strain}|p3\n")
        region_lines.append(f"L001\t{strain}\texemplar\tfull\tSp\tc1\t0\t2\t0:famA\n")
    (d / "config.csv").write_text("".join(cfg_lines))
    (d / "gene_positions.tsv").write_text("".join(gene_lines))
    (d / "cluster.tsv").write_text("".join(cluster_lines))
    (d / "regions.tsv").write_text("".join(region_lines))
    args = ["--regions", str(d / "regions.tsv"), "--config", str(d / "config.csv"),
           "--data_dir", str(d / "data"), "--gff3_dir", str(d / "data" / "gff3"),
           "--gene_positions", str(d / "gene_positions.tsv"),
           "--cluster_tsv", str(d / "cluster.tsv"), "--out_dir", str(d / "gbk"),
           "--max_gap", "100"]
    assert main(args) == 0
    for strain in long_names:
        recs = list(SeqIO.parse(str(d / "gbk" / "L001" / f"{strain}.gbk"), "genbank"))
        assert [r.id for r in recs] == ["b1", "b2", "b3"]
        assert len({r.id for r in recs}) == 3


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


# ---- M7 (final review): clinker_order.txt captures pick order, not the
# by_strain loop's alphabetical processing order. ----

def test_clinker_order_file_reflects_pick_order_not_alphabetical(tmp_path):
    d = tmp_path
    for sub in ("dna", "pep", "gff3"):
        (d / "data" / sub).mkdir(parents=True, exist_ok=True)
    cfg_lines = ["GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"]
    gene_lines = ["Short\tprotein_id\tcontig\tstart\tend\n"]
    cluster_lines = []
    # island_regions.tsv row order is the pick order: S2 (the exemplar)
    # first, S1 (a later "fill" pick) second -- the REVERSE of the
    # by_strain loop's alphabetical processing order.
    region_lines = [REGION_HEADER, "L001\tS2\texemplar\tfull\tSp\tc1\t0\t2\t0:famA\n",
                    "L001\tS1\tfill\tfull\tSp\tc1\t0\t2\t0:famA\n"]
    for strain in ("S1", "S2"):
        (d / "data" / "dna" / f"{strain}.dna.fa").write_text(f">c1\n{'ACGT' * 100}\n")
        (d / "data" / "pep" / f"{strain}.pep.fa").write_text(">p1\nMKV*\n>p2\nMRR\n>p3\nMQQ\n")
        (d / "data" / "gff3" / f"{strain}.gff3").write_text(
            "##gff-version 3\nc1\tsrc\tCDS\t10\t40\t.\t+\t0\tParent=p1\n"
            "c1\tsrc\tCDS\t60\t99\t.\t-\t0\tParent=p2\nc1\tsrc\tCDS\t120\t150\t.\t+\t0\tParent=p3\n")
        cfg_lines.append(f"IN,Sp,{strain},{strain}.pep.fa,{strain}.dna.fa,{strain}.gff3,{strain},t\n")
        gene_lines.append(f"{strain}\tp1\tc1\t10\t40\n{strain}\tp2\tc1\t60\t99\n"
                          f"{strain}\tp3\tc1\t120\t150\n")
        cluster_lines.append(f"famA\t{strain}|p1\nfamB\t{strain}|p2\nfamC\t{strain}|p3\n")
    (d / "config.csv").write_text("".join(cfg_lines))
    (d / "gene_positions.tsv").write_text("".join(gene_lines))
    (d / "cluster.tsv").write_text("".join(cluster_lines))
    (d / "regions.tsv").write_text("".join(region_lines))
    args = ["--regions", str(d / "regions.tsv"), "--config", str(d / "config.csv"),
           "--data_dir", str(d / "data"), "--gff3_dir", str(d / "data" / "gff3"),
           "--gene_positions", str(d / "gene_positions.tsv"),
           "--cluster_tsv", str(d / "cluster.tsv"), "--out_dir", str(d / "gbk")]
    assert main(args) == 0
    order = (d / "gbk" / "L001" / "clinker_order.txt").read_text().splitlines()
    assert order == ["S2.gbk", "S1.gbk"]


def test_clinker_order_file_omits_strains_that_were_skipped(tmp_path):
    # A picked strain missing its DNA/protein/GFF3 file is skipped (an
    # existing WARNING path); clinker_order.txt must not list a .gbk that
    # was never written. S2 gets gene_positions/cluster_tsv rows (so the
    # rank rebuild still matches island_regions.tsv's anchors) but no
    # actual DNA/protein/GFF3 file under data_dir, so it hits the "file not
    # found" skip.
    args = fixture(tmp_path)
    d = tmp_path
    with open(d / "config.csv", "a") as fh:
        fh.write("IN,Sp,S2,S2.pep.fa,S2.dna.fa,S2.gff3,S2,t\n")
    with open(d / "gene_positions.tsv", "a") as fh:
        fh.write("S2\tp1\tc1\t10\t40\nS2\tp2\tc1\t60\t99\nS2\tp3\tc1\t120\t150\n")
    with open(d / "cluster.tsv", "a") as fh:
        fh.write("famA\tS2|p1\nfamB\tS2|p2\nfamC\tS2|p3\n")
    with open(d / "regions.tsv", "a") as fh:
        fh.write("L001\tS2\tfill\tfull\tSp\tc1\t0\t2\t0:famA\n")
    assert main(args) == 0
    order = (d / "gbk" / "L001" / "clinker_order.txt").read_text().splitlines()
    assert order == ["S1.gbk"]


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
