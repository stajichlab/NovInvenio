import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_prefix_fasta import prefix_fasta


def test_prefix_fasta_rewrites_first_token_only(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">protA description here\nMKV\n>protB\nMKL\n")
    out = tmp_path / "out.fa"

    n = prefix_fasta(str(src), str(out), short="AF293", id_sep="|")

    assert n == 2
    text = out.read_text()
    assert text.splitlines() == [
        ">AF293|protA description here",
        "MKV",
        ">AF293|protB",
        "MKL",
    ]


def test_prefix_fasta_custom_separator(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">g1\nAAA\n")
    out = tmp_path / "out.fa"

    prefix_fasta(str(src), str(out), short="S1", id_sep="__")

    assert out.read_text().splitlines()[0] == ">S1__g1"


def test_prefix_fasta_preserves_sequence_lines_untouched(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">g1\nAAAA\nCCCC\n>g2\nGGGG\n")
    out = tmp_path / "out.fa"

    prefix_fasta(str(src), str(out), short="S1")

    lines = out.read_text().splitlines()
    assert lines == [">S1|g1", "AAAA", "CCCC", ">S1|g2", "GGGG"]
