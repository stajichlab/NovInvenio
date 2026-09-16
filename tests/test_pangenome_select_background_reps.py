import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_select_background_reps import select_background_families, write_background_fasta


def test_select_background_families_keeps_only_shell_and_cloud(tmp_path):
    freq = tmp_path / "frequency_table.tsv"
    freq.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t0.9\t90\tcore\n"
        "famB\t0.3\t30\tshell\n"
        "famC\t0.05\t5\tcloud\n"
        "famD\t0.01\t1\tsingleton\n"
    )
    result = select_background_families(str(freq))
    assert result == {"famB", "famC"}


def test_write_background_fasta_filters_by_family_set(tmp_path):
    rep_fasta = tmp_path / "tier1_rep_seq.fasta"
    rep_fasta.write_text(
        ">famA some description\nMSEQA\n"
        ">famB\nMSEQB\n"
        ">famC\nMSEQC\n"
    )
    out_path = tmp_path / "background_reps.fa"
    n = write_background_fasta(str(rep_fasta), {"famB", "famC"}, str(out_path))
    assert n == 2
    content = out_path.read_text()
    assert ">famB" in content and ">famC" in content and ">famA" not in content
