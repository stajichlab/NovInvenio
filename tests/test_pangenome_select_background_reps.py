import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

import pangenome_select_background_reps
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


def test_main_warns_when_nonempty_background_writes_zero_records(tmp_path, monkeypatch):
    # frequency_table has a real shell/cloud background, but the rep FASTA's
    # headers use a different ID convention entirely (simulating an upstream
    # header-format drift) -> write_background_fasta writes 0 records even
    # though the background family set was non-empty. This is a real join
    # failure, not a legitimately-empty background, and must warn loudly.
    freq = tmp_path / "frequency_table.tsv"
    freq.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famB\t0.3\t30\tshell\n"
    )
    rep_fasta = tmp_path / "tier1_rep_seq.fasta"
    rep_fasta.write_text(">totally_different_id\nMSEQA\n")
    out_path = tmp_path / "background_reps.fa"

    argv = [
        "pangenome_select_background_reps.py",
        "--rep_fasta", str(rep_fasta),
        "--frequency_table", str(freq),
        "--output", str(out_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pangenome_select_background_reps.main()
    assert any("0 FASTA records were written" in str(w.message) for w in caught)


def test_main_no_warning_when_background_legitimately_empty(tmp_path, monkeypatch):
    # No shell/cloud families at all (e.g. an all-core frequency table) ->
    # 0 records written is the correct, expected outcome; must not warn.
    freq = tmp_path / "frequency_table.tsv"
    freq.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t0.9\t90\tcore\n"
    )
    rep_fasta = tmp_path / "tier1_rep_seq.fasta"
    rep_fasta.write_text(">famA\nMSEQA\n")
    out_path = tmp_path / "background_reps.fa"

    argv = [
        "pangenome_select_background_reps.py",
        "--rep_fasta", str(rep_fasta),
        "--frequency_table", str(freq),
        "--output", str(out_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pangenome_select_background_reps.main()
    assert not any("0 FASTA records were written" in str(w.message) for w in caught)
