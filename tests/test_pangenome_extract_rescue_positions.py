import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, GENOME_ONLY, ABSENT
from pangenome_extract_rescue_positions import (
    parse_tblastn_best_hit_positions,
    merge_best_hits,
    main,
)


def test_parse_tblastn_best_hit_positions_normalizes_minus_strand():
    lines = [
        # plus strand: sstart < send
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95",
        # minus strand: sstart > send -> start should be min(600,500)=500
        "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t600\t500\t1e-50\t150\t95",
    ]
    hits = parse_tblastn_best_hit_positions(lines, min_pident=90.0, min_qcov=80.0)
    # higher bitscore (200) wins
    assert hits[("famA", "s1")] == ("contig1", 500, 200.0)


def test_merge_best_hits_picks_highest_bitscore_across_files():
    per_file = [
        {("famA", "s1"): ("contig1", 500, 100.0)},
        {("famA", "s1"): ("contig2", 700, 250.0)},
    ]
    merged = merge_best_hits(per_file)
    assert merged[("famA", "s1")] == ("contig2", 700, 250.0)


def test_main_succeeds_with_zero_tblastn_tsv_args(monkeypatch, capsys):
    # Regression test: the workflow can feed a collected channel with
    # .ifEmpty([]) when there are zero GENOME_ONLY calls upstream, producing
    # zero --tblastn_tsv args. This must not be an argparse-required error;
    # it should just emit a header-only output.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1", "s2"])
        pm.set_call("famA", "s1", ABSENT)
        pm.to_tsv(matrix_file)

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_extract_rescue_positions.py",
            "--matrix", str(matrix_file),
            "--output", str(output_file),
        ])

        main()

        assert output_file.exists()
        content = output_file.read_text()
        assert content == "Short\tfamily\tcontig\tstart\n"
        captured = capsys.readouterr()
        assert "0 GENOME_ONLY positions written" in captured.err


def test_main_writes_position_for_genome_only_call(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        matrix_file = tmpdir / "matrix.tsv"
        pm = PresenceMatrix(families=["famA"], strains=["s1"])
        pm.set_call("famA", "s1", GENOME_ONLY)
        pm.to_tsv(matrix_file)

        tblastn_file = tmpdir / "tblastn.tsv"
        tblastn_file.write_text(
            "famA\ts1|contig1\t95.0\t100\t0\t0\t1\t100\t500\t600\t1e-50\t200\t95\n"
        )

        output_file = tmpdir / "output.tsv"

        monkeypatch.setattr(sys, "argv", [
            "pangenome_extract_rescue_positions.py",
            "--matrix", str(matrix_file),
            "--tblastn_tsv", str(tblastn_file),
            "--output", str(output_file),
        ])

        main()

        content = output_file.read_text()
        assert content == "Short\tfamily\tcontig\tstart\ns1\tfamA\tcontig1\t500\n"
