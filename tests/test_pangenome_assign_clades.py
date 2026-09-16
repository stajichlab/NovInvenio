import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_assign_clades import parse_mash_dist_matrix, main


def test_parse_mash_dist_matrix_symmetrizes_and_zeros_diagonal():
    lines = [
        "#query\tA\tB",
        "A\t0.0000\t0.0050",
        "B\t0.0052\t0.0000",
    ]
    names, mat = parse_mash_dist_matrix(lines)
    assert names == ["A", "B"]
    assert mat[0, 0] == 0.0 and mat[1, 1] == 0.0
    assert mat[0, 1] == mat[1, 0]


def test_main_raises_on_duplicate_dna_basenames(tmp_path, monkeypatch):
    # Regression test: two different Short IDs whose samplesheet DNA column
    # points at the same filename collide in path_to_short, which would
    # otherwise silently drop a strain from clade assignment instead of
    # erroring, even though the surrounding code assumes basenames are
    # unique per this script's module docstring.
    dna_dir = tmp_path / "data" / "dna"
    dna_dir.mkdir(parents=True)
    (dna_dir / "strainA.fa").write_text(">c1\nACGT\n")

    config_csv = tmp_path / "config.csv"
    config_csv.write_text(
        "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
        "IN,Sp1,st1,p1.fa,strainA.fa,s1,t1\n"
        "IN,Sp2,st2,p2.fa,strainA.fa,s2,t1\n"
    )
    dist_tsv = tmp_path / "dist.tsv"
    dist_tsv.write_text(
        "#query\tstrainA.fa\tstrainA.fa\n"
        "strainA.fa\t0.0000\t0.0000\n"
        "strainA.fa\t0.0000\t0.0000\n"
    )

    monkeypatch.setattr(sys, "argv", [
        "pangenome_assign_clades.py",
        "--config", str(config_csv),
        "--data_dir", str(tmp_path / "data"),
        "--groups", "IN",
        "--mash_dist_tsv", str(dist_tsv),
        "--output", str(tmp_path / "out.tsv"),
    ])

    with pytest.raises(ValueError, match="Duplicate DNA basenames"):
        main()
