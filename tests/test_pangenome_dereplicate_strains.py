import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_dereplicate_strains import (
    compute_assembly_stats, parse_mash_dist, choose_representatives,
    groups_to_shorts, main,
)


def test_compute_assembly_stats_n50_and_contig_count(tmp_path):
    fa = tmp_path / "strain.dna.fa"
    # two contigs: 100bp and 300bp -> total 400, N50 is the 300bp contig
    fa.write_text(">contig1\n" + "A" * 100 + "\n>contig2\n" + "C" * 300 + "\n")
    stats = compute_assembly_stats(fa)
    assert stats["n_contigs"] == 2
    assert stats["total_length"] == 400
    assert stats["n50"] == 300


def test_parse_mash_dist_groups_below_threshold():
    lines = [
        "#query\ts1\ts2\ts3",
        "s1\t0.0000\t0.0005\t0.0500",
        "s2\t0.0005\t0.0000\t0.0520",
        "s3\t0.0500\t0.0520\t0.0000",
    ]
    groups = parse_mash_dist(lines, threshold=0.001)
    assert {"s1", "s2"} in groups
    assert {"s3"} in groups


def test_choose_representatives_picks_highest_n50():
    groups = [{"s1", "s2"}, {"s3"}]
    assembly_stats = {
        "s1": {"n50": 500_000, "n_contigs": 20, "total_length": 29_000_000},
        "s2": {"n50": 2_000_000, "n_contigs": 8, "total_length": 29_100_000},
        "s3": {"n50": 1_000_000, "n_contigs": 15, "total_length": 28_900_000},
    }
    reps = choose_representatives(groups, assembly_stats)
    assert reps == {"s1": "s2", "s2": "s2", "s3": "s3"}


def test_parse_mash_dist_empty_input():
    assert parse_mash_dist([], threshold=0.001) == []


def test_parse_mash_dist_transitivity():
    lines = [
        "#query\tA\tB\tC",
        "A\t0.0000\t0.0005\t0.0100",
        "B\t0.0005\t0.0000\t0.0005",
        "C\t0.0100\t0.0005\t0.0000",
    ]
    groups = parse_mash_dist(lines, threshold=0.001)
    assert len(groups) == 1
    assert {"A", "B", "C"} in groups


def test_groups_to_shorts_matches_by_basename():
    """Real bug fix: mash reports whatever path string it was invoked with.
    A shared upstream Nextflow MASH_SKETCH process (modules/pangenome/mash.nf)
    stages each strain's genome FASTA into its own work directory, so mash
    reports bare basenames there ('s1.fa'), not the
    '<data_dir>/dna/s1.fa' path this script's own dna_paths dict is built
    from -- path_to_short must be keyed (and looked up) by basename, or every
    real Nextflow-driven run would KeyError here."""
    path_groups = [
        {"s1.fa", "s2.fa"},   # mash's reported names: bare basenames (Nextflow-staged)
        {"s3.fa"},
    ]
    path_to_short = {"s1.fa": "s1", "s2.fa": "s2", "s3.fa": "s3"}
    short_groups = groups_to_shorts(path_groups, path_to_short)
    assert {"s1", "s2"} in short_groups
    assert {"s3"} in short_groups
    assert len(short_groups) == 2


def test_groups_to_shorts_matches_full_paths_via_basename_too():
    """Standalone-CLI fallback (no shared MASH_SKETCH): mash reports the full
    path it was actually given. Basename matching still works since it's a
    superset of the full-path case (as long as filenames are unique)."""
    path_groups = [{"/path/to/data/dna/s1.fa", "/path/to/data/dna/s2.fa"}]
    path_to_short = {"s1.fa": "s1", "s2.fa": "s2"}
    short_groups = groups_to_shorts(path_groups, path_to_short)
    assert {"s1", "s2"} in short_groups


def test_groups_to_shorts_raises_on_unrecognized_path():
    with pytest.raises(KeyError, match="not one of the"):
        groups_to_shorts([{"unknown.fa"}], {"s1.fa": "s1"})


def test_main_raises_on_duplicate_dna_basenames(tmp_path, monkeypatch):
    # Regression test: two different Short IDs pointing at DNA files that
    # share a basename (e.g. copy-paste error in the samplesheet) must
    # raise loudly, not silently drop one strain from dedup grouping (the
    # path_to_short dict build collapses on the collision).
    dna_dir = tmp_path / "data" / "dna"
    dna_dir.mkdir(parents=True)
    (dna_dir / "strainA.fa").write_text(">c1\nACGT\n")
    dup_dir = tmp_path / "data2" / "dna"
    # Same basename "strainA.fa" resolved from a *different* dna filename
    # column value isn't directly reachable via this script's fixed
    # "<data_dir>/dna/<DNA>" convention, so instead simulate the collision
    # the way it actually happens: two Short rows list the SAME DNA filename.
    config_csv = tmp_path / "config.csv"
    config_csv.write_text(
        "GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup\n"
        "IN,Sp1,st1,p1.fa,strainA.fa,s1,t1\n"
        "IN,Sp2,st2,p2.fa,strainA.fa,s2,t1\n"
    )
    dist_tsv = tmp_path / "dist.tsv"
    dist_tsv.write_text("#query\tstrainA.fa\nstrainA.fa\t0.0000\n")

    monkeypatch.setattr(sys, "argv", [
        "pangenome_dereplicate_strains.py",
        "--config", str(config_csv),
        "--data_dir", str(tmp_path / "data"),
        "--mash_dist_tsv", str(dist_tsv),
        "--output", str(tmp_path / "out.tsv"),
    ])

    with pytest.raises(ValueError, match="Duplicate DNA basenames"):
        main()
