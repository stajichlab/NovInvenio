import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

import pangenome_domain_enrichment
from pangenome_domain_enrichment import parse_domtblout, parse_domtblout_accessions, domain_enrichment, bare_pfam_accession


def test_bare_pfam_accession_strips_version_suffix():
    assert bare_pfam_accession("PF13577.9") == "PF13577"
    assert bare_pfam_accession("PF00001.1") == "PF00001"


def test_bare_pfam_accession_passes_through_sentinel():
    assert bare_pfam_accession("-") == "-"


def test_parse_domtblout_filters_by_ievalue(tmp_path):
    dtbl = tmp_path / "test.domtblout"
    # hmmscan domtblout columns (0-indexed): 0 target_name, 1 target_acc,
    # 2 tlen, 3 query_name, 4 query_acc, 5 qlen, 6 full_evalue, 7 full_score,
    # 8 full_bias, 9 dom#, 10 dom_of, 11 c-Evalue, 12 i-Evalue, 13 dom_score,
    # 14 dom_bias, 15-20 hmm/ali/env coords, 21 acc. parse_domtblout reads
    # parts[12] -- the i-Evalue must actually sit there, or this test can't
    # distinguish pass/fail (caught during Opus review: an earlier draft of
    # this fixture put the discriminating value at index 9/10 instead).
    dtbl.write_text(
        "# comment line\n"
        "PF00001 - 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 1.0e-10 100.0 20.0 10 20 10 20 10 20 0.99\n"
        "PF00002 - 100 famA - 50 5.0e-02 100.0 20.0 1 1 5.0e-02 5.0e-02 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )
    hits = parse_domtblout([str(dtbl)], max_ievalue=1e-3)
    assert hits == {"famA": {"PF00001"}}


def test_parse_domtblout_warns_when_nonempty_input_yields_zero_hits(tmp_path):
    # A domtblout with real data lines whose format doesn't match what this
    # parser expects (e.g. an upstream column-layout drift) -> 0 hits parsed
    # out of nonempty input. This should never fail silently.
    dtbl = tmp_path / "test.domtblout"
    dtbl.write_text("this is not a real domtblout line at all\n")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        hits = parse_domtblout([str(dtbl)], max_ievalue=1e-3)
    assert hits == {}
    assert any("0 family/domain hits" in str(w.message) for w in caught)


def test_parse_domtblout_no_warning_on_genuinely_empty_file(tmp_path):
    dtbl = tmp_path / "test.domtblout"
    dtbl.write_text("")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        hits = parse_domtblout([str(dtbl)], max_ievalue=1e-3)
    assert hits == {}
    assert not any("0 family/domain hits" in str(w.message) for w in caught)


def test_parse_domtblout_accessions_maps_domain_name_to_accession(tmp_path):
    dtbl = tmp_path / "test.domtblout"
    dtbl.write_text(
        "SnoaL_2 PF13577.9 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 1.0e-10 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )
    accessions = parse_domtblout_accessions([str(dtbl)], max_ievalue=1e-3)
    assert accessions == {"SnoaL_2": "PF13577.9"}


def test_cli_output_includes_pfam_accession_column(tmp_path, monkeypatch):
    dtbl = tmp_path / "test.domtblout"
    dtbl.write_text(
        "SnoaL_2 PF13577.9 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 1.0e-10 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text(
        "family\tfrequency\tbin\n"
        "famA\t0.5\tshell\n"
        "famB\t0.5\tshell\n"
    )
    significant_islands = tmp_path / "significant_islands.tsv"
    significant_islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\n"
        "2\ts1\t1\tfamA\t1\ttrans\n"
    )
    output = tmp_path / "island_pfam_enrichment.tsv"
    argv = [
        "pangenome_domain_enrichment.py",
        "--significant_islands", str(significant_islands),
        "--domtblout", str(dtbl),
        "--frequency_table", str(frequency_table),
        "--output", str(output),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pangenome_domain_enrichment.main()
    lines = output.read_text().splitlines()
    header = lines[0].split("\t")
    assert header == [
        "domain", "pfam_accession", "pfam_url", "n_with_domain_in_islands", "n_with_domain_in_background",
        "n_island_families", "n_background_families", "fisher_p", "fdr_q",
    ]
    data_row = lines[1].split("\t")
    assert data_row[0] == "SnoaL_2"
    assert data_row[1] == "PF13577.9"
    assert data_row[2] == "https://www.ebi.ac.uk/interpro/entry/pfam/PF13577/"


def test_domain_enrichment_fisher_and_bh(tmp_path):
    # famA (island) has domainX; 3 background families total (famA, famB, famC),
    # only famA has domainX -> strong enrichment signal.
    island_members = {"famA"}
    background = {"famA", "famB", "famC"}
    family_domains = {"famA": {"domainX"}}
    result = domain_enrichment(island_members, background, family_domains)
    assert len(result) == 1
    row = result[0]
    assert row["domain"] == "domainX"
    assert row["n_with_domain_in_islands"] == 1
    assert row["n_with_domain_in_background"] == 1
    assert "fdr_q" in row


def test_domain_enrichment_excludes_island_members_outside_background(capsys):
    # famZ is an island member but a singleton (not in shell+cloud background)
    island_members = {"famA", "famZ"}
    background = {"famA", "famB"}
    family_domains = {"famA": {"domainX"}}
    result = domain_enrichment(island_members, background, family_domains)
    captured = capsys.readouterr()
    assert "outside the eligible" in captured.err
    # famZ excluded -> n_island_families should reflect only famA
    assert all(row["n_island_families"] == 1 for row in result)


def test_cli_domain_evalue_threads_through_to_parse_domtblout(tmp_path, monkeypatch):
    # famA has a domain hit at i-Evalue=5.0e-02 -- filtered out under the
    # default (1e-3) --domain_evalue, but retained under a looser one. This
    # confirms --domain_evalue is actually wired to parse_domtblout's
    # max_ievalue, not just accepted and ignored (the bug this test guards
    # against: the CLI arg previously did not exist at all, and
    # parse_domtblout's max_ievalue was hardcoded to 1e-3 unconditionally).
    dtbl = tmp_path / "test.domtblout"
    dtbl.write_text(
        "PF00001 - 100 famA - 50 5.0e-02 100.0 20.0 1 1 5.0e-02 5.0e-02 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )
    frequency_table = tmp_path / "frequency_table.tsv"
    frequency_table.write_text(
        "family\tfrequency\tbin\n"
        "famA\t0.5\tshell\n"
        "famB\t0.5\tshell\n"
    )
    significant_islands = tmp_path / "significant_islands.tsv"
    significant_islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\n"
        "2\ts1\t1\tfamA\t1\ttrans\n"
    )

    def run(domain_evalue_args):
        output = tmp_path / f"out_{'_'.join(domain_evalue_args) or 'default'}.tsv"
        argv = [
            "pangenome_domain_enrichment.py",
            "--significant_islands", str(significant_islands),
            "--domtblout", str(dtbl),
            "--frequency_table", str(frequency_table),
            "--output", str(output),
            *domain_evalue_args,
        ]
        monkeypatch.setattr(sys, "argv", argv)
        pangenome_domain_enrichment.main()
        return output.read_text()

    default_output = run([])
    assert "PF00001" not in default_output  # filtered out at default 1e-3

    loose_output = run(["--domain_evalue", "0.1"])
    assert "PF00001" in loose_output  # retained once threshold is loosened
