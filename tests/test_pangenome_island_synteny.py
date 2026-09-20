import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).parent.parent / "bin"


def write_inputs(tmp_path):
    islands = tmp_path / "islands_with_domains.tsv"
    islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
        "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
        "n_contigs_in_locus\n"
        "2\tS1\t2\tfamA,famB\t1\tunexplained_physical\tNACHT\t"
        "S1:c1:100-400\tc1\t100\t400\t2\t1\n"
        "1\tS2\t9\tfamC\t1\tunexplained_physical\t-\t"
        "S2:c2:1-900\tc2\t1\t900\t1\t1\n"
    )
    matrix = tmp_path / "presence_matrix.tsv"
    matrix.write_text(
        "family\tS1\tS2\nfamA\tpresent\tpresent\nfamB\tpresent\tabsent\n"
        "famC\tabsent\tpresent\n"
    )
    positions = tmp_path / "family_positions.tsv"
    positions.write_text("Short\tfamily\tcontig\trank\nS1\tfamA\tc1\t1\n"
                         "S1\tfamB\tc1\t2\n")
    return islands, matrix, positions


def run_cli(tmp_path, *extra):
    islands, matrix, positions = write_inputs(tmp_path)
    out = tmp_path / "island_synteny.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out), *extra],
        check=True,
    )
    return out


def payload_of(html_path):
    html = html_path.read_text()
    start = html.index('<script type="application/json"')
    start = html.index(">", start) + 1
    return json.loads(html[start:html.index("</script>", start)])


def test_cli_writes_a_self_contained_page(tmp_path):
    out = run_cli(tmp_path)
    html = out.read_text()
    assert html.lstrip().lower().startswith("<!doctype html")
    assert "demo" in html


def test_cli_excludes_the_single_strain_island(tmp_path):
    payload = payload_of(run_cli(tmp_path))
    assert [i["locus_id"] for i in payload["islands"]] == ["S1:c1:100-400"]
    assert payload["n_islands_excluded"] == 1


def test_cli_orders_columns_by_locus_rank(tmp_path):
    payload = payload_of(run_cli(tmp_path))
    assert payload["islands"][0]["families"] == ["famA", "famB"]


def test_cli_min_strains_is_configurable(tmp_path):
    payload = payload_of(run_cli(tmp_path, "--min_strains", "3"))
    assert payload["islands"] == []


def test_cli_stderr_summary_reports_drawn_excluded_and_truncated(tmp_path):
    """The stderr summary must report all three counts (Ruling R12).

    A summary that reports only "drawn" and "excluded" under-reports on a
    real run where most qualifying islands are cut by --top_islands: on the
    real 529-strain genus_vs_ureesii study, 10,602 islands qualified but only
    50 were drawn, and the old two-count summary never said where the other
    ~10,552 went.
    """
    islands, matrix, positions = write_inputs(tmp_path)
    out = tmp_path / "island_synteny.html"
    proc = subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out),
         "--top_islands", "0"],
        capture_output=True, text=True, check=True,
    )
    # Only S1's island qualifies (S2 is single-strain, excluded); --top_islands 0
    # then truncates that one qualifying island away from being drawn.
    assert "0 islands drawn" in proc.stderr
    assert "1 excluded" in proc.stderr
    assert "1 truncated" in proc.stderr


def test_cli_reads_zst_inputs(tmp_path):
    # Nextflow publishes family_positions.tsv as .zst since PR #113, and
    # stages it as a SYMLINK -- zstd returns an EMPTY stream for a symlink
    # without -f, which open_maybe_compressed handles.
    import shutil
    if shutil.which("zstd") is None:
        import pytest
        pytest.skip("zstd not available")
    islands, matrix, positions = write_inputs(tmp_path)
    subprocess.run(["zstd", "-q", "-f", str(positions),
                    "-o", str(positions) + ".zst"], check=True)
    link = tmp_path / "staged_positions.tsv.zst"
    link.symlink_to(Path(str(positions) + ".zst"))
    out = tmp_path / "out.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(link),
         "--project", "demo", "--output", str(out)], check=True)
    assert payload_of(out)["islands"][0]["families"] == ["famA", "famB"]


def test_cli_with_no_islands_still_writes_a_page(tmp_path):
    islands = tmp_path / "empty.tsv"
    islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
        "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
        "n_contigs_in_locus\n")
    matrix = tmp_path / "m.tsv"
    matrix.write_text("family\tS1\n")
    positions = tmp_path / "p.tsv"
    positions.write_text("Short\tfamily\tcontig\trank\n")
    out = tmp_path / "empty.html"
    proc = subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert payload_of(out)["islands"] == []


# Amendment tests for --domtblout support

def test_cli_with_domtblout_annotates_family_classes(tmp_path):
    """Pass --domtblout and assert non-"unannotated" family_classes."""
    islands, matrix, positions = write_inputs(tmp_path)

    # Build a domtblout fixture with famA annotated with a real domain
    domtblout = tmp_path / "test.domtblout"
    domtblout.write_text(
        "# hmmscan domtblout\n"
        "PF00001 - 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 1.0e-10 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )

    out = tmp_path / "island_synteny.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--domtblout", str(domtblout),
         "--project", "demo", "--output", str(out)],
        check=True,
    )

    payload = payload_of(out)
    # The island has famA and famB; only famA is annotated in domtblout
    island = payload["islands"][0]
    assert len(island["families"]) == 2
    # famA should have a non-"unannotated" family_class
    assert island["family_classes"][0] != "unannotated"
    # famB should be "unannotated" (no entry in domtblout)
    assert island["family_classes"][1] == "unannotated"


def test_cli_domain_evalue_is_threaded_to_parse_domtblout(tmp_path):
    """I1 regression: --pangenome_pfam_domain_evalue must actually gate the
    glyph strip, not fall back to parse_domtblout's hardcoded 1e-3 default.

    The domtblout fixture's domain-level i-Evalue (column index 12) is
    5.0e-2 -- looser than the default 1e-3, so the default CLI call must
    classify famA as unannotated, while passing a looser --domain_evalue
    must pick the hit up.
    """
    islands, matrix, positions = write_inputs(tmp_path)
    domtblout = tmp_path / "loose.domtblout"
    domtblout.write_text(
        "# hmmscan domtblout\n"
        "PF00001 - 100 famA - 50 1.0e-10 100.0 20.0 1 1 1.0e-10 5.0e-2 100.0 20.0 10 20 10 20 10 20 0.99\n"
    )

    default_out = tmp_path / "default.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--domtblout", str(domtblout),
         "--project", "demo", "--output", str(default_out)],
        check=True,
    )
    default_payload = payload_of(default_out)
    assert default_payload["islands"][0]["family_classes"][0] == "unannotated"

    loose_out = tmp_path / "loose.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--domtblout", str(domtblout),
         "--domain_evalue", "1.0",
         "--project", "demo", "--output", str(loose_out)],
        check=True,
    )
    loose_payload = payload_of(loose_out)
    assert loose_payload["islands"][0]["family_classes"][0] != "unannotated"


def test_cli_orders_multicopy_family_by_locus_contig(tmp_path):
    """C1 regression at the CLI level: family_positions.tsv carries one row
    per (strain, family, COPY); a stray paralog on another contig must not
    reorder the columns away from the island's own locus contig.
    """
    islands = tmp_path / "islands_with_domains.tsv"
    islands.write_text(
        "n_strains\texample_strain\tisland_size\tmember_families\t"
        "n_supporting_pairs\tclassifications\tpfam_domains\tlocus_id\t"
        "locus_contig\tlocus_start\tlocus_end\tn_members_with_coordinates\t"
        "n_contigs_in_locus\n"
        "2\tS1\t2\tfamB,famA\t1\tunexplained_physical\t-\t"
        "S1:c1:100-400\tc1\t100\t400\t2\t1\n"
    )
    matrix = tmp_path / "presence_matrix.tsv"
    matrix.write_text(
        "family\tS1\tS2\nfamA\tpresent\tpresent\nfamB\tpresent\tabsent\n"
    )
    positions = tmp_path / "family_positions.tsv"
    positions.write_text(
        "Short\tfamily\tcontig\trank\n"
        "S1\tfamA\tc1\t1\n"
        "S1\tfamA\tc9\t900\n"
        "S1\tfamB\tc1\t2\n"
    )
    out = tmp_path / "out.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out)],
        check=True,
    )
    payload = payload_of(out)
    assert payload["islands"][0]["families"] == ["famA", "famB"]


def test_cli_without_domtblout_gracefully_degrades(tmp_path):
    """Omit --domtblout and assert all families are "unannotated"."""
    islands, matrix, positions = write_inputs(tmp_path)

    out = tmp_path / "island_synteny.html"
    subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(positions),
         "--project", "demo", "--output", str(out)],
        check=True,
    )

    payload = payload_of(out)
    island = payload["islands"][0]
    # With no domtblout, all families should be "unannotated"
    assert all(fc == "unannotated" for fc in island["family_classes"])


def test_cli_nonexistent_family_positions_fails(tmp_path):
    """Passing a nonexistent --family_positions makes CLI fail (non-zero exit).

    This guards against the silent-failure bug in lib/compressed_io.py where
    a missing .zst file would return an empty stream, causing the CLI to
    silently fall back to alphabetical column order instead of locus order.
    A crash is far better than a plausible-looking but wrong answer.
    """
    islands, matrix, _ = write_inputs(tmp_path)
    nonexistent = tmp_path / "does_not_exist.tsv"
    out = tmp_path / "output.html"

    proc = subprocess.run(
        [sys.executable, str(BIN / "pangenome_island_synteny.py"),
         "--islands_with_domains", str(islands),
         "--presence_matrix", str(matrix),
         "--family_positions", str(nonexistent),
         "--project", "demo", "--output", str(out)],
        capture_output=True, text=True)

    # Must exit with non-zero, not silently succeed
    assert proc.returncode != 0, f"Expected failure but got exit 0. stderr: {proc.stderr}"
