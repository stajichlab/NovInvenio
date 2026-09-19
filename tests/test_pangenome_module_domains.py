import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_module_domains import load_family_modules, parse_domtblout  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin" / "pangenome_module_domains.py"


def test_load_family_modules_parses_module_id_and_size(tmp_path):
    f = tmp_path / "family_modules.tsv"
    f.write_text("family\tmodule_id\tmodule_size\nfamA\t0\t2\nfamB\t0\t2\nfamC\t1\t1\n")
    assert load_family_modules(str(f)) == {
        "famA": ("0", 2), "famB": ("0", 2), "famC": ("1", 1),
    }


def test_parse_domtblout_maps_query_to_domain_names_and_accessions(tmp_path):
    domtblout = tmp_path / "pfam.domtblout"
    domtblout.write_text(
        "# comment line\n"
        "NACHT PF05729.19 137 famA - 214 1e-20 60 0 1 1 1e-20 1e-20 60 0 1 10 1 10 1 10 0.9 desc\n"
        "Ank_2 PF12796.14 30 famA - 214 1e-10 40 0 1 1 1e-10 1e-10 40 0 1 10 1 10 1 10 0.9 desc\n"
    )
    hits = parse_domtblout(str(domtblout))
    assert hits["famA"] == {("NACHT", "PF05729.19"), ("Ank_2", "PF12796.14")}


def test_main_writes_per_module_domain_summary(tmp_path):
    family_modules = tmp_path / "family_modules.tsv"
    family_modules.write_text(
        "family\tmodule_id\tmodule_size\n"
        "famA\t0\t2\n"
        "famB\t0\t2\n"
        "famC\t1\t1\n"
    )
    domtblout = tmp_path / "pfam.domtblout"
    domtblout.write_text(
        "NACHT PF05729.19 137 famA - 214 1e-20 60 0 1 1 1e-20 1e-20 60 0 1 10 1 10 1 10 0.9 desc\n"
        "MFS_1 PF07690.19 137 famB - 214 1e-20 60 0 1 1 1e-20 1e-20 60 0 1 10 1 10 1 10 0.9 desc\n"
    )
    output = tmp_path / "module_domains.tsv"

    result = subprocess.run(
        [sys.executable, str(BIN),
         "--family_modules", str(family_modules),
         "--domtblout", str(domtblout),
         "--min_module_size", "2",
         "--output", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    text = output.read_text()
    assert text.startswith("module_id\tmodule_size\tn_families_with_domain\tdomain_names\n")
    assert "0\t2\t2\t" in text
    assert "NACHT" in text and "MFS_1" in text
    # module 1 (size 1) excluded by --min_module_size 2
    assert "\n1\t1\t" not in text


def test_main_handles_empty_family_modules_gracefully(tmp_path):
    # Same graceful-degradation contract as pangenome_detect_trans_modules.py:
    # a study with zero trans edges produces a header-only family_modules.tsv
    # (no rows), which must not crash this downstream step either.
    family_modules = tmp_path / "family_modules.tsv"
    family_modules.write_text("family\tmodule_id\tmodule_size\n")
    domtblout = tmp_path / "pfam.domtblout"
    domtblout.write_text("")
    output = tmp_path / "module_domains.tsv"

    result = subprocess.run(
        [sys.executable, str(BIN),
         "--family_modules", str(family_modules),
         "--domtblout", str(domtblout),
         "--output", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output.read_text() == "module_id\tmodule_size\tn_families_with_domain\tdomain_names\n"
