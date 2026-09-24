import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
sys.path.insert(0, str(REPO / "lib"))
sys.path.insert(0, str(REPO / "tests"))
import uniprot_link as ul  # noqa: E402
from test_uniprot_index import build_fixture_index  # noqa: E402

LINK = REPO / "bin" / "uniprot_link.py"


def test_accession_from_id():
    assert ul.accession_from_id("tr|Q7S6W2|Q7S6W2_NEUCR") == "Q7S6W2"
    assert ul.accession_from_id("sp|P12345-3|X_YEAST") == "P12345"
    assert ul.accession_from_id("Q7S6W2-2") == "Q7S6W2"
    assert ul.accession_from_id("NCU05603T0") is None
    assert ul.accession_from_id("XP_960565.1") is None


def test_choose_order():
    own = {"accession": "B2", "proteome_id": "UP1", "taxid": 1, "reviewed": 0}
    other_rev = {"accession": "A1", "proteome_id": "UP2", "taxid": 2, "reviewed": 1}
    other_unrev = {"accession": "A0", "proteome_id": "UP3", "taxid": 3, "reviewed": 0}
    assert ul.choose([other_rev, own], {1}, None) == (own, "seq_own")
    assert ul.choose([other_unrev, other_rev], {1}, None) == (other_rev, "seq_other")
    assert ul.choose([other_unrev, other_rev], set(), None) == (other_rev, "seq_other")
    assert ul.choose([own, other_rev], {1}, "UP2") == (other_rev, "seq_own")
    assert ul.choose([], {1}, None) == (None, "")


def _run(idx, tmp_path, fasta, *args):
    fa = tmp_path / "p.fa"
    fa.write_text(fasta)
    out = tmp_path / "out.tsv"
    p = subprocess.run([sys.executable, str(LINK), "--index", str(idx), "--protein-fasta", str(fa),
                        *args, "--output", str(out)], capture_output=True, text=True)
    rows = {}
    if p.returncode == 0:
        with open(out) as fh:
            rows = {r["protein_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return p, rows


def test_match_rules_end_to_end(tmp_path):
    idx = build_fixture_index(tmp_path)
    p, rows = _run(idx, tmp_path,
                   ">tr|Q7S6W2|Q7S6W2_NEUCR\nMKVLLAQ\n"   # id
                   ">XP_960565.1\nXXXX\n"                 # refseq
                   ">NCU_gene1\nmkvllaq*\n"               # seq, normalized; own species wins
                   ">NCU_gene2\nMSEQTWO\n"                # seq, own species (only candidate)
                   ">NCU_gene3\nNOMATCH\n",               # none
                   "--short", "Ncra", "--species", "Neurospora crassa")
    assert p.returncode == 0, p.stderr
    assert rows["tr|Q7S6W2|Q7S6W2_NEUCR"]["uniprot_match"] == "id"
    assert rows["tr|Q7S6W2|Q7S6W2_NEUCR"]["uniprot_alphafold_id"] == "Q7S6W2"
    assert rows["XP_960565.1"]["uniprot_match"] == "refseq"
    assert rows["XP_960565.1"]["uniprot_accession"] == "Q7S6W2"
    assert rows["NCU_gene1"]["uniprot_match"] == "seq_own"
    assert rows["NCU_gene1"]["uniprot_accession"] == "Q7S6W2"
    assert rows["NCU_gene1"]["uniprot_n_matches"] == "2"
    assert rows["NCU_gene2"]["uniprot_accession"] == "P00001"
    assert "NCU_gene3" not in rows
    assert list(next(iter(rows.values())).keys()) == list(ul.LINK_COLUMNS)
    assert "id=1 refseq=1 seq_own=2 seq_other=0 none=1" in p.stderr


def test_species_without_own_proteome_falls_to_seq_other(tmp_path):
    idx = build_fixture_index(tmp_path)
    p, rows = _run(idx, tmp_path, ">g1\nMKVLLAQ\n", "--short", "Xxx", "--species", "Other fungus")
    row = rows["g1"]
    assert row["uniprot_match"] == "seq_other"
    assert row["uniprot_accession"] == "P00002"          # reviewed wins the tie
    assert "Saccharomyces cerevisiae" in row["uniprot_match_species"]
    assert "seq_other=1" in p.stderr and "none=0" in p.stderr


def test_explicit_taxid_and_restrict(tmp_path):
    idx = build_fixture_index(tmp_path)
    _, rows = _run(idx, tmp_path, ">g1\nMKVLLAQ\n", "--short", "S", "--species", "Other fungus",
                   "--taxid", "367110")
    assert (rows["g1"]["uniprot_accession"], rows["g1"]["uniprot_match"]) == ("Q7S6W2", "seq_own")
    _, rows = _run(idx, tmp_path, ">g1\nMKVLLAQ\n", "--short", "S", "--species", "Neurospora crassa",
                   "--restrict-proteome", "UP000000002")
    assert (rows["g1"]["uniprot_accession"], rows["g1"]["uniprot_match"]) == ("P00002", "seq_own")


def test_unusable_index_is_a_hard_error(tmp_path):
    p, _ = _run(tmp_path / "nope", tmp_path, ">g1\nMK\n", "--short", "X", "--species", "X y")
    assert p.returncode != 0 and "incomplete" in p.stderr


def test_species_level_taxid_still_counts_same_species_strain_as_own(tmp_path):
    # Config NCBI_TaxID is often species-level (e.g. 5207 C. neoformans) while UniProt
    # reference proteomes are strain-level (235443 H99): a binomial match must still
    # count as own species, not fall to seq_other (which hides gene-database links).
    idx = build_fixture_index(tmp_path)
    _, rows = _run(idx, tmp_path, ">g1\nMKVLLAQ\n", "--short", "Ncra",
                   "--species", "Neurospora crassa", "--taxid", "5141")
    assert (rows["g1"]["uniprot_accession"], rows["g1"]["uniprot_match"]) == ("Q7S6W2", "seq_own")
