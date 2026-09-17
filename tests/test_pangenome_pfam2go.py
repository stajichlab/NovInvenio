import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import pangenome_pfam2go


PFAM2GO_FIXTURE = """\
!date: today
Pfam:PF13577 SnoaL_2 > GO:oxidoreductase activity ; GO:0016491
Pfam:PF13577 SnoaL_2 > GO:metabolic process ; GO:0008152
Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930
"""


def test_parse_pfam2go_maps_bare_accession_to_go_ids(tmp_path):
    p = tmp_path / "pfam2go"
    p.write_text(PFAM2GO_FIXTURE)
    mapping = pangenome_pfam2go.parse_pfam2go(str(p))
    assert mapping["PF13577"] == ["GO:0016491", "GO:0008152"]
    assert mapping["PF00001"] == ["GO:0004930"]


def test_main_strips_version_suffix_before_lookup(tmp_path):
    pfam2go_path = tmp_path / "pfam2go"
    pfam2go_path.write_text(PFAM2GO_FIXTURE)
    enrichment_path = tmp_path / "island_pfam_enrichment.tsv"
    enrichment_path.write_text(
        "domain\tpfam_accession\tpfam_url\tfisher_p\tfdr_q\n"
        "SnoaL_2\tPF13577.9\thttp://example/PF13577/\t1e-5\t2e-5\n"
        "unknown_domain\t-\t-\t1e-2\t2e-2\n"
    )
    output_path = tmp_path / "output.tsv"
    argv = [
        "pangenome_pfam2go.py",
        "--island_pfam_enrichment", str(enrichment_path),
        "--pfam2go", str(pfam2go_path),
        "--output", str(output_path),
    ]
    import sys as _sys
    old_argv = _sys.argv
    _sys.argv = argv
    try:
        pangenome_pfam2go.main()
    finally:
        _sys.argv = old_argv

    with open(output_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["domain"] == "SnoaL_2"
    assert rows[0]["n_go_terms"] == "2"
    assert rows[0]["go_terms"] == "GO:0016491;GO:0008152"
    assert rows[1]["domain"] == "unknown_domain"
    assert rows[1]["n_go_terms"] == "0"
    assert rows[1]["go_terms"] == "-"
