"""Schema validation of the Starship element and SM cluster tables (lib/starship_tables.py)."""
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
import starship_tables as st  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "docs" / "examples" / "starship_overlap"


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_coordinate_elements_load_and_alias(tmp_path):
    p = _w(tmp_path, "e.tsv", "# comment\nelement_id\tstrain\tcontig\tstart\tend\n\nE1\tA\tc1\t5\t10\n")
    els = st.load_elements(p, {"A", "B"})
    assert els == [{"element_id": "E1", "label": "", "confidence": "", "Short": "A", "contig": "c1",
                    "start": 5, "end": 10, "genes": None}]
    p2 = _w(tmp_path, "e2.tsv", "element_id\tShort\tcontig\tstart\tend\nE1\tA\tc1\t5\t10\n")
    assert st.load_elements(p2)[0]["Short"] == "A"


def test_csv_and_gz(tmp_path):
    p = _w(tmp_path, "e.csv", "element_id,strain,contig,start,end,label\nE1,A,c1,5,10,\"Ship, h1\"\n")
    assert st.load_elements(p)[0]["label"] == "Ship, h1"
    gz = tmp_path / "e.tsv.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("element_id\tstrain\tcontig\tstart\tend\nE1\tA\tc1\t5\t10\n")
    assert st.load_elements(gz)[0]["end"] == 10


@pytest.mark.parametrize("body, msg", [
    ("E1\tA\tc1\tx\t10\n", "'x' is not an integer"),
    ("E1\tA\tc1\t20\t10\n", "start 20 > end 10"),
    ("E1\tA\tc1\t0\t10\n", "< 1"),
    ("E1\tZ\tc1\t1\t10\n", "strain 'Z' is not a Short ID"),
    ("E1\tA\tc1\t1\t10\nE1\tA\tc1\t1\t10\n", "duplicate element_id 'E1'"),
    ("E1\tA\t\t1\t10\n", "empty contig"),
])
def test_element_errors(tmp_path, body, msg):
    p = _w(tmp_path, "e.tsv", "element_id\tstrain\tcontig\tstart\tend\n" + body)
    with pytest.raises(st.TableError, match=msg):
        st.load_elements(p, {"A"})


def test_missing_columns(tmp_path):
    p = _w(tmp_path, "e.tsv", "element_id\tstrain\tcontig\tstart\nE1\tA\tc1\t1\n")
    with pytest.raises(st.TableError, match=r"missing column\(s\) \['end'\]"):
        st.load_elements(p)
    p = _w(tmp_path, "n.tsv", "element_id\tcontig\tstart\tend\nE1\tc1\t1\t2\n")
    with pytest.raises(st.TableError, match="missing column 'strain'"):
        st.load_elements(p)
    p = _w(tmp_path, "empty.tsv", "")
    with pytest.raises(st.TableError, match="empty table"):
        st.load_elements(p)


def test_gene_level_elements(tmp_path):
    p = _w(tmp_path, "g.tsv", "strain\tgene_id\telement_id\nA\tg1\tE1\nA\tg2\tE1\nA\tg9\tE2\n")
    els = {e["element_id"]: e for e in st.load_elements(p, {"A"})}
    assert els["E1"]["genes"] == ["g1", "g2"] and els["E1"]["start"] is None
    bad = _w(tmp_path, "b.tsv", "strain\tgene_id\telement_id\nA\tg1\tE1\nB\tg2\tE1\n")
    with pytest.raises(st.TableError, match="two strains"):
        st.load_elements(bad)


def test_sm_clusters(tmp_path):
    p = _w(tmp_path, "s.tsv", "cluster_id\tstrain\tcontig\tstart\tend\ttype\tgenes\tcore_genes\n"
                              "K1\tA\tc1\t1\t100\tNRPS\tg1, g2,\tg1\nK2\tA\tc1\t5\t9\tterpene\t\t\n")
    cs = st.load_sm_clusters(p, {"A"})
    assert cs[0]["genes"] == ["g1", "g2"] and cs[0]["core_genes"] == ["g1"]
    assert cs[1]["genes"] == [] and cs[1]["contig_edge"] == ""
    bad = _w(tmp_path, "b.tsv", "cluster_id\tstrain\tcontig\tstart\tend\ttype\nK1\tA\tc1\t9\t5\t\n")
    with pytest.raises(st.TableError) as exc:
        st.load_sm_clusters(bad)
    assert "start 9 > end 5" in str(exc.value) and "empty type" in str(exc.value)
    nocol = _w(tmp_path, "n.tsv", "cluster_id\tstrain\tcontig\tstart\tend\nK1\tA\tc1\t1\t5\n")
    with pytest.raises(st.TableError, match=r"missing column\(s\) \['type'\]"):
        st.load_sm_clusters(nocol)


def test_error_list_is_capped(tmp_path):
    body = "".join(f"E{i}\tA\tc1\tx\t1\n" for i in range(40))
    p = _w(tmp_path, "e.tsv", "element_id\tstrain\tcontig\tstart\tend\n" + body)
    with pytest.raises(st.TableError, match=r"and 20 more"):
        st.load_elements(p)


def test_samplesheet_strains(tmp_path):
    p = _w(tmp_path, "ss.csv", "GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n"
                               "IN,a,,a.fa,a.fa,a.gff3,A,\nOUT,b,,b.fa,b.fa,b.gff3,B,\nX,c,,c,c,c,C,\n")
    assert st.samplesheet_strains(p) == {"A", "B", "C"}
    assert st.samplesheet_strains(p, ["IN", "OUT"]) == {"A", "B"}


def test_shipped_examples_validate():
    els = st.load_elements(EXAMPLES / "af293_starships_table_s7.elements.tsv")
    assert len(els) == 8 and {e["Short"] for e in els} == {"Asfu_Af293"}
    cs = st.load_sm_clusters(EXAMPLES / "sm_clusters.example.tsv")
    assert cs and all(c["type"] for c in cs)
