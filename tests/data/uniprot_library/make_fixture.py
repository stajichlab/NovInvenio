#!/usr/bin/env python3
"""Regenerate the tiny UniProt library fixture used by tests/test_uniprot_*.py and
-profile test (docs/superpowers/plans/2026-09-23-uniprot-library-index.md, Task 2).
Run from anywhere: writes next to this file."""
import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rec(acc, name, reviewed, taxid, seq, extra_dr=()):
    status = "Reviewed" if reviewed else "Unreviewed"
    dr = "".join(f"DR   {d}\n" for d in extra_dr)
    return (f"ID   {name:<24}{status};{len(seq):>11} AA.\n"
            f"AC   {acc};\n"
            f"DE   SubName: Full=Fixture protein {acc};\n"
            f"OX   NCBI_TaxID={taxid};\n"
            f"RN   [1]\n"
            f"RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].\n"
            f"RX   PubMed=12712197; DOI=10.1038/nature01554;\n"
            f"RT   \"Fixture genome paper.\";\n"
            f"{dr}"
            f"PE   4: Predicted;\n"
            f"SQ   SEQUENCE   {len(seq)} AA;  1 MW;  0000000000000000 CRC64;\n"
            f"     {seq}\n"
            f"//\n")


UP1 = (rec("Q7S6W2", "Q7S6W2_NEUCR", False, 367110, "MKVLLAQ",
           ["RefSeq; XP_960565.1; XM_955472.3.", "VEuPathDB; FungiDB:NCU05603; -.",
            "AlphaFoldDB; Q7S6W2; -.", "PANTHER; PTHR10000; FIXTURE; 1."])
       + rec("P00001", "P00001_NEUCR", True, 367110, "MSEQTWO", ["AlphaFoldDB; P00001; -."]))
UP2 = rec("P00002", "P00002_YEAST", True, 559292, "MKVLLAQ", ["AlphaFoldDB; P00002; -."])

CSV = (
    "proteome_id,tax_id,species_name,oscode,n_canonical,n_isoform,file_prefix,lineage\n"
    '"UP000000001","367110","Neurospora crassa (strain ATCC 24698 / 74-OR23-1A)","NEUCR","2","0","UP000000001_367110","Fungi"\n'
    '"UP000000002","559292","Saccharomyces cerevisiae (strain ATCC 204508 / S288c)","YEAST","1","0","UP000000002_559292","Fungi"\n'
    '"UP000000003","1","Missing fungus","MISS","1","0","UP000000003_1","Fungi"\n'
)

if __name__ == "__main__":
    (HERE / "data").mkdir(exist_ok=True)
    for prefix, text in (("UP000000001_367110", UP1), ("UP000000002_559292", UP2)):
        with gzip.open(HERE / "data" / f"{prefix}.dat.gz", "wt") as fh:
            fh.write(text)
    (HERE / "proteomes.csv").write_text(CSV)
    (HERE / "README").write_text("Fixture UniProt library\nRelease 2099_01, 01-Jan-2099\n")
