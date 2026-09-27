"""bin/pangenome_island_dna_check.py (spec section 4b) on synthetic genomes.
The exemplar EX carries flank + gene A + gene B + flank; KEEP has the same
DNA (a gene-model difference), GONE lacks both genes (a deletion), HALF
lacks gene B. Needs blastn (the NovInvenio pixi env has it)."""
import csv
import gzip
import random
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from pangenome_island_dna_check import main, read_fasta_slices  # noqa: E402

rng = random.Random(7)


def dna(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


LEFT, GENE_A, GENE_B, RIGHT = dna(1000), dna(600), dna(700), dna(1000)
SPACER = dna(100)


def write_case(d: Path, gz: bool = False) -> Path:
    genomes = {
        "EX": LEFT + GENE_A + SPACER + GENE_B + RIGHT,
        "KEEP": LEFT + GENE_A + SPACER + GENE_B + RIGHT,
        "GONE": LEFT + SPACER + RIGHT,
        "HALF": LEFT + GENE_A + SPACER + RIGHT,
        "INFLANK": LEFT + GENE_A + SPACER + GENE_B + RIGHT,
    }
    (d / "dna").mkdir()
    lines = ["GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup"]
    for s, seq in genomes.items():
        name = f"{s}.dna.fa" + (".gz" if gz else "")
        text = f">other desc\nACGT\n>ctg{s} x\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
        if gz:
            with gzip.open(d / "dna" / name, "wt") as fh:
                fh.write(text)
        else:
            (d / "dna" / name).write_text(text)
        lines.append(f"IN,Sp,{s},{s}.pep.fa,{name},{s}.gff3,{s},t")
    lines.append("IN,Sp,NOFILE,x.pep.fa,missing.fa,x.gff3,NOFILE,t")
    (d / "config.csv").write_text("\n".join(lines) + "\n")
    a0, a1 = 1001, 1600
    b0, b1 = 1701, 2400
    rows = ["locus_id\trole\tstrain\tcontig\tstart\tend\tgenes",
            f"L\tquery\tEX\tctgEX\t{a0}\t{b1}\t5:{a0}-{a1};6:{b0}-{b1}",
            f"L\ttarget\tKEEP\tctgKEEP\t1001\t{b1}\t",
            "L\ttarget\tGONE\tctgGONE\t1001\t1100\t",
            "L\ttarget\tHALF\tctgHALF\t1001\t1700\t",
            "L\ttarget\tINFLANK\tctgINFLANK\t1\t3400\t",
            "L\ttarget\tSHORT\tctgGONE\t1001\t1040\t",
            "L\ttarget\tNOPAIR\t-\t0\t0\t",
            "L\ttarget\tNOFILE\tctgX\t1\t500\t"]
    (d / "targets.tsv").write_text("\n".join(rows) + "\n")
    return d / "calls.tsv"


def calls(path: Path) -> dict:
    with open(path) as fh:
        return {(r["strain"], r["col"]): r["status"] for r in csv.DictReader(fh, delimiter="\t")}


def run(d: Path, blastn: str = "blastn") -> dict:
    out = d / "calls.tsv"
    assert main(["--targets", str(d / "targets.tsv"), "--config", str(d / "config.csv"),
                 "--data_dir", str(d), "--blastn", blastn, "--cpus", "2",
                 "--output", str(out)]) == 0
    return calls(out)


needs_blastn = pytest.mark.skipif(not shutil.which("blastn"), reason="blastn not on PATH")


@needs_blastn
def test_same_dna_is_present_and_a_deletion_is_absent(tmp_path):
    write_case(tmp_path)
    got = run(tmp_path)
    assert got[("KEEP", "5")] == "present" and got[("KEEP", "6")] == "present"
    assert got[("GONE", "5")] == "absent" and got[("GONE", "6")] == "absent"
    assert got[("HALF", "5")] == "present" and got[("HALF", "6")] == "absent"
    # Target = flank gene start .. flank gene end (spec 4b): locus DNA inside
    # a flank gene model is found.
    assert got[("INFLANK", "5")] == "present" and got[("INFLANK", "6")] == "present"


@needs_blastn
def test_gzipped_genomes_are_read(tmp_path):
    write_case(tmp_path, gz=True)
    assert run(tmp_path)[("KEEP", "6")] == "present"


def test_short_target_is_absent_and_missing_ones_are_unchecked_without_blastn(tmp_path):
    write_case(tmp_path)
    # Leave only targets that need no alignment; a blastn call would fail.
    rows = (tmp_path / "targets.tsv").read_text().splitlines()
    keep = [r for r in rows if not any(f"\t{s}\t" in r for s in ("KEEP", "GONE", "HALF", "INFLANK"))]
    (tmp_path / "targets.tsv").write_text("\n".join(keep) + "\n")
    got = run(tmp_path, blastn=str(tmp_path / "no-such-blastn"))
    assert got[("SHORT", "5")] == "absent" and got[("SHORT", "6")] == "absent"
    assert got[("NOPAIR", "5")] == "unchecked"
    assert got[("NOFILE", "6")] == "unchecked"


def test_read_fasta_slices_uses_the_first_header_word(tmp_path):
    p = tmp_path / "g.fa"
    p.write_text(">c1 desc\nACGT\nTTGG\n>c2\nAAAA\n")
    assert read_fasta_slices(str(p), {"c1": {(3, 6)}, "c9": {(1, 2)}}) == {("c1", 3, 6): "GTTT"}
