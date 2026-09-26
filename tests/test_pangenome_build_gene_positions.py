import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from pangenome_build_gene_positions import (
    parse_gff3_protein_positions,
    load_protein_ids,
    resolve_positions_for_strain,
    load_protein_gene_names,
)


def test_parse_gff3_protein_positions_ncbi_style_unchanged(tmp_path):
    gff3 = tmp_path / "ncbi.gff3"
    gff3.write_text(
        "chr1\tsrc\tCDS\t100\t200\t.\t+\t0\tID=cds1;protein_id=P1;\n"
        "chr1\tsrc\tCDS\t300\t400\t.\t+\t0\tID=cds2;protein_id=P1;\n"
    )
    positions = parse_gff3_protein_positions(gff3)
    assert positions == {"P1": ("chr1", 100, 400)}


def test_parse_gff3_protein_positions_funannotate_parent_fallback(tmp_path):
    gff3 = tmp_path / "funannotate.gff3"
    gff3.write_text(
        "contig1\tfunannotate\tCDS\t893\t1495\t.\t+\t0\tID=P1.cds;Parent=P1;\n"
        "contig1\tfunannotate\tCDS\t1600\t1800\t.\t+\t0\tID=P1.cds;Parent=P1;\n"
        "contig1\tfunannotate\tCDS\t2000\t2100\t.\t+\t0\tID=P2.cds;Parent=P2;\n"
    )
    positions = parse_gff3_protein_positions(gff3)
    assert positions == {
        "P1": ("contig1", 893, 1800),
        "P2": ("contig1", 2000, 2100),
    }


def test_parse_gff3_protein_positions_handles_multi_parent_cds(tmp_path):
    # GFF3 spec allows Parent=A,B for a CDS row shared by multiple
    # transcripts -- must not treat "A,B" as one literal bogus ID.
    gff3 = tmp_path / "multi_parent.gff3"
    gff3.write_text(
        "contig1\tfunannotate\tCDS\t100\t200\t.\t+\t0\tID=shared.cds;Parent=P1,P2;\n"
    )
    positions = parse_gff3_protein_positions(gff3)
    assert positions == {
        "P1": ("contig1", 100, 200),
        "P2": ("contig1", 100, 200),
    }


def test_load_protein_ids_parses_fasta_headers(tmp_path):
    fasta = tmp_path / "strain.pep.fa"
    fasta.write_text(">P1 some description\nMSEQ\n>P2\nMSEQ2\n")
    assert load_protein_ids(fasta) == {"P1", "P2"}


def test_resolve_positions_for_strain_filters_to_fasta_matches(tmp_path):
    gff3 = tmp_path / "strain.gff3"
    gff3.write_text(
        "contig1\tfunannotate\tCDS\t100\t200\t.\t+\t0\tParent=P1;\n"
        "contig1\tfunannotate\tCDS\t300\t400\t.\t+\t0\tParent=STALE_ID;\n"
    )
    protein_ids = {"P1"}  # STALE_ID isn't a real protein in this strain's FASTA
    result = resolve_positions_for_strain(gff3, protein_ids, short="s1")
    assert result == {"P1": ("contig1", 100, 200)}


def test_resolve_positions_for_strain_hard_errors_on_mostly_unmatched(tmp_path):
    # Simulates a genuine dialect mismatch: GFF3 IDs don't correspond to the
    # protein FASTA at all (e.g. wrong attribute picked up entirely).
    gff3 = tmp_path / "mismatched.gff3"
    gff3.write_text(
        "contig1\tsrc\tCDS\t100\t200\t.\t+\t0\tParent=WRONG1;\n"
        "contig1\tsrc\tCDS\t300\t400\t.\t+\t0\tParent=WRONG2;\n"
    )
    protein_ids = {"P1", "P2", "P3", "P4"}  # none of these match WRONG1/WRONG2
    with pytest.raises(SystemExit):
        resolve_positions_for_strain(gff3, protein_ids, short="s1")


def test_resolve_positions_for_strain_warns_but_continues_above_hard_threshold(tmp_path, capsys):
    # 90/100 proteins resolved (10% unresolved): above the 50% hard-error
    # floor, but above the 2% warn threshold -- should warn, not abort.
    gff3 = tmp_path / "mostly_ok.gff3"
    lines = [f"contig1\tsrc\tCDS\t{100*i}\t{100*i+50}\t.\t+\t0\tParent=P{i};\n" for i in range(1, 91)]
    gff3.write_text("".join(lines))
    protein_ids = {f"P{i}" for i in range(1, 101)}
    result = resolve_positions_for_strain(gff3, protein_ids, short="s1")
    assert len(result) == 90
    captured = capsys.readouterr()
    assert "WARNING" in captured.err


# ---- UniProt FASTA + NCBI GFF3 (issue #187) ------------------------------------
# The FASTA is keyed by UniProt IDs with the locus tag in GN=; the GFF3 is keyed by
# GenBank protein_id= and carries locus_tag= on each CDS. Neither protein_id= nor
# Parent= can match, so the positions must come through GN= <-> locus_tag=.

UNIPROT_FASTA = (
    ">tr|Q4WAA1|Q4WAA1_ASPFU Some protein OS=Aspergillus fumigatus OX=330879 GN=AFUA_1G00100 PE=4 SV=1\nMK\n"
    ">tr|Q4WAA2|Q4WAA2_ASPFU Other protein OS=Aspergillus fumigatus OX=330879 GN=AFUA_1G00200 PE=4 SV=1\nMK\n"
    ">tr|Q4WAA3|Q4WAA3_ASPFU Dup name A OS=Aspergillus fumigatus GN=AFUA_1G00300 PE=4 SV=1\nMK\n"
    ">tr|Q4WAA4|Q4WAA4_ASPFU Dup name B OS=Aspergillus fumigatus GN=AFUA_1G00300 PE=4 SV=1\nMK\n"
)
NCBI_GFF3 = (
    "chr1\tGenbank\tCDS\t100\t200\t.\t+\t0\tID=cds-XP_1;protein_id=XP_1;locus_tag=AFUA_1G00100\n"
    "chr1\tGenbank\tCDS\t250\t300\t.\t+\t0\tID=cds-XP_1;protein_id=XP_1;locus_tag=AFUA_1G00100\n"
    "chr1\tGenbank\tCDS\t500\t600\t.\t-\t0\tID=cds-XP_2;protein_id=XP_2;locus_tag=AFUA_1G00200\n"
    "chr1\tGenbank\tCDS\t900\t950\t.\t+\t0\tID=cds-XP_3;protein_id=XP_3;locus_tag=AFUA_1G00300\n"
)


def _uniprot_case(tmp_path):
    fa = tmp_path / "uniprot.fa"
    fa.write_text(UNIPROT_FASTA)
    gff3 = tmp_path / "ncbi.gff3"
    gff3.write_text(NCBI_GFF3)
    return fa, gff3


def test_load_protein_gene_names_keeps_only_unique_gn(tmp_path):
    fa, _ = _uniprot_case(tmp_path)
    gn = load_protein_gene_names(fa)
    assert gn == {"AFUA_1G00100": "tr|Q4WAA1|Q4WAA1_ASPFU",
                  "AFUA_1G00200": "tr|Q4WAA2|Q4WAA2_ASPFU"}


def test_uniprot_fasta_with_ncbi_gff3_resolves_through_locus_tag(tmp_path, capsys):
    fa, gff3 = _uniprot_case(tmp_path)
    result = resolve_positions_for_strain(gff3, load_protein_ids(fa), short="Afum",
                                          gene_names=load_protein_gene_names(fa))
    assert result == {"tr|Q4WAA1|Q4WAA1_ASPFU": ("chr1", 100, 300),
                      "tr|Q4WAA2|Q4WAA2_ASPFU": ("chr1", 500, 600)}
    assert "locus_tag" in capsys.readouterr().err


def test_locus_tag_fallback_is_not_used_when_protein_ids_match(tmp_path):
    gff3 = tmp_path / "ncbi.gff3"
    gff3.write_text(NCBI_GFF3)
    result = resolve_positions_for_strain(gff3, {"XP_1", "XP_2", "XP_3"}, short="s",
                                          gene_names={"AFUA_1G00100": "other"})
    assert set(result) == {"XP_1", "XP_2", "XP_3"}


def test_still_hard_errors_when_no_dialect_matches(tmp_path):
    fa, gff3 = _uniprot_case(tmp_path)
    with pytest.raises(SystemExit):
        resolve_positions_for_strain(gff3, load_protein_ids(fa), short="Afum",
                                     gene_names={})
