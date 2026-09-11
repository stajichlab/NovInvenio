#!/usr/bin/env python3
"""Cross-walk an NCBI RefSeq-sourced proteome to its UniProt record via the DR RefSeq line.

Reverse of NovInvenio_Investigations' bin/extract_dat_annotations.py: that script
assumes the pipeline's own protein_id IS the UniProt accession (a UniProt-sourced
proteome). Here the pipeline's protein_id is an NCBI RefSeq accession (XP_...), so we
key off the UniProt record's own `DR RefSeq;` cross-reference line instead of its `AC`
line, and re-emit the annotation keyed by the NCBI id so it can be merged into
presence_matrix.function.tsv by protein_id like any other annotation source
(bin/annotate_presence_matrix.py's --uniprot_xref_files).

Coverage is inherently partial: only UniProt records that carry a RefSeq DR line at
all, and whose RefSeq accession+version matches this run's own protein FASTA exactly
(no version-drift reconciliation attempted), are resolved. A species with no UniProt
proteome at all (e.g. Schizophyllum commune, checked 2026-09-10 -- its only UniProt
proteome entry is non-reference and unindexed) simply has no --uniprot_dat_gz in the
config and is skipped upstream; that is not an error. Coverage on other real species
measured 2026-09-10 (configs/agaricomycetes_v1.csv): 89-100%.

Prints a coverage summary to stderr so a silent partial match isn't mistaken for a
full one.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from uniprot_dat import parse_dat_gz  # noqa: E402


def load_refseq_ids(fasta_path: Path) -> list[str]:
    """First whitespace token of every '>' header, in file order."""
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dat-gz", required=True, type=Path, help="UniProt {Proteome}_{taxid}.dat.gz")
    ap.add_argument("--protein-fasta", required=True, type=Path,
                     help="This run's NCBI RefSeq protein FASTA (protein_id space to key output by)")
    ap.add_argument("--short", required=True, help="Config Short code, for the summary line only")
    ap.add_argument("--output", required=True, type=Path, help="Output TSV, keyed by NCBI protein_id")
    args = ap.parse_args()

    refseq_to_record = {}
    n_records = 0
    for rec in parse_dat_gz(args.dat_gz):
        n_records += 1
        for xref in rec["xrefs"].split("|"):
            if xref.startswith("RefSeq:"):
                refseq_to_record[xref.split(":", 1)[1]] = rec
                break  # parse_dat_gz already dedups RefSeq to one line per record

    our_ids = load_refseq_ids(args.protein_fasta)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["protein_id", "uniprot_accession", "uniprot_gene_name", "uniprot_description",
              "uniprot_go_ids", "uniprot_pfam_ids", "uniprot_pfam_names", "uniprot_interpro_ids",
              "uniprot_ec_numbers", "uniprot_alphafold_id", "uniprot_xrefs"]
    n_matched = 0
    with open(args.output, "w") as out:
        out.write("\t".join(fields) + "\n")
        for pid in our_ids:
            rec = refseq_to_record.get(pid)
            if rec is None:
                continue
            n_matched += 1
            row = [pid, rec["accession"], rec["gene_name"], rec["description"],
                   rec["go_ids"], rec["pfam_ids"], rec["pfam_names"],
                   rec["interpro_ids"], rec["ec_numbers"], rec["alphafold_id"],
                   rec["xrefs"]]
            out.write("\t".join(row) + "\n")

    pct = (100.0 * n_matched / len(our_ids)) if our_ids else 0.0
    print(f"{args.short}: {n_matched}/{len(our_ids)} proteins ({pct:.1f}%) matched a "
          f"UniProt record via RefSeq DR cross-reference ({n_records} UniProt records scanned)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
