#!/usr/bin/env python3
"""Match one proteome to UniProt records via a library index (lib/uniprot_index.py).

Match order per protein, first rule that matches wins:
  id        the protein ID is a UniProt accession (sp|ACC|..., tr|ACC|..., bare ACC,
            isoform suffix -N stripped) present in the index
  refseq    the protein ID (version kept) is a RefSeq ID on a record's DR RefSeq line
  seq_own   exact sequence (lib/uniprot_dat.seq_md5) match in the species' own taxid
            (--taxid, else species binomial lookup), or in --restrict-proteome
  seq_other exact sequence match in any other proteome
Ties: own first, then reviewed (Swiss-Prot), then lowest accession.

Writes <Short>.uniprot_link.tsv keyed by protein_id (unmatched proteins omitted)
and a coverage summary on stderr. See
docs/superpowers/specs/2026-09-23-uniprot-library-index-design.md."""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed  # noqa: E402
from uniprot_dat import seq_md5  # noqa: E402
from uniprot_index import IndexFormatError, UniProtIndex  # noqa: E402

ACC_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
LINK_COLUMNS = (
    "protein_id", "uniprot_accession", "uniprot_gene_name", "uniprot_description",
    "uniprot_go_ids", "uniprot_pfam_ids", "uniprot_pfam_names", "uniprot_interpro_ids",
    "uniprot_ec_numbers", "uniprot_alphafold_id", "uniprot_xrefs", "uniprot_match",
    "uniprot_match_species", "uniprot_reviewed", "uniprot_pubs", "uniprot_n_matches",
)
REC_TO_COL = {"gene_name": "uniprot_gene_name", "description": "uniprot_description",
              "go_ids": "uniprot_go_ids", "pfam_ids": "uniprot_pfam_ids",
              "pfam_names": "uniprot_pfam_names", "interpro_ids": "uniprot_interpro_ids",
              "ec_numbers": "uniprot_ec_numbers", "alphafold_id": "uniprot_alphafold_id",
              "xrefs": "uniprot_xrefs", "reviewed": "uniprot_reviewed", "pubs": "uniprot_pubs"}
MATCH_TYPES = ("id", "refseq", "seq_own", "seq_other", "none")


def accession_from_id(pid):
    """UniProt accession in a protein ID, isoform suffix stripped, else None."""
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", pid)
    token = m.group(1) if m else pid.split()[0]
    token = token.split("-")[0]
    return token if ACC_RE.match(token) else None


def choose(cands, own_taxid, restrict):
    """(best candidate, 'seq_own'|'seq_other') from by_md5() hits, or (None, '')."""
    if not cands:
        return None, ""

    def own(c):
        if restrict:
            return c["proteome_id"] == restrict
        return own_taxid is not None and c["taxid"] == own_taxid

    best = sorted(cands, key=lambda c: (not own(c), -int(c["reviewed"]), c["accession"]))[0]
    return best, ("seq_own" if own(best) else "seq_other")


def read_fasta(path):
    pid, chunks = None, []
    with open_maybe_compressed(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(chunks)
                pid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if pid is not None:
        yield pid, "".join(chunks)


def match_proteome(idx, fasta, own_taxid, restrict):
    """{protein_id: (hit dict, match type, n candidates)}, Counter of match types."""
    chosen, counts = {}, Counter()
    for pid, seq in read_fasta(fasta):
        acc = accession_from_id(pid)
        hit = idx.by_accession(acc) if acc else None
        if hit:
            chosen[pid] = (hit, "id", 1)
        else:
            refs = idx.by_refseq(pid)
            hit = idx.by_accession(refs[0]) if refs else None
            if hit:
                chosen[pid] = (hit, "refseq", len(refs))
            else:
                cands = idx.by_md5(seq_md5(seq))
                best, kind = choose(cands, own_taxid, restrict)
                if best:
                    chosen[pid] = (best, kind, len(cands))
        counts[chosen[pid][1] if pid in chosen else "none"] += 1
    return chosen, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--protein-fasta", required=True, type=Path)
    ap.add_argument("--short", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--taxid", type=int, default=None)
    ap.add_argument("--restrict-proteome", default=None, dest="restrict",
                    help="UniProt proteome id (UP...) that counts as own species")
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    try:
        idx = UniProtIndex(a.index)
    except IndexFormatError as e:
        sys.exit(f"ERROR: {e}")
    own_taxid = a.taxid if a.taxid is not None else idx.taxid_for_species(a.species)
    chosen, counts = match_proteome(idx, a.protein_fasta, own_taxid, a.restrict)

    wanted = defaultdict(set)
    for hit, _, _ in chosen.values():
        wanted[hit["proteome_id"]].add(hit["accession"])
    recs = {}
    for prot, accs in wanted.items():
        for acc, r in idx.records(prot, accs).items():
            recs[(prot, acc)] = r

    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LINK_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for pid, (hit, kind, n) in chosen.items():
            r = recs.get((hit["proteome_id"], hit["accession"]))
            if r is None:
                continue
            row = dict.fromkeys(LINK_COLUMNS, "")
            row.update(protein_id=pid, uniprot_accession=hit["accession"], uniprot_match=kind,
                       uniprot_match_species=idx.species_name(hit["proteome_id"]),
                       uniprot_n_matches=str(n))
            for k, col in REC_TO_COL.items():
                row[col] = r[k]
            w.writerow(row)
    n_total = sum(counts.values())
    summary = " ".join(f"{k}={counts[k]}" for k in MATCH_TYPES)
    print(f"{a.short}: {n_total} proteins; {summary} (own taxid {own_taxid})", file=sys.stderr)


if __name__ == "__main__":
    main()
