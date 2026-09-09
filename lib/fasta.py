import re
import sys
from pathlib import Path
from Bio import SeqIO

# mmseqs2 (createdb/easy-cluster) recognizes UniProt-style "sp|ACC|NAME" /
# "tr|ACC|NAME" FASTA deflines and reports the bare accession (the middle
# field) as the sequence ID in its own *_cluster.tsv, rather than the whole
# first-whitespace-delimited header token every other tool in this pipeline
# uses (Biopython's SeqRecord.id, and the "grep '^>' | sed 's/[[:space:]].*//'"
# protein-map extraction in workflows/profile_search.nf's SEED_PROTEIN_MAP).
# BLAST+ independently does the same UniProt-header recognition, so the
# pairwise (--cluster_tool pairwise) path never notices this divergence --
# it only surfaces for --cluster_tool mmseqs's family-profile pathway, whose
# own bin/ scripts (extract_family_seqs.py, profile_to_matrix.py) need to
# apply this same normalization before joining against a mmseqs cluster.tsv
# id. Confirmed empirically against a real mmseqs2 build (see issue #85) --
# not documentation-derived, since mmseqs' own docs don't spell this out.
_MMSEQS_UNIPROT_ID_RE = re.compile(r'^(?:sp|tr)\|([^|]+)\|')


def mmseqs_id(header_id: str) -> str:
    """Normalize a FASTA header's first token to what mmseqs2 would report as
    that sequence's id in its own *_cluster.tsv output -- a no-op for any
    header that isn't UniProt "sp|ACC|NAME"/"tr|ACC|NAME"-shaped."""
    m = _MMSEQS_UNIPROT_ID_RE.match(header_id)
    return m.group(1) if m else header_id


def read_fasta(path: str | Path) -> dict[str, object]:
    """Return {seq_id: SeqRecord} from a FASTA file.

    Tolerant of duplicate IDs (unlike Bio.SeqIO.to_dict, which raises) -- a real
    occurrence when concatenating proteomes from independently-sourced genome
    collections (e.g. two different genomes reusing the same short internal locus
    tag as their protein ID). The first record for a given ID is kept and the
    collision is reported to stderr; silently guessing which record is "right"
    would be worse than a deterministic, visible first-wins policy.
    """
    records: dict[str, object] = {}
    n_dupes = 0
    for rec in SeqIO.parse(str(path), 'fasta'):
        if rec.id in records:
            n_dupes += 1
            continue
        records[rec.id] = rec
    if n_dupes:
        print(f"WARNING: {path}: {n_dupes} duplicate sequence ID(s) -- kept the "
              f"first occurrence of each, dropped the rest", file=sys.stderr)
    return records


def write_fasta(records, path: str | Path) -> None:
    SeqIO.write(records, str(path), 'fasta')


def extract_ids(path: str | Path) -> set[str]:
    """Return the set of sequence IDs in a FASTA file without loading sequences."""
    return {rec.id for rec in SeqIO.parse(str(path), 'fasta')}
