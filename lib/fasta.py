import sys
from pathlib import Path
from Bio import SeqIO


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
