from pathlib import Path
from Bio import SeqIO


def read_fasta(path: str | Path) -> dict[str, object]:
    """Return {seq_id: SeqRecord} from a FASTA file."""
    return SeqIO.to_dict(SeqIO.parse(str(path), 'fasta'))


def write_fasta(records, path: str | Path) -> None:
    SeqIO.write(records, str(path), 'fasta')


def extract_ids(path: str | Path) -> set[str]:
    """Return the set of sequence IDs in a FASTA file without loading sequences."""
    return {rec.id for rec in SeqIO.parse(str(path), 'fasta')}
