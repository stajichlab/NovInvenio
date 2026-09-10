import sys
from pathlib import Path
from Bio import SeqIO

# mmseqs2 (createdb/easy-cluster) recognizes several FASTA defline conventions
# and reports a collapsed field as the sequence ID in its own *_cluster.tsv,
# rather than the whole first-whitespace-delimited header token every other
# tool in this pipeline uses (Biopython's SeqRecord.id, and the
# "grep '^>' | sed 's/[[:space:]].*//'" protein-map extraction in
# workflows/profile_search.nf's SEED_PROTEIN_MAP). BLAST+ independently does
# the same header recognition, so the pairwise (--cluster_tool pairwise) path
# never notices this divergence in its TBLASTN step -- it only surfaces
# wherever a *_cluster.tsv itself is read.
#
# Confirmed empirically against a real mmseqs2 build (issues #85, and this
# session's mmseqs-cluster-ID-restoration design review) -- not
# documentation-derived, since mmseqs's own docs don't spell this out:
#   sp|ACC|NAME, tr|ACC|NAME, gb|ACC|, ref|ACC|, pdb|ACC|CHAIN, bbs|ACC,
#     lcl|ACC, pat|COUNTRY|ACC, cl|ACC          -> second field
#   gnl|db|ACC                                   -> third field
#   pir||ACC, prf||ACC (empty second field)      -> third field
#   gi|N|db|ACC|...  (NCBI's compound defline)    -> fourth field
# Any OTHER pipe-containing header (an unrecognized first field, e.g. this
# repo's own "Pchr|PCH_..." study-specific locus tags) is left completely
# unchanged -- mmseqs does not special-case it, confirmed empirically.
_MMSEQS_SINGLE_PIPE_PREFIXES = {
    'sp', 'tr', 'gb', 'ref', 'pdb', 'bbs', 'lcl', 'pat', 'cl',
}
_MMSEQS_DOUBLE_PIPE_PREFIXES = {'pir', 'prf'}


def mmseqs_id(header_id: str) -> str:
    """Normalize a FASTA header's first token to what mmseqs2 would report as
    that sequence's id in its own *_cluster.tsv output -- a no-op for any
    header whose first field isn't one of mmseqs's recognized prefixes."""
    fields = header_id.split('|')
    if len(fields) < 2:
        return header_id
    prefix = fields[0].lower()
    if prefix == 'gi' and len(fields) >= 4:
        return fields[3]
    if prefix == 'gnl' and len(fields) >= 3:
        return fields[2]
    if prefix in _MMSEQS_DOUBLE_PIPE_PREFIXES and len(fields) >= 3 and fields[1] == '':
        return fields[2]
    if prefix in _MMSEQS_SINGLE_PIPE_PREFIXES:
        return fields[1]
    return header_id


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
