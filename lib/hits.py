import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


def open_input(path):
    """Open a file for reading, transparently decompressing if .gz."""
    p = Path(path)
    if p.suffix == '.gz':
        return gzip.open(p, 'rt')
    return p.open()


@dataclass
class Hit:
    query_id: str
    target_id: str
    evalue: float
    bitscore: float
    query_proteome: str = ''
    target_proteome: str = ''


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_phmmer_tblout(fh) -> Iterator[Hit]:
    """Parse phmmer --tblout output (per-sequence hits, not per-domain)."""
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split()
        # col 0: target_name, col 2: query_name, col 4: full-seq evalue, col 5: score
        yield Hit(
            query_id=cols[2],
            target_id=cols[0],
            evalue=float(cols[4]),
            bitscore=float(cols[5]),
        )


def parse_tabular(fh) -> Iterator[Hit]:
    """Parse diamond / BLAST outfmt-6 with fields: qseqid sseqid evalue bitscore."""
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        yield Hit(
            query_id=cols[0],
            target_id=cols[1],
            evalue=float(cols[2]),
            bitscore=float(cols[3]),
        )


PARSERS = {
    'phmmer':  parse_phmmer_tblout,
    'diamond': parse_tabular,
    'blast':   parse_tabular,
}


# ── Filters ───────────────────────────────────────────────────────────────────

def filter_hits(
    hits: Iterator[Hit],
    evalue_cutoff: float = 1e-5,
    bitscore_cutoff: float | None = None,
) -> Iterator[Hit]:
    for hit in hits:
        if hit.evalue > evalue_cutoff:
            continue
        if bitscore_cutoff is not None and hit.bitscore < bitscore_cutoff:
            continue
        yield hit
