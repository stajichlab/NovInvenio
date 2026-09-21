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
    """query_id/target_id/evalue/bitscore are the only fields every hit file (old or
    new, any of the three tools) is guaranteed to carry -- every caller that predates
    issue #129 uses only these four and must keep working unchanged.

    length/pident/qcov/scov/qlen/slen are populated ONLY for a diamond/blastp hit file
    written with the wider --outfmt (issue #129); they default to None -- never 0, which
    would look like a real "zero coverage" measurement -- for an old 4-column cached
    file, a phmmer hit (--tblout carries no alignment geometry at all -- see issue #150
    for --domtblout follow-up), or a caller that never asked for them. A consumer of
    these fields MUST treat None as "not measured" and fail loudly rather than default
    it to 0/100/pass, per the coverage-floor design in the paralog-rescue spec
    (docs/superpowers/specs/2026-09-20-paralog-novelty-disqualifier-analysis.md, Option 4).
    """
    query_id: str
    target_id: str
    evalue: float
    bitscore: float
    query_proteome: str = ''
    target_proteome: str = ''
    length: int | None = None
    pident: float | None = None
    qcov: float | None = None      # query coverage, percent (diamond qcovhsp / blastp qcovhsp)
    scov: float | None = None      # subject coverage, percent (diamond scovhsp only -- blastp
                                    # has no equivalent keyword)
    qlen: int | None = None
    slen: int | None = None


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
    """Parse diamond / BLAST outfmt-6.

    Reads the base 4-column layout (qseqid sseqid evalue bitscore) that every cached
    hit file is guaranteed to have, plus -- when present -- the wider layout added for
    issue #129: length pident qcovhsp [scovhsp] qlen slen. diamond's outfmt has both
    qcovhsp and scovhsp (10 columns total); blastp has no subject-coverage keyword, so
    its wider layout is one column narrower (9 total) and `scov` stays None.

    A file is distinguished purely by its column COUNT, never a version marker: old
    (4-col) and new (9/10-col) files coexist indefinitely in the same search_cache/,
    since storeDir never invalidates an existing cached file just because the code that
    would produce a new one changed.
    """
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.rstrip('\n').split('\t')
        hit = Hit(
            query_id=cols[0],
            target_id=cols[1],
            evalue=float(cols[2]),
            bitscore=float(cols[3]),
        )
        if len(cols) == 10:      # diamond wide: + length pident qcovhsp scovhsp qlen slen
            hit.length = int(cols[4])
            hit.pident = float(cols[5])
            hit.qcov = float(cols[6])
            hit.scov = float(cols[7])
            hit.qlen = int(cols[8])
            hit.slen = int(cols[9])
        elif len(cols) == 9:     # blastp wide: + length pident qcovhsp qlen slen (no scovhsp)
            hit.length = int(cols[4])
            hit.pident = float(cols[5])
            hit.qcov = float(cols[6])
            hit.qlen = int(cols[7])
            hit.slen = int(cols[8])
        elif len(cols) not in (4,):
            raise ValueError(
                f'unrecognized tabular hit row with {len(cols)} columns (expected 4 '
                f'narrow, 9 blastp-wide, or 10 diamond-wide): {line!r}')
        yield hit


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
