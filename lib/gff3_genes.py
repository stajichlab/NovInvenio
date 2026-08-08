"""
GFF3-derived gene/mRNA chromosome + start-position lookup for the interactive
reports (novelties.html / core.html / losses.html).

Parses only `gene` and `mRNA` feature-type lines (GFF3 column 3), keyed by
their `ID=` attribute, into a single {feature_id: (chrom, start)} index.
Protein IDs in a proteome FASTA don't always match a GFF3 feature ID
directly -- lookup_gene_position() tries a few fallbacks (see its docstring)
before giving up.

Also home to gene_id_from_protein_id()/_GENE_ID_SUFFIX, moved here from
report_data.py: lookup_gene_position() depends on it, and report_data.py
depends on lookup_gene_position() (for the GFF3-lookup wiring), so the
function has to live on whichever side doesn't import the other -- this
module imports nothing from report_data.py, so it lives here.
report_data.py re-imports gene_id_from_protein_id from here for backward
compatibility (existing external callers, tests) rather than duplicating it.

No I/O happens at import time.
"""
import re
from pathlib import Path
from typing import Optional
from urllib.parse import unquote

from hits import open_input

# Trailing transcript/protein suffixes that separate a FungiDB gene ID from the
# per-transcript protein ID used in the proteome FASTAs.
#   Afu1g01620-T-p1      -> Afu1g01620
#   NCU00499-t26_1-p1    -> NCU00499
_GENE_ID_SUFFIX = re.compile(r'-[Tt][^-]*(-p\d+)?$|-p\d+$')

# Just the trailing protein-suffix half of _GENE_ID_SUFFIX -- stripping only
# this recovers a candidate *mRNA* ID (e.g. NCU00499-T1-p1 -> NCU00499-T1),
# for GFF3s that name the mRNA feature with the transcript suffix but never
# give the gene itself a matching top-level ID.
_PROTEIN_ID_SUFFIX = re.compile(r'-p\d+$')

# Feature types (GFF3 column 3) whose position we index. Order doesn't matter
# here -- gene-vs-mRNA precedence on an ID collision is resolved explicitly in
# parse_gff3() below, not by parse order.
_FEATURE_TYPES = ('gene', 'mRNA')

# Subdirectories searched for a GFF3 filename under --data_dir, mirroring
# main.nf's resolve_fa() for Protein/DNA (pep/proteins, dna/genome/scaffolds)
# plus a dedicated gff3/ subdir (CLAUDE.md's "Analysis Config CSV Format").
GFF3_SEARCH_SUBDIRS = ['gff3', 'pep', 'proteins', 'dna', 'genome', 'scaffolds']


def gene_id_from_protein_id(protein_id: str) -> str:
    """Strip transcript/protein suffixes to recover a likely source gene ID."""
    return _GENE_ID_SUFFIX.sub('', protein_id)


def _get_attr(attributes: str, key: str) -> Optional[str]:
    """Pull `key=value` out of a GFF3 column-9 `key=value;key=value;...` field."""
    prefix = key + '='
    for part in attributes.split(';'):
        part = part.strip()
        if part.startswith(prefix):
            value = part[len(prefix):]
            return unquote(value) if '%' in value else value
    return None


def parse_gff3(path) -> dict[str, tuple[str, int]]:
    """Parse a GFF3 file's gene/mRNA features into {feature_id: (chrom, start)}.

    Only `gene` and `mRNA` feature-type lines (column 3) are read. `start`
    (column 4) is used as-is -- GFF3 coordinates are already 1-based, no
    adjustment needed. Transparently handles `.gz` (see hits.open_input()).

    Gene-feature entries always win over an mRNA entry sharing the same ID
    (uncommon, but matches the spec's "gene wins on collision" rule) --
    resolved by keeping the two feature types in separate dicts and merging
    with gene last, rather than relying on file order.
    """
    gene_index: dict[str, tuple[str, int]] = {}
    mrna_index: dict[str, tuple[str, int]] = {}
    with open_input(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqid, _source, ftype, start, _end, _score, _strand, _phase, attrs = cols[:9]
            if ftype not in _FEATURE_TYPES:
                continue
            fid = _get_attr(attrs, 'ID')
            if not fid:
                continue
            try:
                start_i = int(start)
            except ValueError:
                continue
            if ftype == 'gene':
                gene_index[fid] = (seqid, start_i)
            else:
                mrna_index[fid] = (seqid, start_i)
    return {**mrna_index, **gene_index}


def lookup_gene_position(
    protein_id: str, index: dict[str, tuple[str, int]]
) -> Optional[tuple[str, int]]:
    """Resolve a proteome protein_id to a GFF3-derived (chrom, start).

    Annotation ID schemes are inconsistent across proteomes, so this tries,
    in order:
      a. protein_id itself, as a GFF3 feature ID (some funannotate-style
         GFF3s ID the mRNA/CDS directly with the protein ID).
      b. protein_id with only the trailing `-p\\d+` protein suffix stripped,
         as a candidate *mRNA* ID (e.g. NCU00499-T1-p1 -> NCU00499-T1) --
         covers GFF3s where the "gene" has no separately matching top-level
         ID but the mRNA feature keeps the transcript suffix.
      c. gene_id_from_protein_id(protein_id), as a candidate *gene* ID (the
         general funannotate/FungiDB-style transcript+protein suffix strip).
      d. No match -> None. Never raises.
    """
    if protein_id in index:
        return index[protein_id]

    mrna_candidate = _PROTEIN_ID_SUFFIX.sub('', protein_id)
    if mrna_candidate != protein_id and mrna_candidate in index:
        return index[mrna_candidate]

    gene_candidate = gene_id_from_protein_id(protein_id)
    if gene_candidate in index:
        return index[gene_candidate]

    return None


def load_gff3_index(
    path, cache: Optional[dict[str, dict[str, tuple[str, int]]]] = None
) -> dict[str, tuple[str, int]]:
    """parse_gff3(path), memoized in `cache` (keyed by resolved path string).

    Report payload builders can hold thousands of proteins per species but
    only a handful of distinct GFF3 files -- pass the same `cache` dict
    across every lookup for one report build so each file is parsed once.
    """
    key = str(path)
    if cache is None:
        return parse_gff3(path)
    if key not in cache:
        cache[key] = parse_gff3(path)
    return cache[key]


def resolve_gff3_path(data_dir, filename: str):
    """Resolve a GFF3 basename under data_dir, flat then GFF3_SEARCH_SUBDIRS.

    Mirrors main.nf's resolve_fa() for Protein/DNA columns. Returns None
    (never raises) when data_dir/filename is empty or nothing is found --
    a missing/unresolvable GFF3 just means no chrom/start data for that
    species, not a hard error (GFF3 is an optional per-row config column).
    """
    if not data_dir or not filename:
        return None
    base = Path(data_dir)
    candidates = [base / filename] + [base / sub / filename for sub in GFF3_SEARCH_SUBDIRS]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def resolve_gff3_paths(samples, data_dir) -> dict[str, str]:
    """{Sample.short: resolved GFF3 path} for every sample with a GFF3 column
    value that resolves under data_dir. Samples with no GFF3 value, or whose
    GFF3 file can't be found, are simply omitted -- never an error.
    """
    if not data_dir:
        return {}
    paths: dict[str, str] = {}
    for sample in samples:
        gff3 = getattr(sample, 'gff3', '') or ''
        if not gff3:
            continue
        resolved = resolve_gff3_path(data_dir, gff3)
        if resolved:
            paths[sample.short] = str(resolved)
    return paths
