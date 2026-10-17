"""Input schemas for the optional Starship / SM-cluster overlap step of pangenome.nf.

Two user-supplied tables (precomputed outside the pipeline, e.g. starfish element calls
and antiSMASH regions) are matched to a finished pangenome run:

Starship element table (``--pangenome_starship_annotation``), coordinate form
    element_id   unique element name
    strain       the strain's ``Short`` ID as in the samplesheet (``Short`` is accepted
                 as an alias column name)
    contig       contig/chromosome name as in the strain's GFF3 (= gene_positions.tsv)
    start, end   1-based, inclusive, integers, start <= end
    label        optional free text (e.g. "Osiris h4")
    confidence   optional; only the value ``high`` is treated as high confidence

Starship element table, gene-level form (detected when the header has ``gene_id`` and
no ``start``/``end``): one row per gene, columns ``strain``, ``gene_id``, ``element_id``
(+ optional ``label``, ``confidence``). ``gene_id`` is the protein ID as in the strain's
protein FASTA / gene_positions.tsv. The element's span is derived later from its genes'
positions (see ``derive_gene_level_spans``).

SM cluster table (``--pangenome_sm_clusters``)
    cluster_id   unique cluster name
    strain, contig, start, end   as above
    type         cluster type (e.g. antiSMASH region products joined by '+')
    genes        optional, comma-separated protein IDs of the cluster genes; when empty,
                 the cluster genes are the strain's genes in gene_positions that lie
                 fully inside [start, end]
    core_genes   optional, comma-separated protein IDs of the core (biosynthetic) genes
    contig_edge, knowncluster_top_hit, knowncluster_description,
    knowncluster_similarity      optional, passed through to sm_clusters.tsv.gz

Files are tab-separated, or comma-separated when the name ends in ``.csv`` /
``.csv.gz``. ``.gz`` and ``.zst`` are read transparently. Blank lines and lines starting
with ``#`` are skipped.

The same checks are done in Groovy at workflow start (pangenome.nf,
``validate_region_table``); this module repeats them inside the processes and is what the
unit tests cover.
"""
from __future__ import annotations

import csv
from pathlib import Path

from compressed_io import open_maybe_compressed

ELEMENT_REQUIRED = ("element_id", "strain", "contig", "start", "end")
ELEMENT_GENE_REQUIRED = ("strain", "gene_id", "element_id")
SM_REQUIRED = ("cluster_id", "strain", "contig", "start", "end", "type")
SM_PASSTHROUGH = ("contig_edge", "knowncluster_top_hit", "knowncluster_description",
                  "knowncluster_similarity")
STRAIN_ALIASES = ("strain", "Short")
MAX_REPORTED_ERRORS = 20


class TableError(ValueError):
    """A user-supplied table does not match its schema."""


def _delimiter(path: Path) -> str:
    name = path.name
    for suf in (".gz", ".zst"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return "," if name.endswith(".csv") else "\t"


def read_table(path: str | Path) -> tuple[list[str], list[dict]]:
    """Return (header, rows) of a TSV/CSV table, skipping blank and '#' lines."""
    path = Path(path)
    with open_maybe_compressed(path) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise TableError(f"{path}: empty table (no header line)")
    reader = csv.DictReader(lines, delimiter=_delimiter(path))
    header = [h.strip() for h in (reader.fieldnames or [])]
    reader.fieldnames = header
    rows = [{k: (v or "").strip() for k, v in r.items() if k is not None} for r in reader]
    return header, rows


def _strain_column(header: list[str], path: Path) -> str:
    for c in STRAIN_ALIASES:
        if c in header:
            return c
    raise TableError(f"{path}: missing column 'strain' (or its alias 'Short'); found {header}")


def _int(value: str, what: str, errors: list[str]) -> int | None:
    try:
        v = int(value)
    except (TypeError, ValueError):
        errors.append(f"{what}: '{value}' is not an integer")
        return None
    if v < 1:
        errors.append(f"{what}: {v} < 1 (coordinates are 1-based)")
        return None
    return v


def _raise(path: Path, errors: list[str]) -> None:
    if errors:
        shown = errors[:MAX_REPORTED_ERRORS]
        more = len(errors) - len(shown)
        raise TableError(f"{path}: {len(errors)} problem(s):\n  " + "\n  ".join(shown)
                         + (f"\n  ... and {more} more" if more else ""))


def _check_strain(strain: str, known: set[str] | None, where: str, errors: list[str]) -> None:
    if not strain:
        errors.append(f"{where}: empty strain")
    elif known is not None and strain not in known:
        errors.append(f"{where}: strain '{strain}' is not a Short ID of this run's samplesheet")


def is_gene_level(header: list[str]) -> bool:
    return "gene_id" in header and "start" not in header and "end" not in header


def load_elements(path: str | Path, known_strains: set[str] | None = None) -> list[dict]:
    """Validate and load a Starship element table (either form).

    Returns one dict per element with keys element_id, label, confidence, Short, contig,
    start, end, genes. For the gene-level form contig/start/end are None and genes is the
    list of gene IDs; for the coordinate form genes is None.
    """
    path = Path(path)
    header, rows = read_table(path)
    scol = _strain_column(header, path)
    errors: list[str] = []
    if is_gene_level(header):
        missing = [c for c in ELEMENT_GENE_REQUIRED if c != "strain" and c not in header]
        if missing:
            raise TableError(f"{path}: gene-level element table missing column(s) {missing}")
        els: dict[str, dict] = {}
        for i, r in enumerate(rows, start=2):
            where = f"line {i}"
            eid, gid, strain = r.get("element_id", ""), r.get("gene_id", ""), r.get(scol, "")
            if not eid:
                errors.append(f"{where}: empty element_id")
            if not gid:
                errors.append(f"{where}: empty gene_id")
            _check_strain(strain, known_strains, where, errors)
            e = els.setdefault(eid, {"element_id": eid, "label": r.get("label", ""),
                                     "confidence": r.get("confidence", ""), "Short": strain,
                                     "contig": None, "start": None, "end": None, "genes": []})
            if e["Short"] != strain:
                errors.append(f"{where}: element '{eid}' has genes in two strains "
                              f"('{e['Short']}', '{strain}')")
            e["genes"].append(gid)
        _raise(path, errors)
        if not els:
            raise TableError(f"{path}: no element rows")
        return list(els.values())

    missing = [c for c in ELEMENT_REQUIRED if c != "strain" and c not in header]
    if missing:
        raise TableError(f"{path}: missing column(s) {missing}; required: "
                         f"{list(ELEMENT_REQUIRED)} (or gene-level: {list(ELEMENT_GENE_REQUIRED)})")
    out, seen = [], set()
    for i, r in enumerate(rows, start=2):
        where = f"line {i}"
        eid = r.get("element_id", "")
        if not eid:
            errors.append(f"{where}: empty element_id")
        elif eid in seen:
            errors.append(f"{where}: duplicate element_id '{eid}'")
        seen.add(eid)
        strain = r.get(scol, "")
        _check_strain(strain, known_strains, where, errors)
        if not r.get("contig"):
            errors.append(f"{where}: empty contig")
        s = _int(r.get("start", ""), f"{where} start", errors)
        e = _int(r.get("end", ""), f"{where} end", errors)
        if s is not None and e is not None and s > e:
            errors.append(f"{where}: start {s} > end {e}")
        out.append({"element_id": eid, "label": r.get("label", ""),
                    "confidence": r.get("confidence", ""), "Short": strain,
                    "contig": r.get("contig", ""), "start": s, "end": e, "genes": None})
    _raise(path, errors)
    if not out:
        raise TableError(f"{path}: no element rows")
    return out


def _split_ids(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()] if value else []


def load_sm_clusters(path: str | Path, known_strains: set[str] | None = None) -> list[dict]:
    """Validate and load an SM cluster table. Returns one dict per cluster."""
    path = Path(path)
    header, rows = read_table(path)
    scol = _strain_column(header, path)
    missing = [c for c in SM_REQUIRED if c != "strain" and c not in header]
    if missing:
        raise TableError(f"{path}: missing column(s) {missing}; required: {list(SM_REQUIRED)}")
    errors: list[str] = []
    out, seen = [], set()
    for i, r in enumerate(rows, start=2):
        where = f"line {i}"
        cid = r.get("cluster_id", "")
        if not cid:
            errors.append(f"{where}: empty cluster_id")
        elif cid in seen:
            errors.append(f"{where}: duplicate cluster_id '{cid}'")
        seen.add(cid)
        strain = r.get(scol, "")
        _check_strain(strain, known_strains, where, errors)
        if not r.get("contig"):
            errors.append(f"{where}: empty contig")
        if not r.get("type"):
            errors.append(f"{where}: empty type")
        s = _int(r.get("start", ""), f"{where} start", errors)
        e = _int(r.get("end", ""), f"{where} end", errors)
        if s is not None and e is not None and s > e:
            errors.append(f"{where}: start {s} > end {e}")
        # Keep empty slots of an explicit gene list out; an ID that is not in
        # gene_positions is kept and reported as unplaced by the caller.
        c = {"cluster_id": cid, "Short": strain, "contig": r.get("contig", ""), "start": s, "end": e,
             "type": r.get("type", ""), "genes": _split_ids(r.get("genes", "")),
             "core_genes": _split_ids(r.get("core_genes", ""))}
        for k in SM_PASSTHROUGH:
            c[k] = r.get(k, "")
        out.append(c)
    _raise(path, errors)
    if not out:
        raise TableError(f"{path}: no cluster rows")
    return out


def samplesheet_strains(path: str | Path, groups: list[str] | None = None) -> set[str]:
    """Short IDs of a pangenome samplesheet (optionally only rows whose GROUP is in groups)."""
    _, rows = read_table(Path(path))
    return {r["Short"].strip() for r in rows
            if r.get("Short") and (groups is None or (r.get("GROUP") or "").strip() in groups)}
