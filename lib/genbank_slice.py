"""GenBank slices of per-strain locus regions, for clinker (spec section 8).

Pure functions, no I/O at import time; bin/pangenome_island_gbk_slice.py
does the reading and writing.

A region arrives as a family_positions rank range. rank_entries() rebuilds
which gene (or TBLASTN rescue hit) sits at each rank with exactly the rule
bin/pangenome_build_family_positions.py uses: per strain, every
gene_positions row in file order, then every rescue_positions row in file
order, stable-sorted by (contig, start) and enumerated. A rank is therefore
(number of the strain's entries on contigs that sort before this contig) +
(index within this contig).
"""
from __future__ import annotations

import collections
import re

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

_PROTEIN_ID_RE = re.compile(r"(?:^|;)protein_id=([^;\n]+)")
_PARENT_RE = re.compile(r"(?:^|;)Parent=([^;\n]+)")
_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._-]")


def safe_name(name: str) -> str:
    """A file-name-safe form of a strain name."""
    return _UNSAFE_RE.sub("_", name) or "_"


def rank_entries(gene_rows, rescue_rows, wanted: set[tuple[str, str]]) -> dict:
    """{(strain, contig): [(rank, protein_id | None, family | None, start, end), ...]}
    for every (strain, contig) in `wanted`.

    `gene_rows` yields (strain, protein_id, contig, start, end) in
    gene_positions file order; `rescue_rows` yields (strain, family, contig,
    start) in rescue_positions file order. Entries of strains not in
    `wanted` are skipped; entries on other contigs of a wanted strain are
    only counted (they shift the rank offset)."""
    strains = {s for s, _ in wanted}
    counts: dict[tuple[str, str], int] = collections.Counter()
    kept: dict[tuple[str, str], list] = collections.defaultdict(list)
    for strain, pid, contig, start, end in gene_rows:
        if strain not in strains:
            continue
        counts[(strain, contig)] += 1
        if (strain, contig) in wanted:
            kept[(strain, contig)].append((start, pid, None, end))
    for strain, fam, contig, start in rescue_rows:
        if strain not in strains:
            continue
        counts[(strain, contig)] += 1
        if (strain, contig) in wanted:
            kept[(strain, contig)].append((start, None, fam, start))
    by_strain: dict[str, list[str]] = collections.defaultdict(list)
    for strain, contig in counts:
        by_strain[strain].append(contig)
    out = {}
    for strain, contig in wanted:
        offset = sum(counts[(strain, c)] for c in by_strain[strain] if c < contig)
        entries = sorted(kept.get((strain, contig), []), key=lambda e: e[0])
        out[(strain, contig)] = [(offset + i, pid, fam, start, end)
                                 for i, (start, pid, fam, end) in enumerate(entries)]
    return out


def cds_ids(attrs: str) -> list[str]:
    """IDs a GFF3 CDS row belongs to: protein_id= if present, else the
    comma-split Parent= (same rule as
    bin/pangenome_build_gene_positions.py::_resolve_ids)."""
    m = _PROTEIN_ID_RE.search(attrs)
    if m:
        return [m.group(1)]
    m = _PARENT_RE.search(attrs)
    if m:
        return [p for p in m.group(1).split(",") if p]
    return []


def cds_exons(gff3_lines, contigs: set[str], ids: set[str]) -> dict[str, tuple[int, list]]:
    """{protein_id: (strand, [(start, end), ...])} for CDS rows on `contigs`
    whose resolved ID is in `ids`. Coordinates are 1-based inclusive."""
    out: dict[str, list] = {}
    for line in gff3_lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or f[0] not in contigs:
            continue
        for pid in cds_ids(f[8]):
            if pid not in ids:
                continue
            entry = out.setdefault(pid, [-1 if f[6] == "-" else 1, []])
            entry[1].append((int(f[3]), int(f[4])))
    return {pid: (strand, sorted(ex)) for pid, (strand, ex) in out.items()}


def build_record(strain: str, contig: str, contig_seq: str, entries: list,
                 exons: dict, proteins: dict[str, str], family_of: dict[str, str],
                 label_of) -> tuple[SeqRecord, dict]:
    """One GenBank record for a region.

    `entries` are rank_entries() rows inside the region. Genes become a
    `gene` + `CDS` feature pair (the gene feature stops clinker's
    "Could not find parent gene" warning, spike problem 5); rescue hits have
    no gene model and are not drawn. `label_of(protein_id)` gives the
    /locus_tag. Returns (record, summary) with summary keys bp_start,
    bp_end, n_genes, n_rescue, n_missing, labels ({label: family})."""
    starts = [e[3] for e in entries]
    ends = [e[4] for e in entries]
    bp_start, bp_end = min(starts), max(ends)
    sub = contig_seq[bp_start - 1:bp_end]
    rec = SeqRecord(Seq(sub), id=safe_name(strain)[:16], name=safe_name(strain)[:16],
                    description=f"{strain} {contig}:{bp_start}-{bp_end}")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    labels: dict[str, str] = {}
    n_genes = n_rescue = n_missing = 0
    for _rank, pid, fam, _s, _e in entries:
        if pid is None:
            n_rescue += 1
            continue
        if pid not in exons:
            n_missing += 1
            continue
        strand, exs = exons[pid]
        parts = []
        for s, e in exs:
            ns, ne = max(0, s - bp_start), min(len(sub), e - bp_start + 1)
            if ne > ns:
                parts.append(FeatureLocation(ns, ne, strand=strand))
        if not parts:
            n_missing += 1
            continue
        if strand == -1:
            parts = parts[::-1]
        loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        label = label_of(pid)
        family = family_of.get(pid, "")
        quals = {"locus_tag": [label], "note": [f"family={family}"]}
        rec.features.append(SeqFeature(FeatureLocation(loc.start, loc.end, strand=strand),
                                       type="gene", qualifiers={"locus_tag": [label]}))
        cds_q = dict(quals)
        if proteins.get(pid):
            cds_q["translation"] = [proteins[pid].rstrip("*")]
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers=cds_q))
        labels[label] = family
        n_genes += 1
    return rec, {"bp_start": bp_start, "bp_end": bp_end, "n_genes": n_genes,
                 "n_rescue": n_rescue, "n_missing": n_missing, "labels": labels}
