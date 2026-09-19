#!/usr/bin/env python3
"""Per-strain protein_id -> genomic position lookup, built directly from
each strain's own GFF3 CDS records -- the foundation for pair-classification
physical-linkage checks and for checking whether an optional captain-gene
(e.g. fungal Starship DUF3435) hit sits near a co-occurring family pair.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/build_gene_positions.py); `config_parser` is imported directly from
this repo's own lib/ (no cross-repo lookup needed).

Deliberately keyed on PROTEIN ID (the CDS `protein_id=` attribute), not gene
ID: family membership (the tier-1 cluster TSV) is already defined by protein
ID, and every downstream consumer needs "where is this specific protein",
not an intermediate gene locus. A multi-exon CDS repeats the same
protein_id across several rows; this collapses those to one
(min(start), max(end)) span per protein_id, per GFF3 spec's start<=end
regardless of strand.

GFF3 dialect handling (added after the funannotate-annotated Coccidioides
study found 0/535 of its GFF3s carry `protein_id=` at all -- see
NovInvenio_Investigations' 2026-09-15 onboarding spec for the full
diagnosis): `protein_id=` is preferred when present (NCBI-style GFF3,
unchanged behavior); otherwise this falls back to the CDS `Parent=`
attribute (funannotate/generic style), split on comma for the rare case of
a CDS shared by multiple transcripts (GFF3 spec permits `Parent=A,B`).
Neither attribute matching is proof the resolved ID is the one the strain's
own protein FASTA actually uses -- a GFF3 syntactically parsing is not the
same as it being the RIGHT dialect for this strain (a JGI/Ensembl-style
GFF3 can carry a transcript ID distinct from its paired FASTA's protein
IDs) -- so every resolved ID is cross-checked against that strain's actual
protein FASTA headers (`--protein_dir`), and a strain whose resolved IDs
mostly DON'T match its own FASTA hard-errors instead of silently emitting
near-empty position data that collapses everything downstream to
`insufficient_data`.

Usage:
  pangenome_build_gene_positions.py --config config.csv --gff3_dir data_dir/gff3 \\
      --protein_dir data_dir/pep --output gene_positions.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from compressed_io import open_maybe_compressed_write  # noqa: E402

_PROTEIN_ID_RE = re.compile(r"(?:^|;)protein_id=([^;\n]+)")
_PARENT_RE = re.compile(r"(?:^|;)Parent=([^;\n]+)")

# Below this fraction of a strain's proteins actually resolving to a FASTA-
# matching position, treat it as a genuine dialect mismatch and abort rather
# than emit position data that's mostly wrong/missing.
HARD_ERROR_BELOW_MATCH_FRACTION = 0.50
# Above this fraction unresolved (but still above the hard-error floor),
# warn -- a handful of stray annotation-tool quirks, not a structural
# problem, but worth flagging.
WARN_ABOVE_UNRESOLVED_FRACTION = 0.02


def _resolve_ids(attrs: str) -> list[str]:
    """Return the CDS attribute ID(s) this row should be keyed on: prefer
    `protein_id=` (NCBI-style GFF3); otherwise fall back to `Parent=`
    (funannotate/generic style), split on comma since GFF3 permits a CDS
    shared by multiple transcripts (`Parent=A,B`). Returns [] if neither
    attribute is present on this row."""
    m = _PROTEIN_ID_RE.search(attrs)
    if m:
        return [m.group(1)]
    m = _PARENT_RE.search(attrs)
    if m:
        return [pid for pid in m.group(1).split(",") if pid]
    return []


def parse_gff3_protein_positions(gff3_path: str | Path) -> dict[str, tuple[str, int, int]]:
    """Return {protein_id: (contig, start, end)} from a GFF3's CDS records,
    collapsing multi-exon CDS rows that share a protein_id to one
    (min start, max end) span. `protein_id` here is whichever ID
    `_resolve_ids` resolved for that row (protein_id= or Parent= fallback)
    -- this is a raw parse, not yet cross-checked against any protein
    FASTA; see `resolve_positions_for_strain` for that."""
    positions: dict[str, list] = {}
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            contig, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            for protein_id in _resolve_ids(attrs):
                if protein_id not in positions:
                    positions[protein_id] = [contig, start, end]
                else:
                    entry = positions[protein_id]
                    entry[1] = min(entry[1], start)
                    entry[2] = max(entry[2], end)
    return {pid: (contig, start, end) for pid, (contig, start, end) in positions.items()}


def load_protein_ids(protein_fasta_path: str | Path) -> set[str]:
    """{protein_id, ...} parsed from a protein FASTA's header first token
    (before any whitespace) -- this repo's standard FASTA-header-as-ID
    convention, matching how CONCAT_PROTEOMES/CLUSTER_TIER1 key clusters."""
    ids: set[str] = set()
    with open(protein_fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0].rstrip())
    return ids


def resolve_positions_for_strain(
    gff3_path: str | Path,
    protein_ids: set[str],
    short: str,
    hard_error_below: float = HARD_ERROR_BELOW_MATCH_FRACTION,
    warn_above_unresolved: float = WARN_ABOVE_UNRESOLVED_FRACTION,
) -> dict[str, tuple[str, int, int]]:
    """Parse `gff3_path`, then filter to only the IDs that actually appear
    in `protein_ids` (that strain's real protein FASTA headers) -- this is
    what turns "the GFF3 attribute parsed" into "the GFF3 attribute is
    actually this strain's real dialect". Hard-errors (aborts the whole
    run) if the matched fraction falls below `hard_error_below`: this is
    the actual "wrong dialect entirely" signal, not a per-gene data quality
    issue. Warns (but continues) if unresolved fraction exceeds
    `warn_above_unresolved` but stays above the hard-error floor."""
    raw_positions = parse_gff3_protein_positions(gff3_path)
    matched = {pid: pos for pid, pos in raw_positions.items() if pid in protein_ids}
    total = len(protein_ids)
    matched_fraction = (len(matched) / total) if total else 0.0

    if matched_fraction < hard_error_below:
        print(
            f"ERROR: {short}: only {len(matched)}/{total} proteins "
            f"({matched_fraction:.1%}) resolved a GFF3 position matching this "
            f"strain's protein FASTA -- checked protein_id= and Parent= "
            f"fallback, neither corresponds to this strain's actual protein "
            f"IDs. This looks like a genuine GFF3/FASTA dialect mismatch, not "
            f"a handful of odd genes -- aborting rather than emit mostly-empty "
            f"position data that would silently collapse downstream "
            f"co-occurrence calls to insufficient_data.",
            file=sys.stderr,
        )
        sys.exit(1)
    elif matched_fraction < (1 - warn_above_unresolved):
        print(
            f"WARNING: {short}: {total - len(matched)}/{total} proteins "
            f"({1 - matched_fraction:.1%}) had no GFF3 position resolvable "
            f"against this strain's protein FASTA.",
            file=sys.stderr,
        )
    return matched


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--protein_dir", required=True,
                    help="Directory holding each strain's protein FASTA (Sample.protein) "
                         "-- used to cross-check resolved GFF3 IDs against real protein IDs.")
    ap.add_argument("--groups", default="IN,OUT")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from config_parser import parse_config  # noqa: E402

    wanted_groups = set(args.groups.split(","))
    samples = [s for s in parse_config(args.config) if s.group in wanted_groups]
    gff3_dir = Path(args.gff3_dir)
    protein_dir = Path(args.protein_dir)

    n_strains_ok, n_strains_missing_gff3 = 0, 0
    with open_maybe_compressed_write(args.output) as out:
        out.write("Short\tprotein_id\tcontig\tstart\tend\n")
        for s in samples:
            if not s.gff3:
                n_strains_missing_gff3 += 1
                continue
            gff3_path = gff3_dir / s.gff3
            if not gff3_path.exists():
                print(f"WARNING: {gff3_path} not found for {s.short}, skipping", file=sys.stderr)
                n_strains_missing_gff3 += 1
                continue
            protein_fasta_path = protein_dir / s.protein
            if not protein_fasta_path.exists():
                print(f"ERROR: {protein_fasta_path} not found for {s.short} -- "
                      f"cannot cross-check GFF3 IDs without this strain's protein "
                      f"FASTA", file=sys.stderr)
                sys.exit(1)
            protein_ids = load_protein_ids(protein_fasta_path)
            positions = resolve_positions_for_strain(gff3_path, protein_ids, s.short)
            for protein_id, (contig, start, end) in positions.items():
                out.write(f"{s.short}\t{protein_id}\t{contig}\t{start}\t{end}\n")
            n_strains_ok += 1

    print(
        f"build_gene_positions: {n_strains_ok} strains parsed, "
        f"{n_strains_missing_gff3} skipped (no GFF3)", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
