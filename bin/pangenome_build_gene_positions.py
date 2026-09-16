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

Usage:
  pangenome_build_gene_positions.py --config config.csv --gff3_dir data_dir/gff3 \\
      --output gene_positions.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_PROTEIN_ID_RE = re.compile(r"(?:^|;)protein_id=([^;\n]+)")


def parse_gff3_protein_positions(gff3_path: str | Path) -> dict[str, tuple[str, int, int]]:
    """Return {protein_id: (contig, start, end)} from a GFF3's CDS records,
    collapsing multi-exon CDS rows that share a protein_id to one
    (min start, max end) span."""
    positions: dict[str, list] = {}
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            contig, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            m = _PROTEIN_ID_RE.search(attrs)
            if not m:
                continue
            protein_id = m.group(1)
            if protein_id not in positions:
                positions[protein_id] = [contig, start, end]
            else:
                entry = positions[protein_id]
                entry[1] = min(entry[1], start)
                entry[2] = max(entry[2], end)
    return {pid: (contig, start, end) for pid, (contig, start, end) in positions.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--gff3_dir", required=True)
    ap.add_argument("--groups", default="IN,OUT")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
    from config_parser import parse_config  # noqa: E402

    wanted_groups = set(args.groups.split(","))
    samples = [s for s in parse_config(args.config) if s.group in wanted_groups]
    gff3_dir = Path(args.gff3_dir)

    n_strains_ok, n_strains_missing_gff3 = 0, 0
    with open(args.output, "w") as out:
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
            positions = parse_gff3_protein_positions(gff3_path)
            for protein_id, (contig, start, end) in positions.items():
                out.write(f"{s.short}\t{protein_id}\t{contig}\t{start}\t{end}\n")
            n_strains_ok += 1

    print(
        f"build_gene_positions: {n_strains_ok} strains parsed, "
        f"{n_strains_missing_gff3} skipped (no GFF3)", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
