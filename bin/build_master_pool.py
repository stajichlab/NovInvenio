#!/usr/bin/env python3
"""
Build the master species pool consumed by bin/build_targeted_configs.py,
from a BFD samples.csv (ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,
BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,
TRANSL_TABLE,LOCUSTAG) plus its matching repr_assignments.tsv
(_reuse_assignments/repr_assignments.tsv: out,species,is_representative,
representative_out,ani_to_representative,reuse_eligible).

One row per species (never per assembly): the assembly whose dirname
matches that species' repr_assignments.tsv is_representative=True `out`
value is the only one kept. Unlike bin/convert_bfd_samples.py, Lineage
keeps a fixed 7-field width (PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;
FAMILY;GENUS, '' for an unrecorded rank, never dropped) and
Protein/DNA are resolved to absolute paths, not basenames -- see
docs/superpowers/specs/2026-09-05-config-builder-design.md ("Master pool").

Usage:
    bin/build_master_pool.py \\
        --bfd-samples /path/to/Fungi_BFD/samples.csv \\
        --annotation-dir /path/to/Fungi_BFD/genome_annotation \\
        --repr-assignments /path/to/Fungi_BFD_runs/genome_annotation/_reuse_assignments/repr_assignments.tsv \\
        --output config_support/master_pool.csv
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from lineage import RANK_NAMES  # noqa: E402
from master_pool import MASTER_POOL_FIELDS, load_representative_picks  # noqa: E402

LINEAGE_FIELDS = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']
assert LINEAGE_FIELDS == RANK_NAMES  # single source of truth check


def find_annotation(annotation_dir, dirname):
    d = Path(annotation_dir) / dirname / 'predict_results'
    protein = d / f"{dirname}.proteins.fa"
    scaffolds = d / f"{dirname}.scaffolds.fa"
    if protein.exists() and scaffolds.exists():
        return protein.resolve(), scaffolds.resolve()
    return None, None


def render_master_pool(bfd_samples_path, annotation_dir, repr_assignments_path) -> list[dict]:
    with open(bfd_samples_path, newline='') as fh:
        bfd_rows = list(csv.DictReader(fh))

    species_set = {row['SPECIES'].strip() for row in bfd_rows}
    picks = load_representative_picks(repr_assignments_path, expected_species=species_set)

    by_dirname = {
        f"{row['SPECIES'].strip()}_{row['STRAIN'].strip()}".rstrip('_').replace(' ', '_'): row
        for row in bfd_rows
    }

    out_rows = []
    missing = []
    for species, dirname in sorted(picks.items()):
        row = by_dirname.get(dirname)
        if row is None:
            missing.append(f"{species}: representative dirname {dirname!r} not found in {bfd_samples_path}")
            continue
        protein, scaffolds = find_annotation(annotation_dir, dirname)
        if protein is None:
            missing.append(f"{species}: no predict_results/{dirname}.proteins.fa+scaffolds.fa under {annotation_dir}")
            continue
        lineage = [row[f].strip() for f in LINEAGE_FIELDS]
        out_rows.append({
            'Species': species,
            'Strain': row['STRAIN'].strip(),
            'ProteinPath': str(protein),
            'DNAPath': str(scaffolds),
            'Lineage': ';'.join(lineage),
            'NCBI_TaxID': row['NCBI_TAXONID'].strip(),
        })

    if missing:
        raise SystemExit("Could not build master pool rows:\n  " + "\n  ".join(missing))
    return out_rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bfd-samples', required=True, help='Fungi_BFD/samples.csv')
    p.add_argument('--annotation-dir', required=True, help='Fungi_BFD/genome_annotation directory')
    p.add_argument('--repr-assignments', required=True, help='_reuse_assignments/repr_assignments.tsv')
    p.add_argument('--output', required=True, help='Destination master pool CSV')
    args = p.parse_args()

    rows = render_master_pool(args.bfd_samples, args.annotation_dir, args.repr_assignments)

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=MASTER_POOL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output}: {len(rows)} species", file=sys.stderr)


if __name__ == '__main__':
    main()
