#!/usr/bin/env python3
"""
Convert a BFD Fungi_BFD samples.csv (ASMID,SPECIES_IN,STRAIN,BIOPROJECT,
NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,
SPECIES,TRANSL_TABLE,LOCUSTAG) into the master samples CSV format expected by
bin/make_config.py (Species,Strain,Protein,DNA,Short,Lineage).

For each row, locates the funannotate output directory under
--annotation-dir/<Species>_<Strain>/predict_results/ (spaces -> underscores)
and reads its <dirname>.proteins.fa / <dirname>.scaffolds.fa. Rows with no
matching annotation directory are skipped (reported on stderr) rather than
guessing filenames that don't exist.

Short IDs are derived from the first 4 letters of the genus + first 4 of the
species epithet, disambiguated with a numeric suffix on collision.

Output is an intermediate sample pool, not a run config — write --output into
config_support/ (not configs/), then pass it as --samples to make_config.py to
produce the actual configs/<name>.csv run config.

Example (Chaetothyriales, broadened to --class-name so make_config.py has
sibling orders to auto-derive an outgroup from):

    python3 bin/convert_bfd_samples.py \\
        --input /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \\
        --annotation-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation \\
        --class-name Eurotiomycetes \\
        --link-dir data \\
        --output config_support/Chaetothyriales_samples.csv

    python3 bin/make_config.py \\
        --ingroup-taxon Chaetothyriales \\
        --max-per-outgroup-taxon 3 \\
        --seed 42 \\
        --output configs/Chaetothyriales.csv \\
        --data-dir data \\
        --samples config_support/Chaetothyriales_samples.csv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

LINEAGE_FIELDS = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']


def make_short(species, used):
    parts = species.split()
    genus = re.sub(r'[^A-Za-z]', '', parts[0])[:4] if parts else 'Sp'
    epithet = re.sub(r'[^A-Za-z]', '', parts[1])[:4] if len(parts) > 1 else ''
    base = (genus[:1].upper() + genus[1:].lower() + epithet.lower())[:8] or 'Sp'
    short = base
    n = 2
    while short in used:
        short = f"{base[:8 - len(str(n))]}{n}"
        n += 1
    used.add(short)
    return short


def find_annotation(annotation_dir, species, strain):
    dirname = f"{species}_{strain}".strip().rstrip('_').replace(' ', '_') if strain else species.replace(' ', '_')
    d = Path(annotation_dir) / dirname / 'predict_results'
    protein = d / f"{dirname}.proteins.fa"
    scaffolds = d / f"{dirname}.scaffolds.fa"
    if protein.exists() and scaffolds.exists():
        return dirname, protein, scaffolds
    return None, None, None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='BFD samples.csv')
    p.add_argument('--annotation-dir', required=True, help='BFD genome_annotation/ directory')
    p.add_argument('--order', help='Restrict to rows whose ORDER column matches this value')
    p.add_argument('--class-name', dest='class_name', help='Restrict to rows whose CLASS column matches this value')
    p.add_argument('--link-dir', help='If given, symlink matched Protein/DNA files into <link-dir>/pep and <link-dir>/dna')
    p.add_argument('--output', required=True, help='Destination master samples CSV (intermediate — put under config_support/, not configs/)')
    args = p.parse_args()

    used_shorts = set()
    out_rows = []
    missing = []

    with open(args.input, newline='') as fh:
        for row in csv.DictReader(fh):
            if args.order and row['ORDER'] != args.order:
                continue
            if args.class_name and row['CLASS'] != args.class_name:
                continue
            species = row['SPECIES'].strip()
            strain = row['STRAIN'].strip()
            dirname, protein, scaffolds = find_annotation(args.annotation_dir, species, strain)
            if dirname is None:
                missing.append(f"{species} ({strain}) [{row['ASMID']}]")
                continue

            lineage = ['Fungi'] + [row[f].strip() for f in LINEAGE_FIELDS if row[f].strip()]
            short = make_short(species, used_shorts)

            protein_name = protein.name
            dna_name = scaffolds.name
            if args.link_dir:
                pep_dir = Path(args.link_dir) / 'pep'
                dna_dir = Path(args.link_dir) / 'dna'
                pep_dir.mkdir(parents=True, exist_ok=True)
                dna_dir.mkdir(parents=True, exist_ok=True)
                pep_link = pep_dir / protein_name
                dna_link = dna_dir / dna_name
                if not pep_link.exists():
                    pep_link.symlink_to(protein.resolve())
                if not dna_link.exists():
                    dna_link.symlink_to(scaffolds.resolve())

            out_rows.append({
                'Species': species,
                'Strain': strain,
                'Protein': protein_name,
                'DNA': dna_name,
                'Short': short,
                'Lineage': ';'.join(lineage),
            })

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['Species', 'Strain', 'Protein', 'DNA', 'Short', 'Lineage'])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {args.output}: {len(out_rows)} samples", file=sys.stderr)
    if missing:
        print(f"Skipped {len(missing)} rows with no matching annotation directory under {args.annotation_dir}:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)


if __name__ == '__main__':
    main()
