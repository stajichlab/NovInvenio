#!/usr/bin/env python3
"""
Convert 1KFG (1000 Fungal Genomes, JGI MycoCosm) proteomes under
/bigdata/stajichlab/shared/projects/1KFG/genomes/final_combine/ into the master
samples CSV format expected by bin/make_config.py (Species,Strain,Protein,DNA,
Short,Lineage) -- same downstream contract as bin/convert_bfd_samples.py, for a
different upstream genome source.

1KFG's final_combine/ is a flat directory of JGI-style filenames
(<Genus>_<species>_<strain-or-portal-id>.<version>.aa.fasta under pep/, the
same basename with .fasta instead of .aa.fasta under DNA/) with no taxonomy of
its own. Taxonomy is recovered by joining each 1KFG species (parsed
genus/epithet from the filename) against Fungi_BFD's samples.csv on
GENUS/SPECIES -- BFD is NCBI-assembly-based and 1KFG is JGI-based, so the two
collections only partially overlap by exact genus+epithet match (~70% in
practice); unmatched species are skipped (reported on stderr) rather than
guessing a lineage.

Output is an intermediate sample pool, not a run config -- write --output into
config_support/ (not configs/), then pass it as --samples to make_config.py to
produce the actual configs/<name>.csv run config. See bin/make_config.py's
--max-per-outgroup-taxon/--random/--seed for stratified random sampling within
target phyla/subphyla once this pool exists.

Example:

    python3 bin/convert_1kfg_samples.py \\
        --pep-dir /bigdata/stajichlab/shared/projects/1KFG/genomes/final_combine/pep \\
        --dna-dir /bigdata/stajichlab/shared/projects/1KFG/genomes/final_combine/DNA \\
        --bfd-samples /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \\
        --link-dir data \\
        --output config_support/1KFG_samples.csv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

LINEAGE_FIELDS = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']


def load_bfd_taxonomy(bfd_samples_path):
    """Return dict: (genus.lower(), species_epithet.lower()) -> BFD row dict.

    BFD's SPECIES column holds the full binomial ("Genus species"), not just the
    epithet -- split it apart to match against 1KFG's parsed (genus, epithet).
    """
    taxonomy = {}
    with open(bfd_samples_path, newline='') as fh:
        for row in csv.DictReader(fh):
            genus = row['GENUS'].strip()
            species_full = row['SPECIES'].strip()
            toks = species_full.split(' ', 1)
            epithet = toks[1] if len(toks) > 1 else ''
            key = (genus.lower(), epithet.lower())
            if key and key not in taxonomy:
                taxonomy[key] = row
    return taxonomy


def parse_1kfg_basename(basename):
    """Return (genus, species_epithet, display_species, strain) from a 1KFG
    filename basename like 'Chaetomium_globosum_CBS_148.51' or
    'Neurospora_crassa_OR74A' or 'Acidomyces_richmondensis_BFW.Aciri1_iso'.

    Convention: first two underscore-separated tokens are Genus_species; the
    strain/portal-id tokens are the remaining underscore-separated tokens, with
    the version/portal-code suffix after the first '.' stripped from the
    display strain (kept in the on-disk filename regardless -- Protein/DNA
    columns always use the real basename).
    """
    parts = basename.split('_')
    if len(parts) < 2:
        return None
    genus = parts[0]
    epithet_raw = parts[1].split('.')[0]
    strain_parts = parts[2:]
    strain = '_'.join(strain_parts).split('.')[0] if strain_parts else ''
    display_species = f"{genus} {epithet_raw}"
    return genus, epithet_raw, display_species, strain


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


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pep-dir', required=True, dest='pep_dir',
                   help='1KFG final_combine/pep directory (plain .aa.fasta files)')
    p.add_argument('--dna-dir', required=True, dest='dna_dir',
                   help='1KFG final_combine/DNA directory (plain .fasta files)')
    p.add_argument('--bfd-samples', required=True, dest='bfd_samples',
                   help='Fungi_BFD samples.csv, for taxonomy join by genus+species')
    p.add_argument('--phylum', help='Restrict to rows whose BFD PHYLUM matches this value')
    p.add_argument('--subphylum', help='Restrict to rows whose BFD SUBPHYLUM matches this value')
    p.add_argument('--class-name', dest='class_name', help='Restrict to rows whose BFD CLASS matches this value')
    p.add_argument('--order', help='Restrict to rows whose BFD ORDER matches this value')
    p.add_argument('--link-dir', help='If given, symlink matched Protein/DNA files into <link-dir>/pep and <link-dir>/dna')
    p.add_argument('--output', required=True, help='Destination master samples CSV (intermediate -- put under config_support/, not configs/)')
    args = p.parse_args()

    taxonomy = load_bfd_taxonomy(args.bfd_samples)

    pep_dir = Path(args.pep_dir)
    dna_dir = Path(args.dna_dir)
    pep_files = sorted(f for f in pep_dir.glob('*.aa.fasta'))

    used_shorts = set()
    out_rows = []
    unmatched = []
    filtered_out = 0

    for pep_path in pep_files:
        basename = pep_path.name[:-len('.aa.fasta')]
        parsed = parse_1kfg_basename(basename)
        if parsed is None:
            unmatched.append(f"{basename} [unparseable filename]")
            continue
        genus, epithet, display_species, strain = parsed
        key = (genus.lower(), epithet.lower())
        row = taxonomy.get(key)
        if row is None:
            unmatched.append(f"{display_species} ({basename}) [no BFD taxonomy match]")
            continue

        if args.phylum and row['PHYLUM'] != args.phylum:
            filtered_out += 1
            continue
        if args.subphylum and row['SUBPHYLUM'] != args.subphylum:
            filtered_out += 1
            continue
        if args.class_name and row['CLASS'] != args.class_name:
            filtered_out += 1
            continue
        if args.order and row['ORDER'] != args.order:
            filtered_out += 1
            continue

        dna_path = dna_dir / f"{basename}.fasta"
        if not dna_path.exists():
            unmatched.append(f"{display_species} ({basename}) [no matching DNA file]")
            continue

        lineage = ['Fungi'] + [row[f].strip() for f in LINEAGE_FIELDS if row[f].strip()]
        short = make_short(display_species, used_shorts)

        if args.link_dir:
            out_pep_dir = Path(args.link_dir) / 'pep'
            out_dna_dir = Path(args.link_dir) / 'dna'
            out_pep_dir.mkdir(parents=True, exist_ok=True)
            out_dna_dir.mkdir(parents=True, exist_ok=True)
            pep_link = out_pep_dir / pep_path.name
            dna_link = out_dna_dir / dna_path.name
            if not pep_link.exists():
                pep_link.symlink_to(pep_path.resolve())
            if not dna_link.exists():
                dna_link.symlink_to(dna_path.resolve())

        out_rows.append({
            'Species': display_species,
            'Strain': strain,
            'Protein': pep_path.name,
            'DNA': dna_path.name,
            'Short': short,
            'Lineage': ';'.join(lineage),
        })

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['Species', 'Strain', 'Protein', 'DNA', 'Short', 'Lineage'])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {args.output}: {len(out_rows)} samples "
          f"({len(unmatched)} unmatched, {filtered_out} filtered by taxon restriction)",
          file=sys.stderr)
    if unmatched:
        print("Unmatched (no BFD taxonomy join, unparseable filename, or missing DNA):", file=sys.stderr)
        for m in unmatched:
            print(f"  {m}", file=sys.stderr)


if __name__ == '__main__':
    main()
