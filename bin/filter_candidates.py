#!/usr/bin/env python3
"""
Filter presence_matrix.tsv (or .function.tsv) to lineage-specific candidates:
  - present in >= INGROUP_MIN fraction of ingroup proteomes
  - absent (0) in ALL outgroup proteomes

Ingroup/outgroup membership is read from the analysis config CSV (GROUP column).
"""
import argparse
import csv
import sys


def load_groups(config_csv):
    """Return (ingroup_ids, outgroup_ids) as lists of Short IDs."""
    ingroup, outgroup = [], []
    with open(config_csv, newline='') as fh:
        for row in csv.DictReader(fh):
            short = row['Short'].strip()
            if row['GROUP'].strip().upper() == 'IN':
                ingroup.append(short)
            else:
                outgroup.append(short)
    return ingroup, outgroup


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.tsv or presence_matrix.function.tsv')
    ap.add_argument('--config', required=True,
                    help='Analysis config CSV (GROUP, Short columns required)')
    ap.add_argument('--ingroup_min', type=float, default=1.0,
                    help='Minimum fraction of ingroup proteomes that must have a hit '
                         '(default: 1.0 = 100%%)')
    ap.add_argument('--output', default='-',
                    help='Output TSV (default: stdout)')
    args = ap.parse_args()

    if not (0.0 < args.ingroup_min <= 1.0):
        sys.exit('--ingroup_min must be in (0, 1]')

    ingroup_ids, outgroup_ids = load_groups(args.config)

    fout = open(args.output, 'w', newline='') if args.output != '-' else sys.stdout

    with open(args.matrix, newline='') as fin:
        reader = csv.DictReader(fin, delimiter='\t')

        # Validate that the column IDs from the config exist in the matrix
        header = reader.fieldnames
        missing = [c for c in ingroup_ids + outgroup_ids if c not in header]
        if missing:
            sys.exit(f'ERROR: columns not found in matrix: {missing}\n'
                     f'Matrix columns: {header}')

        writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()

        n_total = n_kept = 0
        for row in reader:
            n_total += 1

            ingroup_present = sum(int(row[c]) for c in ingroup_ids)
            ingroup_frac = ingroup_present / len(ingroup_ids)

            outgroup_present = sum(int(row[c]) for c in outgroup_ids)

            if ingroup_frac >= args.ingroup_min and outgroup_present == 0:
                writer.writerow(row)
                n_kept += 1

    if fout is not sys.stdout:
        fout.close()

    print(
        f'Kept {n_kept} / {n_total} rows '
        f'(ingroup >= {args.ingroup_min:.0%}, outgroup = 0)',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
