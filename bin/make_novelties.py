#!/usr/bin/env python3
"""
Produce per-ingroup-species novelty candidate files.

For each ingroup species (identified by its Short ID in the config CSV),
writes novelties.<SHORT>.tsv containing proteins that are:
  1. Present in that species (presence matrix column == 1)
  2. Present in >= ingroup_min of ingroup proteomes
  3. Absent from all outgroup proteomes (all outgroup columns == 0)
  4. Absent from all outgroup genomes (no significant TBLASTN hit in tblastn_summary)

The annotated presence matrix (presence_matrix.function.tsv) must already
contain the Best_Swissprot and Pfam_Names columns added by
annotate_presence_matrix.py.  If tblastn_summary.tsv is absent, step 4 is
skipped with a warning.

When --cluster_tsv (mmseqs easy-cluster *_cluster.tsv) is supplied, each row
also gets family_id/family_size/family_members columns identifying the
mmseqs cluster its candidate belongs to — the same gene family recurring
across other ingroup species shows up under one family_id rather than as
unrelated rows in each species' table.
"""
import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from clusters import build_families, read_cluster_tsv
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES


def load_groups(config_csv):
    """Return (ingroup_ids, outgroup_ids) as lists of Short IDs.

    'IN'/'DISCOVERY_TARGET' are treated as the ingroup band, 'OUT'/'DISCOVERY_OUT'/
    'NEAR_INGROUP'/'BROAD_OUTGROUP' as the outgroup band (see lib/config_parser.py
    INGROUP_ROLES/OUTGROUP_ROLES) — this is the same coarse banding the report payload
    builders use, so novelty_discovery configs (DISCOVERY_TARGET/DISCOVERY_OUT) produce
    novelties.<SHORT>.tsv just like classic IN/OUT configs.
    """
    ingroup, outgroup = [], []
    with open(config_csv, newline='') as fh:
        for row in csv.DictReader(fh):
            short = row['Short'].strip()
            group = row['GROUP'].strip().upper()
            if group in INGROUP_ROLES:
                ingroup.append(short)
            elif group in OUTGROUP_ROLES:
                outgroup.append(short)
    return ingroup, outgroup


def load_tblastn_summary(path):
    """Return {protein_id: set_of_outgroup_ids_with_hits} from tblastn_summary.tsv."""
    tblastn: dict[str, set] = {}
    if not path or not os.path.exists(path):
        return tblastn
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        genome_cols = [c for c in (reader.fieldnames or []) if c != 'protein_id']
        for row in reader:
            hits = {g for g in genome_cols if row.get(g, '0') == '1'}
            if hits:
                tblastn[row['protein_id']] = hits
    return tblastn


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix', required=True,
                    help='presence_matrix.function.tsv (annotated, with Best_Swissprot/Pfam_Names)')
    ap.add_argument('--tblastn_summary',
                    help='tblastn_summary.tsv from summarize_tblastn.py (optional)')
    ap.add_argument('--cluster_tsv',
                    help='mmseqs easy-cluster *_cluster.tsv (rep -> member) for gene-family '
                         'grouping across ingroup species (optional)')
    ap.add_argument('--config', required=True,
                    help='Analysis config CSV (GROUP, Short, Species columns required)')
    ap.add_argument('--ingroup_min', type=float, default=0.75,
                    help='Min fraction of ingroup proteomes that must contain a hit (default: 0.75)')
    ap.add_argument('--output_dir', default='.',
                    help='Directory for novelties.<SHORT>.tsv output files (default: .)')
    ap.add_argument('--skip_tblastn_filter', action='store_true',
                    help='Do not exclude proteins with TBLASTN outgroup hits (hits still reported in tblastn_outgroup_hits column)')
    args = ap.parse_args()

    if not (0.0 < args.ingroup_min <= 1.0):
        sys.exit('--ingroup_min must be in (0, 1]')

    ingroup_ids, outgroup_ids = load_groups(args.config)
    tblastn = load_tblastn_summary(args.tblastn_summary)

    families = {}
    if args.cluster_tsv and os.path.exists(args.cluster_tsv):
        families = build_families(read_cluster_tsv(args.cluster_tsv))
    member_to_rep = {m: rep for rep, members in families.items() for m in members}

    if args.tblastn_summary and not tblastn and os.path.exists(args.tblastn_summary):
        print('WARNING: tblastn_summary loaded but contains no hits', file=sys.stderr)
    elif not args.tblastn_summary:
        print('WARNING: --tblastn_summary not provided; TBLASTN filter skipped', file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)

    # Load entire annotated matrix into memory (presence matrices are not huge)
    with open(args.matrix, newline='') as fin:
        reader = csv.DictReader(fin, delimiter='\t')
        header = reader.fieldnames or []

        # Validate expected columns
        missing = [c for c in ingroup_ids + outgroup_ids if c not in header]
        if missing:
            sys.exit(f'ERROR: columns not found in matrix: {missing}\nMatrix columns: {header}')

        rows = list(reader)

    # Desired output columns: protein_id, source_proteome (if present), per-species presence
    # columns, then annotation columns
    presence_cols = ingroup_ids + outgroup_ids
    annotation_cols = [c for c in header if c not in presence_cols]
    # Reorder: annotation first, then ingroup presence, then outgroup presence
    out_fields = (annotation_cols + ['family_id', 'family_size', 'family_members']
                  + ingroup_ids + outgroup_ids + ['tblastn_outgroup_hits'])

    counts: dict[str, tuple[int, int]] = {}

    for short in ingroup_ids:
        n_total = n_kept = 0
        out_path = os.path.join(args.output_dir, f'novelties.{short}.tsv')

        with open(out_path, 'w', newline='') as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter='\t',
                                    extrasaction='ignore', lineterminator='\n')
            writer.writeheader()

            for row in rows:
                n_total += 1

                # Must originate from the focal species
                if row.get('source_proteome', '') != short:
                    continue

                # Must meet ingroup_min fraction across all ingroup species
                ingroup_present = sum(int(row.get(c, '0')) for c in ingroup_ids)
                if ingroup_present / len(ingroup_ids) < args.ingroup_min:
                    continue

                # Must be absent from all outgroup proteomes
                if any(row.get(c, '0') == '1' for c in outgroup_ids):
                    continue

                # Must have no TBLASTN hits in any outgroup genome (unless filter is skipped)
                pid = row.get('protein_id', '')
                tblastn_hit_genomes = tblastn.get(pid, set())
                if tblastn_hit_genomes and not args.skip_tblastn_filter:
                    continue

                row['tblastn_outgroup_hits'] = ','.join(sorted(tblastn_hit_genomes))

                rep = member_to_rep.get(pid)
                if rep:
                    members = families[rep]
                    row['family_id'] = rep
                    row['family_size'] = len(members)
                    row['family_members'] = ','.join(m for m in members if m != pid)
                else:
                    row['family_id'] = ''
                    row['family_size'] = 1
                    row['family_members'] = ''

                writer.writerow(row)
                n_kept += 1

        counts[short] = (n_kept, n_total)
        print(f'{short}: wrote {n_kept} novelty candidates → {out_path}', file=sys.stderr)

    print(
        f'\nSummary: {sum(v[0] for v in counts.values())} total novelties across '
        f'{len(ingroup_ids)} ingroup species',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
