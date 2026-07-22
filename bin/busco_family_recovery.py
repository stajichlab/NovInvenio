#!/usr/bin/env python3
"""BUSCO single-copy family-recovery metric for the family-profile sweep (ADR-0002 Q8).

The primary *curation-free* family-quality metric (issue #6): a single-copy BUSCO
ortholog is, by definition, one gene present once per species, so a well-behaved gene
clustering should recover the complete copies of each such BUSCO as exactly **one gene
family** — not split across families (under-clustering) and not merged with an unrelated
BUSCO (over-clustering). This measures both failure modes directly from the mmseqs cluster
membership, with no hand curation.

Inputs:
  * BUSCO `full_table.tsv` files, one per ingroup proteome (from `busco -m proteins -l
    fungi_odb12`, or any BUSCO run). Columns (tab, `#`-comment header): Busco id, Status,
    Sequence(=protein_id), ... . Only Status == Complete rows are used (single-copy;
    Duplicated/Fragmented/Missing are ignored). Pass as `SHORT=path` so each table's
    proteins are attributed to a proteome. Alternatively pass a pre-consolidated
    `--busco-map` (busco_id<TAB>species<TAB>protein_id[<TAB>status]).
  * The run's mmseqs `*_cluster.tsv` (rep<TAB>member) and `families.tsv` (profiled,
    >=min-member families) — the same pair the presence matrix is built from.

A BUSCO is *eligible* when it is Complete in >= --min-species proteomes (an ortholog you
would expect the clustering to group). For each eligible BUSCO its complete proteins are
mapped to family representatives:
  * all in one family, none dropped        -> recovered (single family)
  * spread over >1 family                  -> split  (under-clustering)
  * some members not in any profiled family -> partial (dropped as singletons)
recovery_rate = recovered / eligible. Over-clustering is reported separately as the number
of families that contain proteins from >1 distinct BUSCO.
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

BUSCO_COLS = ('busco_id', 'status', 'sequence')


def parse_busco_full_table(path, species):
    """Yield (busco_id, species, protein_id) for Complete rows of a BUSCO full_table.tsv."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            busco_id, status, sequence = f[0], f[1], f[2]
            if status != 'Complete':
                continue
            # BUSCO may append :start-end to the sequence id; strip to the bare protein id.
            protein_id = sequence.split(':')[0]
            if protein_id:
                yield busco_id, species, protein_id


def load_busco_map(map_path=None, tables=None):
    """busco_id -> {species: protein_id} for Complete single-copy BUSCOs.

    From a consolidated --busco-map (busco_id<TAB>species<TAB>protein_id[...]) or from raw
    per-species full_table.tsv files given as SHORT=path.
    """
    mapping: dict[str, dict[str, str]] = defaultdict(dict)
    if map_path:
        with open(map_path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 3 or f[0].startswith('#'):
                    continue
                # a status column, if present, must be Complete
                if len(f) >= 4 and f[3] and f[3] != 'Complete':
                    continue
                mapping[f[0]][f[1]] = f[2]
    for spec in tables or []:
        if '=' not in spec:
            sys.exit(f'--tables entries must be SHORT=path (got: {spec!r})')
        short, path = spec.split('=', 1)
        for busco_id, species, protein_id in parse_busco_full_table(path, short):
            mapping[busco_id][species] = protein_id
    return mapping


def load_profiled_reps(families_tsv):
    reps = set()
    with open(families_tsv) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1]:
                reps.add(parts[1])
    return reps


def load_member_to_rep(cluster_tsv, keep_reps):
    member_to_rep = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep in keep_reps:
                member_to_rep[member] = rep
    return member_to_rep


def score_recovery(busco_map, member_to_rep, min_species):
    """Per-BUSCO recovery + over-clustering. Returns (per_busco rows, summary dict)."""
    rep_to_buscos: dict[str, set] = defaultdict(set)
    per_busco = []
    eligible = recovered = split = partial = 0

    for busco_id, by_species in sorted(busco_map.items()):
        if len(by_species) < min_species:
            continue
        eligible += 1
        reps = set()
        n_unmapped = 0
        for protein_id in by_species.values():
            rep = member_to_rep.get(protein_id)
            if rep is None:
                n_unmapped += 1
            else:
                reps.add(rep)
                rep_to_buscos[rep].add(busco_id)

        if len(reps) == 1 and n_unmapped == 0:
            verdict = 'recovered'
            recovered += 1
        elif len(reps) > 1:
            verdict = 'split'
            split += 1
        else:
            verdict = 'partial'
            partial += 1

        per_busco.append({
            'busco_id': busco_id,
            'n_species': len(by_species),
            'n_families': len(reps),
            'n_unmapped': n_unmapped,
            'verdict': verdict,
            'families': ';'.join(sorted(reps)),
        })

    overmerged = sum(1 for reps in rep_to_buscos.values() if len(reps) > 1)
    summary = {
        'eligible_buscos': eligible,
        'recovered_single_family': recovered,
        'split_across_families': split,
        'partial_unmapped': partial,
        'recovery_rate': (recovered / eligible) if eligible else None,
        'overmerged_families': overmerged,
    }
    return per_busco, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tables', nargs='*', default=[],
                    help='BUSCO full_table.tsv per proteome, as SHORT=path')
    ap.add_argument('--busco-map', default=None, dest='busco_map',
                    help='consolidated busco_id<TAB>species<TAB>protein_id[<TAB>status] TSV')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep<TAB>member)')
    ap.add_argument('--families', required=True, help='families.tsv (profiled families)')
    ap.add_argument('--min-species', type=int, default=2, dest='min_species',
                    help='BUSCO must be Complete in >= this many species to be eligible')
    ap.add_argument('--output', required=True, help='per-BUSCO results TSV')
    ap.add_argument('--summary', default=None,
                    help='summary TSV (default: *.summary.tsv next to --output)')
    args = ap.parse_args()

    if not args.tables and not args.busco_map:
        sys.exit('provide --tables SHORT=path ... and/or --busco-map')

    busco_map = load_busco_map(args.busco_map, args.tables)
    profiled_reps = load_profiled_reps(args.families)
    member_to_rep = load_member_to_rep(args.cluster_tsv, profiled_reps)
    per_busco, summary = score_recovery(busco_map, member_to_rep, args.min_species)

    fields = ['busco_id', 'n_species', 'n_families', 'n_unmapped', 'verdict', 'families']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(per_busco)

    summary_path = args.summary or (str(Path(args.output).with_suffix('')) + '.summary.tsv')
    with open(summary_path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['metric', 'value'])
        for k, v in summary.items():
            w.writerow([k, '' if v is None else v])

    rate = summary['recovery_rate']
    print(f"BUSCO recovery: "
          f"{'n/a' if rate is None else f'{rate:.3f}'} "
          f"({summary['recovered_single_family']}/{summary['eligible_buscos']} eligible); "
          f"{summary['split_across_families']} split, {summary['partial_unmapped']} partial, "
          f"{summary['overmerged_families']} over-merged families", file=sys.stderr)


if __name__ == '__main__':
    main()
