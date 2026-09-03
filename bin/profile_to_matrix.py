#!/usr/bin/env python3
"""Build the presence/absence matrix + candidates from family HMM search results.

The family-profile pathway (ADR-0002) replaces the pairwise search: the ingroup is
clustered into gene families, each family becomes an HMM, and every proteome (ingroup
AND outgroup — ADR-0002 Q1) is scanned with `hmmsearch --domtblout`. This script turns
those per-proteome domtblout files, plus the cluster membership, into the SAME artifacts
`build_presence_matrix.py` emits, so the rest of the pipeline is reused unchanged:

  * presence_matrix.tsv : columns `protein_id, source_proteome, <sorted proteome shorts>`
                          with 0/1 cells (one row per ingroup family member protein).
  * candidates.txt      : `source_proteome::protein_id`, one per line.

Presence is uniform and alignment-based (ADR-0002 Q1/Q5): a family is present in a
proteome iff some target sequence there has a hit for the family HMM with full-sequence
`E < --evalue` AND profile-coverage ≥ --min-coverage (merged HMM-coordinate span across
that target's domains, divided by the HMM length -- never mixing one target's E-value
with another's coverage, nor pooling coverage across different targets; see
lib/family_presence.py's parse_domtblout, shared with the novelty_discovery pathway).
Cluster membership only seeds the family; the member's own source proteome is always
marked present (it is a member by construction, mirroring build_presence_matrix.py's
"source proteome present by definition").

Family granularity (min members) is applied upstream in extract_family_seqs.py; here a
family is any cluster whose representative appears in --cluster-tsv AND was profiled
(i.e. is present in --families, the families.tsv index). The candidate keep-rule mirrors
build_presence_matrix.py exactly: query_frac >= --ingroup-min-frac AND other_frac <=
--other-max-frac, with --query-group selecting the direction (IN = novelty, OUT = loss).
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import INGROUP_ROLES, OUTGROUP_ROLES, parse_config  # noqa: E402
from family_presence import parse_domtblout  # noqa: E402


def load_protein_map(path):
    """protein_id -> proteome Short, from a two-column TSV (no header)."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            pid, short = line.split('\t')[:2]
            mapping[pid] = short
    return mapping


def load_families(families_tsv):
    """representative_id set for families that were actually profiled (≥ min members)."""
    reps = set()
    with open(families_tsv) as fh:
        header = fh.readline()  # family_index, representative_id, n_members
        del header
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                reps.add(parts[1])
    return reps


def load_cluster_members(cluster_tsv, keep_reps):
    """rep_id -> [member_id, ...] for reps in keep_reps (profiled families only)."""
    members = defaultdict(list)
    with open(cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            rep, member = line.split('\t')[:2]
            if rep in keep_reps:
                members[rep].append(member)
    return members


def short_from_domtblout_name(path):
    """`<Short>.family.domtblout` -> Short (strip the known suffix, else use the stem)."""
    name = Path(path).name
    for suffix in ('.family.domtblout', '.domtblout'):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--domtblout', nargs='+', required=True,
                    help='Per-proteome hmmsearch domtblout files, named <Short>.family.domtblout')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep_id<TAB>member_id)')
    ap.add_argument('--families', required=True,
                    help='families.tsv from extract_family_seqs.py (profiled families)')
    ap.add_argument('--protein-map', required=True, dest='protein_map',
                    help='protein_id<TAB>proteome_short TSV for the ingroup')
    ap.add_argument('--config', required=True, help='Analysis description CSV')
    ap.add_argument('--evalue', type=float, default=1e-3,
                    help='Full-sequence E-value ceiling for presence (default 1e-3)')
    ap.add_argument('--min-coverage', type=float, default=0.5, dest='min_coverage',
                    help='Minimum profile coverage for presence (default 0.5)')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75,
                    dest='ingroup_min_frac',
                    help="Presence threshold within --query-group")
    ap.add_argument('--query-group', choices=['IN', 'OUT'], default='IN',
                    dest='query_group',
                    help='IN = novelty direction (default); OUT = loss direction')
    ap.add_argument('--other-max-frac', type=float, default=0.0,
                    dest='other_max_frac',
                    help='Max fraction of the other group a candidate may be present in')
    ap.add_argument('--output-matrix', required=True)
    ap.add_argument('--output-candidates', required=True)
    args = ap.parse_args()

    samples = parse_config(args.config)
    # Coarse banding (see config_parser.INGROUP_ROLES/OUTGROUP_ROLES): NEAR_INGROUP/
    # BROAD_OUTGROUP/DISCOVERY_OUT fold into the outgroup band here, same as
    # make_novelties.py and build_presence_matrix.py, so those samples still get a
    # matrix column instead of being silently dropped.
    ingroup_ids = {s.short for s in samples if s.group in INGROUP_ROLES}
    outgroup_ids = {s.short for s in samples if s.group in OUTGROUP_ROLES}
    all_ids = ingroup_ids | outgroup_ids
    query_ids = ingroup_ids if args.query_group == 'IN' else outgroup_ids
    other_ids = outgroup_ids if args.query_group == 'IN' else ingroup_ids

    protein_to_proteome = load_protein_map(args.protein_map)
    profiled_reps = load_families(args.families)
    cluster_members = load_cluster_members(args.cluster_tsv, profiled_reps)

    # family rep -> set of proteome shorts where the family HMM is present.
    family_presence: dict[str, set] = defaultdict(set)
    for dom_path in args.domtblout:
        short = short_from_domtblout_name(dom_path)
        for rep in parse_domtblout(dom_path, args.evalue, args.min_coverage):
            family_presence[rep].add(short)

    sorted_ids = sorted(all_ids)
    columns = ['protein_id', 'source_proteome'] + sorted_ids

    rows = []
    for rep, members in cluster_members.items():
        present = set(family_presence.get(rep, set()))
        for member in members:
            source = protein_to_proteome.get(member)
            if source is None:
                # Member not in the ingroup map — clustering seeded from the ingroup,
                # so every member should map. Skip defensively rather than mislabel.
                continue
            row = {'protein_id': member, 'source_proteome': source}
            member_present = present | {source}  # source present by construction
            for sp in sorted_ids:
                row[sp] = int(sp in member_present)
            rows.append(row)

    matrix = pd.DataFrame(rows, columns=columns)
    matrix.to_csv(args.output_matrix, sep='\t', index=False)

    n_query = len(query_ids)
    query_cols = sorted(query_ids)
    other_cols = sorted(other_ids)
    if matrix.empty or n_query == 0:
        candidates = []
    else:
        query_count = matrix[query_cols].sum(axis=1)
        other_count = matrix[other_cols].sum(axis=1) if other_cols else 0
        other_frac = (other_count / len(other_cols)) if other_cols else 0
        keep = (query_count / n_query >= args.ingroup_min_frac) & (other_frac <= args.other_max_frac)
        kept = matrix[keep]
        candidates = (kept['source_proteome'] + '::' + kept['protein_id']).tolist()

    with open(args.output_candidates, 'w') as fh:
        if candidates:
            fh.write('\n'.join(candidates) + '\n')

    direction = 'ingroup' if args.query_group == 'IN' else 'outgroup'
    print(f"Candidates: {len(candidates)} / {len(matrix)} {direction} family-member "
          f"proteins pass thresholds", file=sys.stderr)


if __name__ == '__main__':
    main()
