#!/usr/bin/env python3
"""
Summarise per-outgroup TBLASTN results into a protein × genome presence matrix.

Each input TSV file should be named <SHORT>.tblastn.tsv where SHORT is the
outgroup species identifier matching the config CSV.  The cluster TSV
(mmseqs easy-cluster output) is used to expand representative-level hits to
all cluster members so the output can be joined to the per-protein presence
matrix.

Output TSV columns: protein_id, <SHORT1>, <SHORT2>, ...
  Value 1 = significant TBLASTN hit found for that protein in that genome.
  Value 0 = no hit.
"""
import argparse
import os
import sys
from collections import defaultdict


def parse_cluster_tsv(path):
    """Return {member_id: rep_id} from mmseqs easy-cluster TSV."""
    member_to_rep = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
    return member_to_rep


def parse_tblastn_hits(tsv_path, evalue_cutoff):
    """Return set of rep IDs with at least one hit below evalue_cutoff."""
    reps_with_hits = set()
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            qseqid = parts[0]
            try:
                evalue = float(parts[2])
            except ValueError:
                continue
            if evalue <= evalue_cutoff:
                reps_with_hits.add(qseqid)
    return reps_with_hits


def genome_short(tsv_path):
    """Extract <SHORT> from a filename like <SHORT>.tblastn.tsv."""
    basename = os.path.basename(tsv_path)
    # Strip .tblastn.tsv or .tblastn.tsv.gz suffix
    for suffix in ('.tblastn.tsv.gz', '.tblastn.tsv'):
        if basename.endswith(suffix):
            return basename[: -len(suffix)]
    # Fallback: strip everything from the first dot
    return basename.split('.')[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hits', nargs='+', required=True,
                    help='TBLASTN TSV files (<SHORT>.tblastn.tsv), one per outgroup genome')
    ap.add_argument('--cluster_tsv', required=True,
                    help='mmseqs easy-cluster *_cluster.tsv mapping rep → member')
    ap.add_argument('--evalue', type=float, default=1e-5,
                    help='E-value cutoff for counting a TBLASTN hit (default: 1e-5)')
    ap.add_argument('--output', required=True,
                    help='Output TSV: protein_id × outgroup genome hit matrix')
    args = ap.parse_args()

    member_to_rep = parse_cluster_tsv(args.cluster_tsv)
    # All proteins appearing as members (covers reps too, since reps are their own members)
    all_proteins = sorted(set(member_to_rep.keys()))

    # Build reverse map: rep → all its members
    rep_to_members: dict[str, list[str]] = defaultdict(list)
    for member, rep in member_to_rep.items():
        rep_to_members[rep].append(member)

    # Parse each TBLASTN file
    genome_ids = []
    reps_hit_per_genome: dict[str, set] = {}
    for tsv in args.hits:
        gid = genome_short(tsv)
        genome_ids.append(gid)
        reps_hit_per_genome[gid] = parse_tblastn_hits(tsv, args.evalue)

    if not genome_ids:
        sys.exit('ERROR: no TBLASTN hit files provided')

    # Expand rep hits to all cluster members
    # protein_hits[protein_id][genome_id] = 1 or 0
    protein_hit: dict[str, dict[str, int]] = {p: {g: 0 for g in genome_ids} for p in all_proteins}

    for gid, hit_reps in reps_hit_per_genome.items():
        for rep in hit_reps:
            # Mark rep itself and all its cluster members
            targets = rep_to_members.get(rep, [rep])
            for member in targets:
                if member in protein_hit:
                    protein_hit[member][gid] = 1
                # If member is not in the cluster map (shouldn't happen), still record
                else:
                    protein_hit.setdefault(member, {g: 0 for g in genome_ids})[gid] = 1

    with open(args.output, 'w') as fout:
        fout.write('protein_id\t' + '\t'.join(genome_ids) + '\n')
        for protein in sorted(protein_hit.keys()):
            vals = '\t'.join(str(protein_hit[protein][g]) for g in genome_ids)
            fout.write(f'{protein}\t{vals}\n')

    n_hit = sum(
        1 for p in protein_hit
        if any(protein_hit[p][g] for g in genome_ids)
    )
    print(
        f'Wrote {len(protein_hit)} proteins × {len(genome_ids)} genomes; '
        f'{n_hit} proteins have at least one outgroup TBLASTN hit.',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
