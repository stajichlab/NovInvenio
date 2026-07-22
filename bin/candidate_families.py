#!/usr/bin/env python3
"""Family-as-cluster: derive VALIDATE/ANNOTATE/REPORT cluster inputs from the profile
pathway's existing gene families, skipping the redundant candidate re-clustering (ADR-0002
Q7, issue #11).

When --cluster_tool mmseqs, MMSEQS_FAMILY_CLUSTER has already clustered the whole ingroup
into gene families. Re-running mmseqs on just the candidates (the generic CLUSTER workflow)
is redundant and re-defines the family unit from a fresh clustering. Instead, this reuses
the profile HMM's own families, restricted to the ones that contain a candidate, and emits
the same three artifacts the CLUSTER workflow does:

  * cluster_tsv     : the candidate families' rows of families_cluster.tsv (rep<TAB>member) —
                      TBLASTN member-expansion + REPORT family grouping.
  * representatives : the family reps of candidate-containing families (a subset of
                      families_rep_seq.fasta) — TBLASTN probes.
  * candidates_fa   : the candidate protein sequences (from the seed-group FASTA) — ANNOTATE
                      input.

A "candidate family" is any family whose representative or any member appears in
candidates.txt. Candidates that fall outside every profiled family (should not happen — a
candidate is a profiled-family member by construction — but possible for a singleton-dropped
edge case) are still written to candidates_fa so annotation is not lost; they simply carry
no family (matching the pairwise pathway's singleton behaviour).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import read_fasta, write_fasta  # noqa: E402


def load_candidate_ids(candidates_txt):
    """protein ids from candidates.txt lines of `source_proteome::protein_id`."""
    ids = set()
    with open(candidates_txt) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.add(line.split('::', 1)[1] if '::' in line else line)
    return ids


def load_family_membership(cluster_tsv):
    """Return (member_to_rep, rep_to_members) from families_cluster.tsv (rep<TAB>member)."""
    member_to_rep = {}
    rep_to_members = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            member_to_rep[member] = rep
            rep_to_members.setdefault(rep, []).append(member)
    return member_to_rep, rep_to_members


def candidate_family_reps(candidate_ids, member_to_rep):
    """Reps of families that contain at least one candidate."""
    reps = set()
    for cid in candidate_ids:
        rep = member_to_rep.get(cid)
        if rep is not None:
            reps.add(rep)
    return reps


def write_cluster_tsv(reps, rep_to_members, out_path):
    with open(out_path, 'w') as fh:
        for rep in sorted(reps):
            for member in rep_to_members.get(rep, []):
                fh.write(f'{rep}\t{member}\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--candidates', required=True, help='candidates.txt (source::protein)')
    ap.add_argument('--family-cluster-tsv', required=True, dest='family_cluster_tsv',
                    help='families_cluster.tsv (rep<TAB>member) from MMSEQS_FAMILY_CLUSTER')
    ap.add_argument('--family-reps', required=True, dest='family_reps',
                    help='families_rep_seq.fasta (one rep sequence per family)')
    ap.add_argument('--seed-fasta', required=True, dest='seed_fasta',
                    help='concatenated seed-group proteome FASTA (candidate sequences); the '
                         'ingroup for novelty, the outgroup for the loss mirror')
    ap.add_argument('--out-cluster-tsv', required=True, dest='out_cluster_tsv')
    ap.add_argument('--out-representatives', required=True, dest='out_representatives')
    ap.add_argument('--out-candidates-fa', required=True, dest='out_candidates_fa')
    args = ap.parse_args()

    candidate_ids = load_candidate_ids(args.candidates)
    member_to_rep, rep_to_members = load_family_membership(args.family_cluster_tsv)
    reps = candidate_family_reps(candidate_ids, member_to_rep)

    # cluster_tsv: candidate families only.
    write_cluster_tsv(reps, rep_to_members, args.out_cluster_tsv)

    # representatives: the reps of candidate families (subset of families_rep_seq.fasta).
    rep_records = read_fasta(args.family_reps)
    write_fasta([rep_records[r] for r in sorted(reps) if r in rep_records],
                args.out_representatives)

    # candidates_fa: the candidate protein sequences (from the seed-group proteomes).
    seed_records = read_fasta(args.seed_fasta)
    cand_records = [seed_records[c] for c in sorted(candidate_ids) if c in seed_records]
    write_fasta(cand_records, args.out_candidates_fa)

    missing = len(candidate_ids) - len(cand_records)
    print(f'candidate families: {len(reps)} reps over {len(candidate_ids)} candidates; '
          f'{len(cand_records)} sequences written'
          + (f' ({missing} candidate ids not found in seed FASTA)' if missing else ''),
          file=sys.stderr)


if __name__ == '__main__':
    main()
