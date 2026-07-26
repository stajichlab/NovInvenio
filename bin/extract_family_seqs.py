#!/usr/bin/env python3
"""Split an mmseqs cluster into per-family member FASTAs (--min-members <= size <=
--max-members).

Reads an mmseqs `*_cluster.tsv` (two columns: representative_id<TAB>member_id) and
the concatenated FASTA the clustering was run on, and writes one FASTA per family
whose size is in [--min-members, --max-members] into --outdir. Singletons (and any
family below the threshold) are skipped — see ADR-0002 Q6: a family profile needs ≥2
observed members, and a true species-unique orphan cannot be a clade novelty anyway.

Families ABOVE --max-members are also skipped (issue #22 follow-up): a family this
large is typically an ancient, broadly-conserved multi-copy superfamily (e.g. an
aldehyde-dehydrogenase-like domain) that is present across most of the sampled tree
by construction — it cannot pass the novelty predicate (absent from every outgroup)
or the loss predicate (absent from most of the ingroup) either way, so profiling it
buys no signal. It IS expensive: famsa's alignment cost grows superlinearly with
member count, and a handful of 200+-member families landing in the same BUILD_CHUNK
task risked exceeding the SLURM time limit. Skipped families are logged to
--oversized-report for visibility, not silently dropped.

Output filenames are `fam_<index>.faa`; a `families.tsv` index
(family_index<TAB>representative_id<TAB>n_members) is written alongside so the
representative id (used as the HMM name downstream) is recoverable.
"""
import argparse
import sys
from collections import OrderedDict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import read_fasta  # noqa: E402


def load_clusters(cluster_tsv):
    """Return OrderedDict rep_id -> [member_id, ...] preserving first-seen order."""
    clusters: OrderedDict[str, list] = OrderedDict()
    with open(cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            rep, member = line.split('\t')[:2]
            clusters.setdefault(rep, []).append(member)
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                    help='mmseqs *_cluster.tsv (rep_id<TAB>member_id)')
    ap.add_argument('--fasta', required=True,
                    help='Concatenated FASTA the clustering was run on')
    ap.add_argument('--min-members', type=int, default=2, dest='min_members',
                    help='Minimum members for a family to be emitted (default 2)')
    ap.add_argument('--max-members', type=int, default=None, dest='max_members',
                    help='Maximum members for a family to be emitted (default: no limit). '
                         'Families above this are skipped as almost-certainly-uninteresting '
                         'multi-copy superfamilies and reported to --oversized-report')
    ap.add_argument('--outdir', required=True, help='Directory for per-family FASTAs')
    ap.add_argument('--oversized-report', default=None, dest='oversized_report',
                    help='Optional TSV (representative_id, n_members) listing families '
                         'skipped for exceeding --max-members')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.fasta)
    clusters = load_clusters(args.cluster_tsv)

    index_rows = []
    oversized_rows = []
    fam_idx = 0
    for rep, members in clusters.items():
        if len(members) < args.min_members:
            continue
        if args.max_members is not None and len(members) > args.max_members:
            oversized_rows.append((rep, len(members)))
            continue
        fam_idx += 1
        fam_path = outdir / f'fam_{fam_idx:06d}.faa'
        with open(fam_path, 'w') as out:
            for m in members:
                if m not in records:
                    # A member id absent from the FASTA means a mismatch between the
                    # cluster tsv and the sequences it was built from — fail loud.
                    sys.exit(f"ERROR: cluster member '{m}' not found in {args.fasta}")
                out.write(f'>{m}\n{records[m].seq}\n')
        index_rows.append((fam_idx, rep, len(members)))

    with open(outdir / 'families.tsv', 'w') as fh:
        fh.write('family_index\trepresentative_id\tn_members\n')
        for idx, rep, n in index_rows:
            fh.write(f'fam_{idx:06d}\t{rep}\t{n}\n')

    if args.oversized_report:
        with open(args.oversized_report, 'w') as fh:
            fh.write('representative_id\tn_members\n')
            for rep, n in oversized_rows:
                fh.write(f'{rep}\t{n}\n')

    max_note = f", {len(oversized_rows)} oversized (>{args.max_members}) skipped" \
        if args.max_members is not None else ""
    print(f"Wrote {len(index_rows)} families (≥{args.min_members} members) to {outdir}"
          f"{max_note}", file=sys.stderr)


if __name__ == '__main__':
    main()
