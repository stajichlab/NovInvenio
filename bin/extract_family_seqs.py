#!/usr/bin/env python3
"""Split an mmseqs cluster into per-family member FASTAs (families ≥ --min-members).

Reads an mmseqs `*_cluster.tsv` (two columns: representative_id<TAB>member_id) and
the concatenated FASTA the clustering was run on, and writes one FASTA per family
that has at least --min-members members into --outdir. Singletons (and any family
below the threshold) are skipped — see ADR-0002 Q6: a family profile needs ≥2
observed members, and a true species-unique orphan cannot be a clade novelty anyway.

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
    ap.add_argument('--outdir', required=True, help='Directory for per-family FASTAs')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.fasta)
    clusters = load_clusters(args.cluster_tsv)

    index_rows = []
    fam_idx = 0
    for rep, members in clusters.items():
        if len(members) < args.min_members:
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

    print(f"Wrote {len(index_rows)} families (≥{args.min_members} members) to {outdir}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
