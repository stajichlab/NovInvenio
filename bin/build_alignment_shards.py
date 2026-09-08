#!/usr/bin/env python3
"""
Build per-genome, gzip-compressed JSON "alignment shards" for the interactive
report's TBLASTN alignment popup (a docs/-only feature -- see CLAUDE.md's
report constraints for why the offline results/ copy never gets this).

Only actual novelty/loss candidates are archived (bin/build_presence_matrix.py's
candidates.txt / loss_candidates.txt), not the full all-vs-all TBLASTN search
space. TBLASTN only ever queries cluster representative proteins
(modules/tblastn.nf's reps_fa input) -- a non-representative candidate's
alignment is really its representative's, recorded via an `aligned_as` field
so the report never silently misattributes it.

Sharding is one gzip JSON file per query genome (<genome_short>.json.gz), not
one file per candidate, to keep the published file count in the dozens (one
per genome in the config) instead of hundreds/thousands. A manifest.json
alongside the shards records a schema version and generation timestamp so the
client can detect a stale/mismatched shard before trusting it.

Shard schema: {protein_id: [hit, ...]} -- a list, never a single object, since
one candidate can have multiple HSPs/scaffolds against the same genome. Each
hit: {genome, sseqid, evalue, bitscore, pident, length, qstart, qend, sstart,
send, sframe, qseq, sseq[, aligned_as]}.

Usage:
  build_alignment_shards.py \
      --hits results/<project>/tblastn/*.tblastn.tsv \
      --candidates results/<project>/candidates.txt \
      --cluster_tsv results/<project>/clusters/clusters_cluster.tsv \
      --evalue 1e-5 \
      --project pezizo5 \
      --outdir results/<project>/alignments
"""
import argparse
import gzip
import json
import os
import sys
import time
from collections import defaultdict

SCHEMA_VERSION = 1

# Column order emitted by modules/tblastn.nf's -outfmt 6 string (see issue #70):
# qseqid sseqid evalue bitscore pident length qstart qend sstart send sframe qseq sseq
_COLS = (
    'qseqid', 'sseqid', 'evalue', 'bitscore', 'pident', 'length',
    'qstart', 'qend', 'sstart', 'send', 'sframe', 'qseq', 'sseq',
)


def parse_cluster_tsv(path):
    """{member_id: rep_id} from an mmseqs easy-cluster TSV (rep\tmember per
    line, including a rep\trep self-line)."""
    member_to_rep = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                member_to_rep[parts[1]] = parts[0]
    return member_to_rep


def load_candidate_ids(candidates_txt):
    """candidates.txt/loss_candidates.txt lines are '<Short>::<protein_id>' --
    only the bare protein_id is what the TBLASTN tsvs and cluster_tsv key on."""
    ids = set()
    with open(candidates_txt) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.add(line.split('::', 1)[-1])
    return ids


def genome_short(tsv_path):
    """Extract <SHORT> from a filename like <SHORT>.tblastn.tsv[.gz]."""
    basename = os.path.basename(tsv_path)
    for suffix in ('.tblastn.tsv.gz', '.tblastn.tsv'):
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
    return basename.split('.')[0]


def _parse_row(parts):
    if len(parts) < len(_COLS):
        return None
    return dict(zip(_COLS, parts))


def build_genome_shard(tsv_path, genome_id, evalue_cutoff, member_to_rep, candidate_ids):
    """{protein_id: [hit, ...]} for every candidate reachable (directly, or
    via its cluster representative) from this genome's TBLASTN hits."""
    rep_hits = defaultdict(list)
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            row = _parse_row(line.split('\t'))
            if row is None:
                continue
            try:
                evalue = float(row['evalue'])
            except ValueError:
                continue
            if evalue > evalue_cutoff:
                continue
            rep_hits[row['qseqid']].append(row)

    rep_to_members = defaultdict(list)
    for member, rep in member_to_rep.items():
        rep_to_members[rep].append(member)

    protein_hits = defaultdict(list)
    for rep_id, hits in rep_hits.items():
        members = rep_to_members.get(rep_id, [rep_id])
        for member in members:
            if member not in candidate_ids:
                continue
            for row in hits:
                entry = {
                    'genome': genome_id,
                    'sseqid': row['sseqid'],
                    'evalue': float(row['evalue']),
                    'bitscore': float(row['bitscore']),
                    'pident': float(row['pident']),
                    'length': int(row['length']),
                    'qstart': int(row['qstart']), 'qend': int(row['qend']),
                    'sstart': int(row['sstart']), 'send': int(row['send']),
                    'sframe': int(row['sframe']),
                    'qseq': row['qseq'], 'sseq': row['sseq'],
                }
                if member != rep_id:
                    entry['aligned_as'] = rep_id
                protein_hits[member].append(entry)
    return dict(protein_hits)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hits', nargs='+', required=True,
                    help='TBLASTN TSV files (<SHORT>.tblastn.tsv), one per genome')
    ap.add_argument('--candidates', required=True,
                    help='candidates.txt or loss_candidates.txt')
    ap.add_argument('--cluster_tsv', required=True,
                    help='mmseqs easy-cluster *_cluster.tsv mapping rep -> member')
    ap.add_argument('--evalue', type=float, default=1e-5,
                    help='E-value cutoff for archiving a hit (default: 1e-5)')
    ap.add_argument('--project', required=True,
                    help='Project name, recorded in manifest.json')
    ap.add_argument('--outdir', required=True,
                    help='Output directory for <genome>.json.gz shards + manifest.json')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    candidate_ids = load_candidate_ids(args.candidates)
    member_to_rep = parse_cluster_tsv(args.cluster_tsv)

    shard_files = []
    total_hits = 0
    for tsv in sorted(args.hits):
        gid = genome_short(tsv)
        protein_hits = build_genome_shard(tsv, gid, args.evalue, member_to_rep, candidate_ids)
        if not protein_hits:
            continue
        total_hits += sum(len(hits) for hits in protein_hits.values())

        shard_name = f'{gid}.json.gz'
        payload = json.dumps(protein_hits, separators=(',', ':')).encode('utf-8')
        # mtime=0 keeps the gzip header byte-identical across reruns with the
        # same content, so a rebuilt shard with unchanged data doesn't look
        # like a new file to anything hashing/diffing the published tree.
        with open(os.path.join(args.outdir, shard_name), 'wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as fh:
                fh.write(payload)
        shard_files.append(shard_name)

    manifest = {
        'project': args.project,
        'schema_version': SCHEMA_VERSION,
        'generated_at': int(time.time()),
        'evalue_cutoff': args.evalue,
        'shards': sorted(shard_files),
    }
    with open(os.path.join(args.outdir, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2)
        fh.write('\n')

    print(
        f'Wrote {len(shard_files)} genome shard(s), {total_hits} candidate hit(s), to {args.outdir}',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
