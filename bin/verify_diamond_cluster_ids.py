#!/usr/bin/env python3
"""Verify diamond cluster's *_cluster.tsv keeps the full FASTA header token.

Unlike mmseqs easy-cluster (see bin/restore_mmseqs_cluster_ids.py and
lib/fasta.py's mmseqs_id() docstring), `diamond cluster` is documented (and
empirically observed against raw UniProt-style headers -- see
.living/learnings.md's 2026-09-08 entry) to keep the whole first-token
header verbatim as its sequence id, rather than collapsing sp|ACC|NAME-style
deflines to a bare accession. That means diamond needs no restoration step
the way mmseqs does -- but that expectation had never been checked against
this pipeline's own header convention (bin/pangenome_prefix_fasta.py's
Short<id_sep>orig_id prefix, including the untested case where orig_id
itself contains id_sep characters, e.g. a UniProt original header prefixed
to Short|sp|ACC|NAME). This script is the fail-loud safety net that checks
the expectation actually held for a given run, instead of assuming it (see
docs/adr/0003-diamond-tier1-clustering-backend.md). It never rewrites
*_cluster.tsv -- if every id checks out, the file is untouched.

Fails loudly on:
  - a cluster.tsv id matching no header in the input FASTA at all.
  - a cluster.tsv id matching a header only after removing an mmseqs-style
    sp|/tr|/gb|/... wrapper -- i.e. diamond appears to have collapsed a
    header the way mmseqs does, which is a distinct, more alarming failure
    than a plain unknown id and gets its own message.

Usage:
  verify_diamond_cluster_ids.py \
      --input-fasta <the FASTA diamond clustered> \
      --cluster-tsv <diamond cluster's *_cluster.tsv, read-only>
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import mmseqs_id  # noqa: E402


def read_full_ids(fasta_path: str) -> set[str]:
    """Full (first-whitespace-token) header ids from a FASTA file."""
    ids: set[str] = set()
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                ids.add(line[1:].split()[0])
    return ids


def build_collapsed_lookup(full_ids: set[str]) -> dict[str, str]:
    """{mmseqs-collapsed_id: full_id} for every full_id that mmseqs_id()
    would actually shorten -- used only to recognize an id that looks like
    an mmseqs-style collapse, not to accept it."""
    lookup: dict[str, str] = {}
    for full_id in full_ids:
        collapsed = mmseqs_id(full_id)
        if collapsed != full_id:
            lookup[collapsed] = full_id
    return lookup


def verify_cluster_tsv(cluster_tsv_path: str, full_ids: set[str]) -> None:
    """Fail loud (sys.exit) on the first id in cluster_tsv_path that isn't a
    verbatim member of full_ids. Never modifies cluster_tsv_path."""
    collapsed_lookup = build_collapsed_lookup(full_ids)
    with open(cluster_tsv_path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip('\n')
            if not line:
                continue
            for x in line.split('\t'):
                if x in full_ids:
                    continue
                if x in collapsed_lookup:
                    sys.exit(
                        f"ERROR: {cluster_tsv_path} line {lineno}: id '{x}' does not match any "
                        f"FASTA header verbatim, but matches '{collapsed_lookup[x]}' after "
                        "stripping an mmseqs-style sp|/tr|/... wrapper -- diamond appears to "
                        "have collapsed this header the way mmseqs does, which this pipeline's "
                        "diamond backend does not expect and has no restoration step for. Do "
                        "not treat this as a plain unknown id -- investigate the diamond "
                        "version/invocation before trusting this cluster.tsv."
                    )
                sys.exit(
                    f"ERROR: {cluster_tsv_path} line {lineno}: id '{x}' does not match any "
                    f"header in the clustered FASTA -- diamond cluster output cannot be "
                    "trusted for this run"
                )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-fasta', required=True, dest='input_fasta',
                     help='The exact FASTA fed to diamond cluster')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                     help="diamond cluster's *_cluster.tsv, verified read-only")
    args = ap.parse_args()

    full_ids = read_full_ids(args.input_fasta)
    verify_cluster_tsv(args.cluster_tsv, full_ids)

    print(f"Verified {args.cluster_tsv}: every id matches a full FASTA header", file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
