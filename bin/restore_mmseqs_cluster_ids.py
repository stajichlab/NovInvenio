#!/usr/bin/env python3
"""Restore mmseqs easy-cluster's collapsed FASTA header IDs in *_cluster.tsv.

mmseqs easy-cluster recognizes several FASTA defline conventions (UniProt's
sp|ACC|NAME/tr|ACC|NAME, and several NCBI-style ones -- see lib/fasta.py's
mmseqs_id() docstring) and reports a collapsed field as the sequence ID in
its own *_cluster.tsv output, while its *_rep_seq.fasta/*_all_seqs.fasta
outputs from the SAME run keep the full header token. Every other artifact
in this pipeline (candidates.txt, presence_matrix.tsv, novelties.<Short>.tsv,
TBLASTN's qseqid, the report payload) treats the full header token as
protein_id -- so any code joining a *_cluster.tsv row against protein_id
silently matches nothing for a UniProt-sourced study.

This script corrects *_cluster.tsv in place, immediately after mmseqs
produces it (wired into modules/mmseqs_cluster.nf and
modules/mmseqs_family_cluster.nf), using the exact FASTA fed to mmseqs to
reconstruct {collapsed_id: full_header} and rewrite every cluster.tsv value
back to the full header -- before any of this pipeline's ~16-22 cluster.tsv
consumers ever reads the file. See
docs/superpowers/specs/2026-09-09-mmseqs-cluster-id-restoration-design.md.

Fails loudly (never silently mismatches) on:
  - two distinct FASTA headers collapsing to the same id (a genuine
    accession collision in the input, most likely when concatenating
    multiple independently-sourced proteomes for family clustering)
  - a cluster.tsv value matching neither a known full header nor its
    computed collapsed form (an mmseqs output this script cannot reconcile)

Usage:
  restore_mmseqs_cluster_ids.py \
      --input-fasta <the FASTA mmseqs clustered> \
      --cluster-tsv <mmseqs easy-cluster's *_cluster.tsv, corrected in place>
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from fasta import mmseqs_id  # noqa: E402


def build_collapsed_to_full(fasta_path: str) -> dict[str, str]:
    """{collapsed_id: full_id} for every record in fasta_path, using the same
    header the FASTA itself carries (first whitespace-delimited token) as
    full_id. Fails loud on a collapsed-id collision between two distinct
    full ids."""
    collapsed_to_full: dict[str, str] = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            full_id = line[1:].split()[0]
            collapsed_id = mmseqs_id(full_id)
            if collapsed_id in collapsed_to_full and collapsed_to_full[collapsed_id] != full_id:
                sys.exit(
                    f"ERROR: {fasta_path}: both '{collapsed_to_full[collapsed_id]}' and "
                    f"'{full_id}' collapse to the same mmseqs id '{collapsed_id}' -- "
                    f"cannot unambiguously restore cluster.tsv"
                )
            collapsed_to_full[collapsed_id] = full_id
    return collapsed_to_full


def restore_cluster_tsv(cluster_tsv_path: str, collapsed_to_full: dict[str, str]) -> str:
    """Return the corrected cluster.tsv content (rep\tmember per line)."""
    full_ids = set(collapsed_to_full.values())
    out_lines = []
    with open(cluster_tsv_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            fixed = []
            for x in parts:
                if x in collapsed_to_full:
                    fixed.append(collapsed_to_full[x])
                elif x in full_ids:
                    fixed.append(x)  # already full -- idempotent no-op
                else:
                    sys.exit(
                        f"ERROR: {cluster_tsv_path}: id '{x}' matches neither a known "
                        f"full header nor a computed mmseqs id -- cannot restore"
                    )
            out_lines.append('\t'.join(fixed))
    return '\n'.join(out_lines) + ('\n' if out_lines else '')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-fasta', required=True, dest='input_fasta',
                     help='The exact FASTA fed to mmseqs easy-cluster')
    ap.add_argument('--cluster-tsv', required=True, dest='cluster_tsv',
                     help='mmseqs easy-cluster *_cluster.tsv, corrected in place')
    args = ap.parse_args()

    collapsed_to_full = build_collapsed_to_full(args.input_fasta)
    corrected = restore_cluster_tsv(args.cluster_tsv, collapsed_to_full)
    with open(args.cluster_tsv, 'w') as fh:
        fh.write(corrected)

    print(f"Restored {len(collapsed_to_full)} header(s) in {args.cluster_tsv}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
