#!/usr/bin/env python3
"""Two-tier, swappable (mmseqs2 | diamond) clustering backend for the
pangenome-profiling subworkflow. Tier 1 is the allele/ortholog unit used by
every downstream step (presence, frequency, co-occurrence, synteny); tier 2
re-clusters tier-1 representatives loosely, purely as a superfamily
annotation label -- it is never used for presence/frequency calls.

Ported from NovInvenio_Investigations' Afumigatus_pangenome study
(bin/cluster_backend.py). Per NEXTFLOW_MIGRATION_NOTES.md section B, the
study version hardcoded its identity/coverage thresholds and wrote a fixed
``tmp_mmseqs`` scratch directory to the current working directory with no
``--threads``/``-p`` plumbing; this version takes all three as CLI options
(each Nextflow process passes its own `task.cpus` and $SCRATCH/$TMPDIR) so a
different species/dataset can use different identity thresholds without
editing this script, and two concurrent invocations in the same work
directory (there should never be two in the SAME process, but a
$SCRATCH-per-task convention still argues for not hardcoding a literal
relative name) don't collide.

Usage:
  pangenome_cluster_backend.py mmseqs-tier1 --fasta all_strains.fa --out_prefix tier1 \\
      --min_seq_id 0.9 --cov 0.8 --threads 8
  pangenome_cluster_backend.py mmseqs-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2 \\
      --min_seq_id 0.4 --cov 0.8 --threads 8
  pangenome_cluster_backend.py diamond-tier1 --fasta all_strains.fa --out_prefix tier1 \\
      --approx_id 90 --member_cover 80 --threads 8
  pangenome_cluster_backend.py diamond-tier2 --fasta tier1_rep_seq.fasta --out_prefix tier2 \\
      --approx_id 40 --member_cover 80 --threads 8

Input FASTA requirement: every strain's (isoform-collapsed) proteome
concatenated, with each header Short-PREFIXED as
``><Short><id_sep><original_protein_id>`` (see pangenome_prefix_fasta.py).
The cluster TSV carries nothing but sequence IDs, so that prefix is the
only way pangenome_build_presence_matrix.py can recover which strain a
family member came from.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))
from pangenome_matrix import build_families, read_cluster_tsv  # noqa: E402


def _scratch_tmp_dir(prefix: str) -> str:
    """A tmp directory under $SCRATCH (node-local, if set), $TMPDIR, or /tmp,
    never a fixed cwd-relative name -- see this script's module docstring."""
    base = os.environ.get("SCRATCH") or os.environ.get("TMPDIR") or "/tmp"
    return tempfile.mkdtemp(prefix=f"{prefix}_", dir=base)


def run_mmseqs_cluster(
    fasta: Path, out_prefix: str, min_seq_id: float, cov: float,
    cluster_reassign: bool, threads: int,
) -> Path:
    tmp_dir = _scratch_tmp_dir("mmseqs")
    cmd = [
        "mmseqs", "easy-cluster", str(fasta), out_prefix, tmp_dir,
        "--min-seq-id", str(min_seq_id), "-c", str(cov), "--cov-mode", "0",
        "--threads", str(threads),
    ]
    if cluster_reassign:
        cmd.append("--cluster-reassign")
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def run_diamond_cluster(
    fasta: Path, out_prefix: str, approx_id: float, member_cover: float, threads: int,
) -> Path:
    cmd = [
        "diamond", "cluster", "-d", str(fasta), "-o", f"{out_prefix}_cluster.tsv",
        "--approx-id", str(approx_id), "--member-cover", str(member_cover),
        "-p", str(threads),
    ]
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def two_tier_families(tier1_cluster_tsv: str | Path, tier2_cluster_tsv: str | Path) -> dict:
    """Combine tier-1 family membership with the tier-2 superfamily label
    for each tier-1 representative. Returns
    {tier1_rep: {"members": [...], "superfamily": tier2_rep_of_tier1_rep}}."""
    tier1_families = build_families(read_cluster_tsv(tier1_cluster_tsv))
    tier1_rep_to_tier2_rep = read_cluster_tsv(tier2_cluster_tsv)
    return {
        rep: {
            "members": members,
            "superfamily": tier1_rep_to_tier2_rep.get(rep, rep),
        }
        for rep, members in tier1_families.items()
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    mmseqs_defaults = {"mmseqs-tier1": (0.9, 0.8, True), "mmseqs-tier2": (0.4, 0.8, False)}
    for name, (default_id, default_cov, default_reassign) in mmseqs_defaults.items():
        p = sub.add_parser(name)
        p.add_argument("--fasta", required=True, type=Path)
        p.add_argument("--out_prefix", required=True)
        p.add_argument("--min_seq_id", type=float, default=default_id)
        p.add_argument("--cov", type=float, default=default_cov)
        p.add_argument("--cluster_reassign", type=int, default=int(default_reassign),
                        help="1/0 -- mmseqs --cluster-reassign (default matches the "
                             "tier's original study default)")
        p.add_argument("--threads", type=int, default=1)

    diamond_defaults = {"diamond-tier1": (90, 80), "diamond-tier2": (40, 80)}
    for name, (default_approx_id, default_member_cover) in diamond_defaults.items():
        p = sub.add_parser(name)
        p.add_argument("--fasta", required=True, type=Path)
        p.add_argument("--out_prefix", required=True)
        p.add_argument("--approx_id", type=float, default=default_approx_id)
        p.add_argument("--member_cover", type=float, default=default_member_cover)
        p.add_argument("--threads", type=int, default=1)

    args = ap.parse_args()

    if args.command in mmseqs_defaults:
        run_mmseqs_cluster(
            args.fasta, args.out_prefix, min_seq_id=args.min_seq_id, cov=args.cov,
            cluster_reassign=bool(args.cluster_reassign), threads=args.threads,
        )
    else:
        run_diamond_cluster(
            args.fasta, args.out_prefix, approx_id=args.approx_id,
            member_cover=args.member_cover, threads=args.threads,
        )


if __name__ == "__main__":
    main()
