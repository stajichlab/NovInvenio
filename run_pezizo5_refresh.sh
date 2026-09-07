#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out logs/pezizo5_refresh.log
module load nextflow

# Fresh (non -resume) re-run of the pezizo5 PAIRWISE novelty/loss pipeline
# (--cluster_tool pairwise, the default -- no mmseqs family clustering, independent
# of the deep_broad_1kfg/sordariales_shallow mmseqs sweeps), to pick up the H2
# paralog-cutoff fix (a49ee0c, 2026-09-03). The original results/pezizo5 run
# (2026-07-18/19) predates that fix, and modules/build_presence_matrix.nf's
# process script was NOT touched by the fix commit -- same CLI flags, same
# rendered script text -- so a plain `-resume` would silently reuse the stale
# pre-fix cached presence-matrix output rather than recomputing (see
# .living/learnings.md's 2026-09-04 entry on this exact caching gotcha).
#
# Isolation:
#   - Launched from its own .nf_launch/pezizo5_refresh/ directory (same pattern
#     bin/run_param_sweep.sh uses) so its Nextflow session/work/ dir cannot
#     collide with the two live mmseqs sweeps' isolated .nf_launch/{clade}/
#     sessions, or with the repo-root .nextflow/ history from prior runs.
#   - --outdir/--config/etc. are all converted to absolute paths first so they
#     still resolve correctly after the `cd` into the isolated launch dir.
#   - --project pezizo5 is kept (required to reuse results/pezizo5/search_cache/
#     -- storeDir is keyed by ${outdir}/${project}/search_cache, not by launch
#     dir), so the expensive phmmer pairwise+self searches are NOT re-run --
#     only BUILD_PRESENCE_MATRIX and everything downstream recomputes.
#   - The pre-fix results/pezizo5/ was copied to results/pezizo5_prefix_h2/
#     beforehand so old-vs-new can be diffed once this finishes.
#
# No -resume: intentional, to force BUILD_PRESENCE_MATRIX (and downstream) to
# actually recompute with the corrected code, since its process script is
# unchanged and would otherwise cache-hit the stale pre-fix task.

set -euo pipefail

REPO_ROOT=$(pwd)
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

CONFIG=$(abspath configs/pezizo5_fungi.csv)
DATA_DIR=$(abspath data)
PFAM=$(abspath db/pfam/38.2/Pfam-A.hmm)
SWISSPROT=$(abspath db/uniprot/uniprot_sprot.fasta.dmnd)
MODELORGS=$(abspath configs/modelorgs.yaml)
OUTDIR=$(abspath results)
SITE_CONFIG=$(abspath conf/ucr_hpcc_slurm.config)

LAUNCH_DIR="$REPO_ROOT/.nf_launch/pezizo5_refresh"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

nextflow run "$REPO_ROOT/main.nf" \
    --config "$CONFIG" \
    --data_dir "$DATA_DIR" \
    --run_tool phmmer \
    --pfam_hmm "$PFAM" \
    --swissprot_dmnd "$SWISSPROT" \
    --modelorgs_config "$MODELORGS" \
    --outdir "$OUTDIR" \
    --project pezizo5 \
    -profile slurm \
    -c "$SITE_CONFIG"
