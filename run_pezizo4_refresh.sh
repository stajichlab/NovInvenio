#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out logs/pezizo4_refresh.log

module load nextflow

# Fresh (non -resume) re-run of the pezizo4 pairwise novelty/loss pipeline, to pick up
# the presence-calling fixes merged in PR #59 (H2 paralog-cutoff fix a49ee0c,
# 2026-09-03). The original results/pezizo4 run predates that
# fix. Not using -resume: modules/build_presence_matrix.nf's process script was NOT
# changed by the fix commits (same CLI flags, same rendered script text), so
# a plain -resume would silently reuse the stale pre-fix cached output rather
# than recomputing -- see .living/learnings.md's 2026-09-04 entry on this
# caching gotcha, and its 2026-09-05 entry on not editing this same working
# tree while a job like this is running.
#
# --project pezizo4 is kept so this reuses results/pezizo4/search_cache/
# (storeDir is keyed by outdir/project/search_cache), so the expensive
# phmmer pairwise + self searches are NOT re-run -- only
# BUILD_PRESENCE_MATRIX and everything downstream recomputes.
#
# Isolated launch dir (.nf_launch/pezizo4_refresh/) so this session/work/ dir
# cannot collide with any other concurrently-running Nextflow session.

set -euo pipefail

REPO_ROOT=$(pwd)
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

CONFIG=$(abspath configs/pezizo4_asco.csv)
DATA_DIR=$(abspath data)
PFAM=$(abspath db/pfam/38.2/Pfam-A.hmm)
SWISSPROT=$(abspath db/uniprot/uniprot_sprot.fasta.dmnd)
MODELORGS=$(abspath configs/modelorgs.yaml)
OUTDIR=$(abspath results)
SITE_CONFIG=$(abspath conf/ucr_hpcc_slurm.config)

LAUNCH_DIR="$REPO_ROOT/.nf_launch/pezizo4_refresh"
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
    --project pezizo4 \
    -profile slurm \
    -c "$SITE_CONFIG"
