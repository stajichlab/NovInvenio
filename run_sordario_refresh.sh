#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out logs/sordario_refresh.log

module load nextflow

# Fresh (non -resume) re-run of the sordario novelty_discovery two-phase pipeline, to pick up
# the presence-calling fixes merged in PR #59 (H2 paralog-cutoff fix a49ee0c,
# the circular family-calibration fix, 2026-09-03). The original results/sordario run predates that
# fix. Not using -resume: workflows/novelty_discovery.nf's NOVELTY_PRESENCE_MATRIX process's script was NOT
# changed by the fix commits (same CLI flags, same rendered script text), so
# a plain -resume would silently reuse the stale pre-fix cached output rather
# than recomputing -- see .living/learnings.md's 2026-09-04 entry on this
# caching gotcha, and its 2026-09-05 entry on not editing this same working
# tree while a job like this is running.
#
# --project sordario is kept so this reuses results/sordario/search_cache/
# (storeDir is keyed by outdir/project/search_cache), so the expensive
# phmmer pairwise + self searches are NOT re-run -- only
# BUILD_PRESENCE_MATRIX and everything downstream recomputes.
#
# Isolated launch dir (.nf_launch/sordario_refresh/) so this session/work/ dir
# cannot collide with any other concurrently-running Nextflow session.

set -euo pipefail

REPO_ROOT=$(pwd)
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

CONFIG=$(abspath configs/sordario.csv)
DATA_DIR=$(abspath data)
PFAM=$(abspath db/pfam/38.2/Pfam-A.hmm)
SWISSPROT=$(abspath db/uniprot/uniprot_sprot.fasta.dmnd)
MODELORGS=$(abspath configs/modelorgs.yaml)
OUTDIR=$(abspath results)
SITE_CONFIG=$(abspath conf/ucr_hpcc_slurm.config)

LAUNCH_DIR="$REPO_ROOT/.nf_launch/sordario_refresh"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

nextflow run "$REPO_ROOT/main.nf" \
    --config "$CONFIG" \
    --data_dir "$DATA_DIR" \
    --run_tool phmmer \
    --cluster_tool novelty_discovery \
    --pfam_hmm "$PFAM" \
    --swissprot_dmnd "$SWISSPROT" \
    --modelorgs_config "$MODELORGS" \
    --outdir "$OUTDIR" \
    --project sordario \
    -profile slurm \
    -c "$SITE_CONFIG"
