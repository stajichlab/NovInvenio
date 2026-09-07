#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out logs/zoosporic_opisthokont_dikarya.log

module load nextflow

# First run of the zoosporic-fungi-vs-animals-vs-Dikarya novelty/loss study
# (configs/zoosporic_opisthokont_dikarya.csv, batch spec
# configs/batches/zoosporic_opisthokont_dikarya.yaml): tests for genes shared
# between zoosporic fungi (Batrachochytrium dendrobatidis, Spizellomyces
# punctatus, Chytriomyces hyalinus, Synchytrium endobioticum, Catenaria
# anguillulae, Paraphysoderma sedebokerense) and animals/their closest
# non-animal relative (Drosophila melanogaster, Caenorhabditis elegans, Mus
# musculus, Monosiga brevicollis), but absent from Dikarya (Saccharomyces
# cerevisiae, Neurospora crassa, Aspergillus nidulans, Schizosaccharomyces
# pombe, Coprinopsis cinerea).
#
# --data_dir points at the per-batch link-dir (data/zoosporic_opisthokont_dikarya/),
# built by bin/build_targeted_configs.py's --link-dir -- NOT the flat
# top-level data/ the other run_*_refresh.sh scripts use, since this config's
# Protein/DNA columns are basenames resolved against that link-dir's pep/
# and dna/ subdirs, not the repo-wide flat layout.
#
# Isolated launch dir (.nf_launch/zoosporic_opisthokont_dikarya/) so this
# session/work/ dir cannot collide with any other concurrently-running
# Nextflow session.

set -euo pipefail

REPO_ROOT=$(pwd)
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

CONFIG=$(abspath configs/zoosporic_opisthokont_dikarya.csv)
DATA_DIR=$(abspath data/zoosporic_opisthokont_dikarya)
PFAM=$(abspath db/pfam/38.2/Pfam-A.hmm)
SWISSPROT=$(abspath db/uniprot/uniprot_sprot.fasta.dmnd)
MODELORGS=$(abspath configs/modelorgs.yaml)
OUTDIR=$(abspath results)
SITE_CONFIG=$(abspath conf/ucr_hpcc_slurm.config)

LAUNCH_DIR="$REPO_ROOT/.nf_launch/zoosporic_opisthokont_dikarya"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

nextflow run "$REPO_ROOT/main.nf" \
    --config "$CONFIG" \
    --data_dir "$DATA_DIR" \
    --run_tool diamond \
    --pfam_hmm "$PFAM" \
    --swissprot_dmnd "$SWISSPROT" \
    --modelorgs_config "$MODELORGS" \
    --cluster_tool pairwise \
    --outdir "$OUTDIR" \
    --project zoosporic_opisthokont_dikarya \
    -profile slurm \
    -c "$SITE_CONFIG"
