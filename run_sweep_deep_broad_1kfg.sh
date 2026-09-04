#!/usr/bin/bash
#SBATCH -p batch -c 4 --mem 16gb --time=6-00:00:00 --out logs/sweep_deep_broad_1kfg.log
module load nextflow

# Deep-divergence, large-panel coverage-floor sweep (2026-09-03, config 3 per the
# broader-grid plan in todo/validate-hmm-presence-coverage-broader-sweep.md): ingroup =
# 105 species across 7 "filamentous Ascomycota" classes (Sordariomycetes/Eurotiomycetes/
# Dothideomycetes/Leotiomycetes/Pezizomycetes/Lecanoromycetes/Xylonomycetes), outgroup =
# stratified random sample (5/phylum, seed 42) across 11 other fungal phyla. 130
# proteomes total -- real ADR-0002/Chaetothyriales scale, the "hard" end of the
# divergence-depth axis and the main stress test for whether the pezizo5 sweep's
# recommendation (decision #13) generalizes.
#
# --time=6-00:00:00: this driver just orchestrates (bash loop + nextflow head process);
# the actual mmseqs clustering / family-HMM build / hmmsearch scatter-gather work is
# dispatched as separate SLURM sub-jobs via -profile slurm and governed by
# nextflow.config's own per-process time limits (up to --max_time, 24h). Given up to 128
# grid points could in principle run here (though MIN_SEQ_IDS_LIST/COVS_LIST/
# HMM_EVALUES_LIST below scope it to the same 4-point hmm_cov x hmm_residues grid as
# pezizo5) and each point's own nextflow run can itself take hours at this proteome
# count, this generously exceeds what a single scoped 4-point sweep should need (batch
# partition allows up to 30 days; DefaultTime alone is already 7 days) so the driver
# itself is never the bottleneck.

export CONFIG=configs/deep_broad_1kfg.csv
export DATA_DIR=data
export CLADE=deep_broad_1kfg
export PFAM=db/pfam/38.2/Pfam-A.hmm
export SWISSPROT=db/uniprot/uniprot_sprot.fasta.dmnd
export MODELORGS=configs/modelorgs.yaml
export CONTROLS=configs/controls/deep_broad_1kfg.controls.csv
export PROFILE="-profile slurm"
export OUTDIR=results
export SWEEP_DIR=results/sweep_deep_broad_1kfg

# BUSCO tables derived at runtime from the config CSV + busco_1kfg/ outputs (one batch
# array job BUSCO-ran all 130 species in this config).
BUSCO_TABLES=""
BUSCO_OUTGROUP_TABLES=""
while IFS=, read -r group species strain protein dna short taxongroup; do
  [ "$group" = "GROUP" ] && continue
  table="busco_1kfg/${short}.busco/run_fungi_odb12/full_table.tsv"
  if [ "$group" = "IN" ]; then
    BUSCO_TABLES="$BUSCO_TABLES ${short}=${table}"
  else
    BUSCO_OUTGROUP_TABLES="$BUSCO_OUTGROUP_TABLES ${short}=${table}"
  fi
done < "$CONFIG"
export BUSCO_TABLES
export BUSCO_OUTGROUP_TABLES

export MIN_SEQ_IDS_LIST="0.3"
export COVS_LIST="0.8"
export HMM_EVALUES_LIST="1e-3"
export HMM_COVS_LIST="0.5 0.3"
export MIN_RESIDUES_LIST="0 100"

bash bin/run_param_sweep.sh
