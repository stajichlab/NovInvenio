#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --out logs/sweep_pezizo5_coverage.log
module load nextflow

# Scoped hmm_presence_cov / hmm_presence_min_residues sweep (2026-09-03 follow-up to
# ADR-0002 Q8, decision #12 in .living/decisions.md). Holds the clustering dimensions
# (family_min_seq_id, family_cov, hmm_presence_evalue) fixed at the already-shipped
# defaults (0.3 / 0.8 / 1e-3) rather than re-sweeping the full ADR-0002 grid -- see
# bin/run_param_sweep.sh's grid-size NOTE. 2x2 = 4 grid points, not 128.
#
# NOTE: no BUSCO_OUTGROUP_TABLES are available yet for pezizo5 (only the 5 IN-group
# species have been BUSCO-run, under busco_pezizo5/) -- presence_recovery, the metric
# that actually targets hmm_presence_cov/hmm_presence_min_residues, will come back
# unmeasured this round. busco_recovery (clustering quality, insensitive to these two
# params) and n_novelties/n_families still get reported. Running BUSCO on an OUTGROUP
# proteome (e.g. Spom, Scer) and re-running with BUSCO_OUTGROUP_TABLES set is a
# worthwhile follow-up before trusting the recommended default this run prints.

export CONFIG=configs/pezizo5.csv
export DATA_DIR=data
export CLADE=pezizo5_coverage
export PFAM=db/pfam/38.2/Pfam-A.hmm
export SWISSPROT=db/uniprot/uniprot_sprot.fasta.dmnd
export MODELORGS=configs/modelorgs.yaml
export CONTROLS=configs/controls/pezizo5.controls.csv
export PROFILE="-profile slurm"
export OUTDIR=results
export SWEEP_DIR=results/sweep_pezizo5_coverage
export BUSCO_TABLES="Amega=busco_pezizo5/Amega.busco/run_fungi_odb12/full_table.tsv \
Ncra=busco_pezizo5/Ncra.busco/run_fungi_odb12/full_table.tsv \
Afum=busco_pezizo5/Afum.busco/run_fungi_odb12/full_table.tsv \
Ztri=busco_pezizo5/Ztri.busco/run_fungi_odb12/full_table.tsv \
Cimm=busco_pezizo5/Cimm.busco/run_fungi_odb12/full_table.tsv"
export BUSCO_OUTGROUP_TABLES=""

export MIN_SEQ_IDS_LIST="0.3"
export COVS_LIST="0.8"
export HMM_EVALUES_LIST="1e-3"
export HMM_COVS_LIST="0.5 0.3"
export MIN_RESIDUES_LIST="0 100"

bash bin/run_param_sweep.sh
