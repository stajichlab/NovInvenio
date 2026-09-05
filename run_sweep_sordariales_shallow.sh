#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --time=2-00:00:00 --out logs/sweep_sordariales_shallow.log
module load nextflow

# "Shallow divergence" broader-grid validation point (2026-09-03), alongside the pezizo5
# (medium) and deep_broad_1kfg (deep+large) configs -- see
# todo/validate-hmm-presence-coverage-broader-sweep.md and .living/decisions.md #12/#13.
# Ingroup = Sordariales (Chaetomium/Neurospora/Podospora/Sordaria, 4 species, one order);
# outgroup = sibling Sordariomycetes orders (Hypocreales/Glomerellales/Coniochaetales/
# Diaporthales/Xylariales/Ophiostomatales, capped 3/order). 14 species total -- real
# orthologs should mostly still align well at this shallow a divergence, so this is
# partly a sanity check that loosening hmm_cov/hmm_residues doesn't start admitting
# false positives once coverage would already have been fine.
#
# BUSCO_TABLES/BUSCO_OUTGROUP_TABLES built dynamically from the config CSV itself
# (busco_1kfg/<Short>.busco/, from the 130-species union BUSCO array job, issue: see
# run_busco_1kfg_array.sbatch) so there's nothing to keep in sync by hand.

export CONFIG=configs/sordariales_shallow_1kfg.csv
export DATA_DIR=data
export CLADE=sordariales_shallow
export PFAM=db/pfam/38.2/Pfam-A.hmm
export SWISSPROT=db/uniprot/uniprot_sprot.fasta.dmnd
export MODELORGS=configs/modelorgs.yaml
export CONTROLS=configs/controls/sordariales_shallow.controls.csv
export PROFILE="-profile slurm"
export OUTDIR=results
export SWEEP_DIR=results/sweep_sordariales_shallow

export BUSCO_TABLES=$(awk -F',' 'NR>1 && $1=="IN"{printf "%s=busco_1kfg/%s.busco/run_fungi_odb12/full_table.tsv ", $6, $6}' "$CONFIG")
export BUSCO_OUTGROUP_TABLES=$(awk -F',' 'NR>1 && $1=="OUT"{printf "%s=busco_1kfg/%s.busco/run_fungi_odb12/full_table.tsv ", $6, $6}' "$CONFIG")

export MIN_SEQ_IDS_LIST="0.3"
export COVS_LIST="0.8"
export HMM_EVALUES_LIST="1e-3"
export HMM_COVS_LIST="0.5 0.3"
export MIN_RESIDUES_LIST="0 100"

bash bin/run_param_sweep.sh
