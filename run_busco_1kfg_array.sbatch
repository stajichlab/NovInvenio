#!/usr/bin/bash -l
#SBATCH -p short -c 8 --mem 16gb --time=2:00:00 --array=1-130%25
#SBATCH --job-name=busco1kfg
#SBATCH --output=logs/busco_1kfg_%a.log
# Array job: BUSCO fungi_odb12 for the 130-species union of configs/sordariales_shallow_1kfg.csv
# and configs/deep_broad_1kfg.csv (2026-09-03 broader-grid validation, see
# todo/validate-hmm-presence-coverage-broader-sweep.md). One task per line of
# config_support/1KFG_busco_species.tsv (Short<TAB>ProteinFilename).
module load busco

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" config_support/1KFG_busco_species.tsv)
SHORT=$(echo "$LINE" | cut -f1)
PROTEIN=$(echo "$LINE" | cut -f2)

if [ -z "$SHORT" ]; then
  echo "No line ${SLURM_ARRAY_TASK_ID} in species list" >&2
  exit 1
fi

busco -m proteins -l fungi_odb12 -c 8 --offline \
    -i "data/pep/${PROTEIN}" -o "${SHORT}.busco" \
    --out_path busco_1kfg \
    --download_path /srv/projects/db/BUSCO/v12
