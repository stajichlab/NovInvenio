module load nextflow
nextflow run main.nf -profile slurm --data_dir data --config configs/pezio4_asco.csv --project neolecta -resume
