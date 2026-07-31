#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --out logs/sordario_pairwise_diamond.log
module load nextflow

nextflow run main.nf --config configs/sordario_pairwise.csv --data_dir data --run_tool diamond \
	--cluster_tool pairwise -profile slurm \
	--pfam_hmm db/pfam/38.2/Pfam-A.hmm \
	--swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd --project sordario_pairwise_diamond
