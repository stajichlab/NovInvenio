#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --out logs/sordario.log
module load nextflow

nextflow run main.nf -resume --config configs/sordario.csv --data_dir data --run_tool phmmer \
	--cluster_tool novelty_discovery -profile slurm \
	--pfam_hmm db/pfam/38.2/Pfam-A.hmm \
	--swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd --project sordario
