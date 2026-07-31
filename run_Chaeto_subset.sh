#!/usr/bin/bash
#SBATCH -p batch -c 2 --mem 8gb --out logs/chaeto_subset.log
module load nextflow

nextflow run main.nf -resume --config configs/Chaetothyriales_subset.csv --data_dir data --run_tool phmmer \
	--cluster_tool mmseqs \
	--pfam_hmm db/pfam/38.2/Pfam-A.hmm \
	--swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd -profile slurm --project Chaeto_subset
