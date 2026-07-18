module load nextflow

nextflow run main.nf -resume --config configs/pezio5_fungi.csv --data_dir data --run_tool phmmer --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
	--swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd --modelorgs_config configs/modelorgs.yaml -profile slurm --project pezizo5

