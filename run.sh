module load nextflow singularity

# Shared SIF cache so concurrent SLURM array tasks reuse one pulled/converted
# image instead of each task pulling its own copy. Only set a default if the
# user hasn't already configured one in their environment.
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/jstajich/singularity_cache}"

nextflow run main.nf -resume --config configs/nirr.csv --data_dir data --run_tool phmmer --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
	--swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd --modelorgs_config configs/modelorgs.yaml \
	-profile slurm,singularity -c conf/ucr_hpcc_slurm.config --project neolecta

