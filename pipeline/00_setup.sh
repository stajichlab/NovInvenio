#!/usr/bin/bash -l
#SBATCH -c 4 --mem 8gb 
# Setup Swissprot, Pfam database

module load hmmer
module load ncbi-blast
module load diamond
DB=db
mkdir -p $DB

PFAM=pfam
PFAM_VERSION=38.2
URL=https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
FNAME=$(basename $URL .gz)
mkdir -p $DB/$PFAM/${PFAM_VERSION}
TARGET=$DB/$PFAM/${PFAM_VERSION}/${FNAME}
if [ ! -f $TARGET ]; then
	curl $URL | pigz -dc > $TARGET
	hmmpress $TARGET
fi

SWISSPROT=uniprot
URL=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mkdir -p $DB/$SWISSPROT
FNAME=$(basename $URL .gz)
TARGET=$DB/$SWISSPROT/$FNAME
if [ ! -f $TARGET ]; then
	curl $URL | pigz -dc > $TARGET
	makeblastdb -dbtype prot -title swissprot -in $TARGET
	diamond makedb  --in $TARGET --db $TARGET.dmnd
fi

MODELORG=modelorgs
mkdir -p $DB/$MODELORG
pushd $DB/$MODELORG

URLS=(
    https://fungidb.org/a/service/raw-files/release-68/NcrassaOR74A/fasta/data/FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
    https://fungidb.org/a/service/raw-files/release-68/AfumigatusAf293/fasta/data/FungiDB-68_AfumigatusAf293_AnnotatedProteins.fasta
    https://fungidb.org/a/service/raw-files/release-68/ScerevisiaeS288C/fasta/data/FungiDB-68_ScerevisiaeS288C_AnnotatedProteins.fasta
    https://fungidb.org/a/service/raw-files/release-68/Spombe972h/fasta/data/FungiDB-68_Spombe972h_AnnotatedProteins.fasta

)

for url in "${URLS[@]}"; do
    fname=$(basename "$url")
    if [ ! -f "$fname" ]; then
        curl -L -O "$url"
    fi
done

# Build diamond databases for FungiDB protein sets
for fasta in FungiDB-*_AnnotatedProteins.fasta; do
    if [ ! -f "${fasta}.dmnd" ]; then
        diamond makedb --in "$fasta" --db "$fasta" --threads "$SLURM_CPUS_ON_NODE"
    fi
done

# Map project Ncra proteins (JGI IDs) to FungiDB NCU IDs
NCRA_QUERY=../../data/pep/Neurospora_crassa_OR74A.proteins.fa
NCRA_DB=FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
NCRA_HITS=Ncra_vs_FungiDB_Ncra.diamond.tsv
if [ ! -f "$NCRA_HITS" ]; then
    diamond blastp \
        --query "$NCRA_QUERY" \
        --db "${NCRA_DB}.dmnd" \
        --outfmt 6 qseqid sseqid evalue bitscore pident \
        --evalue 1e-5 \
        --max-target-seqs 1 \
        --threads "${SLURM_CPUS_ON_NODE:-4}" \
        --out "$NCRA_HITS"
fi

popd
