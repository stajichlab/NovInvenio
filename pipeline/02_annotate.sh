#!/usr/bin/bash -l
#SBATCH -c 16 --mem 32gb --time 2:00:00 --out logs/annotate_run.log
# Annotate presence_matrix.tsv with functional information.
# Run after the Nextflow pipeline has produced results/<PROJECT>/presence_matrix.tsv
# and results/<PROJECT>/candidates.fa.
#
# Usage: sbatch pipeline/02_annotate.sh <project> [modelorgs_config]
#   e.g. sbatch pipeline/02_annotate.sh neolecta configs/modelorgs.yaml

module load hmmer
module load diamond

PROJECT=${1:?Usage: $0 <project> [modelorgs_config]}
MODELORGS_CONFIG=${2:-configs/modelorgs.yaml}
RESULTS=results/$PROJECT
CANDIDATES=$RESULTS/candidates.fa
MATRIX=$RESULTS/presence_matrix.tsv

DB=db
PFAM_HMM=$(ls $DB/pfam/*/Pfam-A.hmm 2>/dev/null | head -1)
SPROT_DMND=$DB/uniprot/uniprot_sprot.fasta.dmnd

if [ ! -f "$CANDIDATES" ]; then
    echo "ERROR: $CANDIDATES not found — run the Nextflow pipeline first." >&2
    exit 1
fi
if [ ! -f "$MATRIX" ]; then
    echo "ERROR: $MATRIX not found — run the Nextflow pipeline first." >&2
    exit 1
fi

# --- Pfam hmmsearch ---
PFAM_HITS=$RESULTS/candidates.pfam.tblout
if [ ! -f "$PFAM_HITS" ] && [ -n "$PFAM_HMM" ]; then
    hmmsearch \
        --tblout "$PFAM_HITS" \
        --noali \
        --cpu "$SLURM_CPUS_ON_NODE" \
        "$PFAM_HMM" \
        "$CANDIDATES" \
        > /dev/null
fi

# --- SwissProt diamond blastp ---
SPROT_HITS=$RESULTS/candidates.swissprot.tsv
if [ ! -f "$SPROT_HITS" ] && [ -f "$SPROT_DMND" ]; then
    diamond blastp \
        --query "$CANDIDATES" \
        --db "$SPROT_DMND" \
        --outfmt 6 qseqid sseqid stitle evalue bitscore \
        --evalue 1e-5 \
        --max-target-seqs 1 \
        --threads "$SLURM_CPUS_ON_NODE" \
        --out "$SPROT_HITS"
fi

# --- Annotate presence matrix ---
MORGS_FLAG=""
[ -f "$MODELORGS_CONFIG" ] && MORGS_FLAG="--modelorgs_config $MODELORGS_CONFIG"

python bin/annotate_presence_matrix.py \
    --matrix "$MATRIX" \
    $MORGS_FLAG \
    ${PFAM_HITS:+--pfam_hits "$PFAM_HITS"} \
    ${SPROT_HITS:+--swissprot_hits "$SPROT_HITS"} \
    --output "$RESULTS/presence_matrix.function.tsv"

echo "Done: $RESULTS/presence_matrix.function.tsv"
