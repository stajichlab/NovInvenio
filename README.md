# NovInvenio

Identifies lineage-specific (novel) genes: proteins present in ≥N% of an ingroup but absent from all outgroup proteomes *and* outgroup genomes.

## Pipeline overview

```
Config CSV
    │
    ├─► SEARCH   — pairwise proteome searches (phmmer / diamond / BLAST)
    │               → presence_matrix.tsv, candidates.txt
    │
    ├─► CLUSTER  — mmseqs2 easy-cluster of candidate proteins
    │               → candidates.fa, cluster reps, clusters_cluster.tsv
    │
    ├─► VALIDATE — TBLASTN of cluster reps vs outgroup genomes
    │               → tblastn_summary.tsv (protein × genome hit matrix)
    │
    ├─► ANNOTATE — Pfam hmmsearch + SwissProt diamond + model organism gene names
    │               → presence_matrix.function.tsv
    │                 (adds gene_name, product_description, Best_Swissprot, Pfam_Names)
    │
    └─► SUMMARIZE — per-species novelty files
                    → novelties.<SHORT>.tsv for each ingroup species
```

## Quick start

### 1. Set up databases (one time)

```bash
bash pipeline/00_setup.sh
```

This downloads Pfam-A, SwissProt (diamond-formatted), and any model-organism
gene name files into `db/`.

### 2. Configure model organism gene names (optional)

Edit `configs/modelorgs.yaml` to list the ingroup species you want gene names
for.  Each entry maps a `Short` identifier (from your analysis CSV) to a gene
names file and describes how protein IDs map to gene IDs.  See the comments in
`configs/modelorgs.yaml` for details and examples.

### 3. Run the pipeline

```bash
nextflow run main.nf \
    --config configs/pezio4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --pfam_hmm db/pfam/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

### 4. Resume after failure

```bash
nextflow run main.nf -resume --config configs/... --data_dir ...
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--config` | *(required)* | Analysis config CSV (see format below) |
| `--data_dir` | *(required)* | Directory containing FASTA files listed in the CSV |
| `--run_tool` | `phmmer` | Search tool: `phmmer`, `diamond`, or `blast` |
| `--evalue` | `1e-5` | E-value cutoff for all search steps |
| `--ingroup_min_frac` | `0.75` | Min fraction of ingroup proteomes that must contain a hit |
| `--pfam_hmm` | `null` | Path to Pfam-A.hmm; skips Pfam annotation if unset |
| `--swissprot_dmnd` | `null` | Path to SwissProt `.dmnd` database; skips if unset |
| `--modelorgs_config` | `null` | YAML listing model organisms for gene name lookup (see `configs/modelorgs.yaml`) |
| `--outdir` | `results` | Root output directory |

## Config CSV format

```csv
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
OUT,Neolecta irregularis,,Neolecta_irregularis.proteins.fa,Neolecta_irregularis.scaffolds.fa,Nirr,Taphrinomycotina
IN,Neurospora crassa,OR74A,Ncrassa.pep.fa,Ncrassa.dna.fa,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af293,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
```

- `GROUP`: `IN` (ingroup) or `OUT` (outgroup)
- `Short`: ≤8-char unique identifier used in all output filenames
- `Protein`, `DNA`: FASTA basenames resolved relative to `--data_dir`
- Config filename (without `.csv`) becomes the results subdirectory

## Model organism config (configs/modelorgs.yaml)

Gene name lookup is driven by a separate YAML file.  This decouples annotation
from the analysis CSV and lets you reuse the same model organism config across
multiple runs.

```yaml
model_organisms:
  # FungiDB organism requiring diamond mapping to cross-reference IDs
  - short: Ncra
    gene_names_csv: db/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv
    id_transform: diamond_fasta
    diamond_hits: db/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv
    protein_fasta: db/modelorgs/FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
    fasta_gene_field: gene

  # FungiDB organism where protein ID → gene ID via suffix stripping
  - short: Afum
    gene_names_csv: db/modelorgs/Aspergillus_fumigatus_gene_names_FungiDB.csv
    id_transform: strip
    strip_pattern: "-T-p\\d+$"
```

**Supported `id_transform` values:**

| Value | Description |
|---|---|
| `direct` | Protein ID is used as the gene lookup key unchanged |
| `strip` | Regex `strip_pattern` is removed from protein ID to get gene key |
| `diamond_fasta` | protein_id → diamond_hits → ref protein_id → FASTA header gene field |

Gene name CSV column names default to FungiDB format (`Gene ID`,
`Gene Name or Symbol`, `Product Description`) and can be overridden per entry
with `gene_id_col`, `gene_name_col`, `product_col`.  Set `csv_delimiter: "\t"`
for tab-separated files.

## Outputs

All outputs are written to `results/<project>/`:

| File | Description |
|---|---|
| `presence_matrix.tsv` | Gene × proteome presence/absence matrix |
| `candidates.txt` | Protein IDs meeting ingroup/outgroup criteria |
| `candidates.fa` | FASTA of all candidate proteins |
| `clusters/` | mmseqs2 cluster output (TSV, rep FASTA, all-seqs FASTA) |
| `tblastn/` | Per-outgroup TBLASTN TSVs |
| `tblastn_summary.tsv` | Protein × outgroup genome hit matrix |
| `candidates.pfam.tblout` | Pfam hmmsearch tblout (if `--pfam_hmm` set) |
| `candidates.swissprot.tsv` | SwissProt diamond hits (if `--swissprot_dmnd` set) |
| `presence_matrix.function.tsv` | Annotated matrix (gene_name, Best_Swissprot, Pfam_Names) |
| `novelties.<SHORT>.tsv` | Per-species novelty candidates (one file per ingroup species) |

### novelties.\<SHORT\>.tsv columns

Each file contains proteins that are:
- Present in species `<SHORT>`
- Present in ≥ `ingroup_min_frac` of all ingroup species
- Absent from all outgroup proteomes (protein-level search)
- Absent from all outgroup genomes (TBLASTN validation)

Columns include: `protein_id`, `source_proteome`, per-species presence columns,
`gene_name`, `product_description`, `function_source`, `Best_Swissprot`,
`Pfam_Names`.

## Running on SLURM

```bash
nextflow run main.nf \
    -profile slurm \
    --config configs/pezio4_asco.csv \
    --data_dir /bigdata/stajichlab/shared/data/fungi/proteomes \
    --run_tool diamond \
    --ingroup_min_frac 0.8 \
    --pfam_hmm db/pfam/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

## Development

See `CLAUDE.md` for architecture details, module/workflow conventions, and
how to add new analysis stages.

### Running tests

```bash
pytest tests/ -v
```

### Built-in test profile

```bash
nextflow run main.nf -profile test
```
