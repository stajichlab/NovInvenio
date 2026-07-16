# NovInvenio

Identifies lineage-specific (novel) genes: proteins present in ≥N% of an ingroup but absent from all outgroup proteomes *and* outgroup genomes.

## Pipeline overview

```
Config CSV
    │
    ├─► SEARCH   — pairwise proteome searches (phmmer / diamond / BLAST)
    │             + self-vs-self search per ingroup species → paralog_cutoffs.tsv
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
    ├─► SUMMARIZE — per-species novelty files
    │               → novelties.<SHORT>.tsv for each ingroup species
    │
    └─► REPORT   — interactive HTML report of the candidates and the matrix
                    → report.html (self-contained; opens offline in a browser)
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
    --data_dir data \
    --run_tool diamond \
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
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
| `--evalue` | `1e-5` | Fallback e-value for proteins with no detectable within-proteome paralog |
| `--parse_evalue` | `0.01` | Noise ceiling applied when parsing raw pairwise hits; the authoritative per-gene cutoff comes from the self-vs-self paralog search |
| `--ingroup_min_frac` | `0.75` | Min fraction of ingroup proteomes that must contain a hit |
| `--pfam_hmm` | `null` | Path to Pfam-A.hmm; skips Pfam annotation if unset |
| `--swissprot_dmnd` | `null` | Path to SwissProt `.dmnd` database; skips if unset |
| `--modelorgs_config` | `null` | YAML listing model organisms for gene name lookup (see `configs/modelorgs.yaml`) |
| `--report_sequences` | `novelties` | Which proteins carry a sequence in `report.html`: `novelties`, `all`, or `none`. Sequences dominate the file size |
| `--project` | *(auto)* | Output subdirectory name; defaults to config CSV basename |
| `--outdir` | `results` | Root output directory |
| `--hmm_mpi` | `false` | Run hmmsearch in MPI mode (requires MPI-enabled HMMER) |
| `--hmm_mpi_tasks` | `null` | Number of MPI tasks; defaults to `max_cpus` when `--hmm_mpi true` |

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
| `self_hits/<SHORT>.paralog_cutoffs.tsv` | Per-gene paralog e-value cutoffs from self-vs-self search (columns: `protein_ID`, `paralog_protein_ID`, `bitscore`, `evalue`) |
| `presence_matrix.tsv` | Gene × proteome presence/absence matrix (presence scored against paralog cutoffs) |
| `candidates.txt` | Protein IDs meeting ingroup/outgroup criteria |
| `candidates.fa` | FASTA of all candidate proteins |
| `clusters/` | mmseqs2 cluster output (TSV, rep FASTA, all-seqs FASTA) |
| `tblastn/` | Per-outgroup TBLASTN TSVs |
| `tblastn_summary.tsv` | Protein × outgroup genome hit matrix |
| `candidates.pfam.tblout` | Pfam hmmsearch tblout (if `--pfam_hmm` set) |
| `candidates.swissprot.tsv` | SwissProt diamond hits (if `--swissprot_dmnd` set) |
| `presence_matrix.function.tsv` | Annotated matrix (gene_name, Best_Swissprot, Pfam_Names) |
| `novelties.<SHORT>.tsv` | Per-species novelty candidates (one file per ingroup species) |
| `report.html` | Interactive report — browse the candidates and the presence matrix ([below](#interactive-report-reporthtml)) |

### novelties.\<SHORT\>.tsv columns

Each file contains proteins that are:
- Present in species `<SHORT>`
- Present in ≥ `ingroup_min_frac` of all ingroup species
- Absent from all outgroup proteomes (protein-level search)
- Absent from all outgroup genomes (TBLASTN validation)

Columns include: `protein_id`, `source_proteome`, per-species presence columns,
`gene_name`, `product_description`, `function_source`, `Best_Swissprot`,
`Pfam_Names`.

## Interactive report (report.html)

The pipeline writes `results/<project>/report.html` — a single file for browsing the
novelty candidates and the presence/absence matrix. It has no external dependencies and
makes no network requests, so copy it anywhere and open it directly:

```bash
# from your laptop
scp cluster:/path/to/NovInvenio/results/pezio4_asco/report.html .
open report.html          # macOS  (Linux: xdg-open, or just drag into a browser)
```

No web server is needed — `file://` works, including offline.

### What you can do in it

- **Heatmap** — one row per protein, one column per proteome. Blue marks presence from
  the protein search; green marks a TBLASTN hit in an outgroup genome. Ingroup, outgroup
  and genome columns are separate labelled blocks, so a clean candidate reads as a solid
  blue block on the left with nothing to the right. Hover a cell for the species and the
  evidence; click a row (or use ↑/↓) to pin it in the detail panel.
- **Detail panel** — annotation plus links out: Pfam domains → InterPro, SwissProt →
  UniProt and AlphaFold, model-organism genes → FungiDB, and BLASTP at NCBI with the
  sequence pre-filled. Copy FASTA copies the single protein.
- **Filters** — search (ID, gene, product, Pfam), source proteome, annotation source,
  novelty-only, has-Pfam, no-TBLASTN-hit, and sorting. They scope both tabs at once.
- **Table tab** — the same rows with every value as text, for reading without hovering.
- **Download TSV / FASTA** — exports whatever the filters currently select, so you can
  take a subset straight into another tool.

> **Tip:** `SUMMARIZE` currently runs with `--skip_tblastn_filter`, so proteins with
> TBLASTN hits in outgroup genomes are *kept* and flagged rather than dropped. Tick
> **No TBLASTN hit** to see the subset that is also absent at the nucleotide level.

### Regenerating without re-running the pipeline

Useful after re-annotating, or to embed every sequence rather than just the novelties:

```bash
bin/make_report.py \
    --matrix results/pezio4_asco/presence_matrix.function.tsv \
    --config configs/pezio4_asco.csv \
    --tblastn_summary results/pezio4_asco/tblastn_summary.tsv \
    --novelties results/pezio4_asco/novelties.*.tsv \
    --candidates_fa results/pezio4_asco/candidates.fa \
    --output results/pezio4_asco/report.html
```

`--tblastn_summary`, `--novelties` and `--candidates_fa` are all optional: without
`--novelties` the novelty status is recomputed from the matrix (it reproduces the
published tables exactly), and without `--tblastn_summary` the genome columns are simply
omitted.

`--sequences` controls which rows carry a sequence, and sequences are the bulk of the
file (~4.9 MB with them, ~2.2 MB without). They are read from `--candidates_fa` and from
the `protein_sequence` column of the novelties tables, so `all` only embeds more than
`novelties` when the novelty set is a strict *subset* of `candidates.fa` — that is, when
the TBLASTN filter is active. With `--skip_tblastn_filter` wired in `SUMMARIZE` today the
two sets are identical, so `all` and `novelties` give the same result. `--sequences none`
drops them entirely: the BLASTP link, Copy FASTA and the sequence view disappear from the
detail panel, and Download FASTA reports that it has nothing to write.

Gene names and Pfam/SwissProt annotation are read from the **matrix**, so if the report
shows an empty Gene column, re-run `ANNOTATE` with `--modelorgs_config` rather than
looking for the problem in the report.

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

### MPI mode for hmmsearch

When running the downstream `HMMSEARCH` process on SLURM with an MPI-enabled
HMMER build (e.g. from conda-forge), pass `--hmm_mpi true` to use
`mpirun -np N hmmsearch --mpi` instead of `hmmsearch --cpu N`.  The SLURM
profile automatically requests `--ntasks N` instead of `--cpus-per-task N`.

```bash
nextflow run main.nf \
    -profile slurm \
    --config configs/pezio4_asco.csv \
    --data_dir /bigdata/stajichlab/shared/data/fungi/proteomes \
    --run_tool diamond \
    --hmm_mpi true \
    --hmm_mpi_tasks 16
```

If `--hmm_mpi_tasks` is omitted it defaults to `--max_cpus` (32 by default).

### Database paths

The pipeline resolves `--pfam_hmm`, `--swissprot_dmnd`, and
`--modelorgs_config` to absolute paths at startup and automatically creates a
`db/` symlink inside every task work directory.  You can pass relative paths
(e.g. `db/pfam/Pfam-A.hmm`) from the project root and they will resolve
correctly regardless of where Nextflow stages the task.

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
