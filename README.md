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
    ├─► ANNOTATE — Pfam hmmscan + SwissProt diamond + model organism gene names
    │               → presence_matrix.function.tsv
    │                 (adds gene_name, product_description, Best_Swissprot, Pfam_Names)
    │
    ├─► SUMMARIZE — per-species novelty files
    │               → novelties.<SHORT>.tsv for each ingroup species
    │
    ├─► LOSS_* — a mirror of SEARCH→…→ANNOTATE with the ingroup/outgroup roles
    │             swapped (outgroup as query) → loss_presence_matrix.function.tsv,
    │               loss_tblastn_summary.tsv — candidate lineage-specific gene losses
    │
    └─► REPORT   — self-contained interactive HTML reports (open offline in a browser)
                    → novelties.html, core.html, losses.html
                      then COLLATE_REPORTS copies all three into view/<project>/
                      with a report.html landing page describing the run
```

The `SEARCH → CLUSTER` step (how the presence/absence matrix is built) has **two
interchangeable implementations** selected with `--cluster_tool` — see the next section.
Everything downstream is identical either way.

## Two analysis approaches: pairwise vs family-profile clustering

The presence/absence matrix is the heart of the pipeline, and there are two ways to build
it. Pick one with `--cluster_tool`. Both emit the *same* matrix/candidates contract, so
validation, annotation, and the reports are byte-for-byte the same pipeline downstream.

### `--cluster_tool pairwise` (default) — exact, best for small clades

Searches every ingroup proteome against every other proteome (phmmer/diamond/blast),
calibrating each gene against a self-vs-self paralog search. Exact and well-calibrated,
but the job count grows as `|ingroup| × (N−1)` — great up to a few dozen proteomes, slow
beyond.

```bash
nextflow run main.nf \
    --config configs/pezizo4_asco.csv \
    --data_dir data \
    --run_tool diamond \                 # phmmer | diamond | blast
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

### `--cluster_tool mmseqs` — family profiles, scales to hundreds of proteomes

Clusters the seed group into gene families once (mmseqs), turns each family into a profile
HMM (famsa + hmmbuild), then scans every proteome with the family-HMM database. Presence is
read from family membership plus a uniform HMM hit rule, collapsing the redundancy — roughly
`N` hmmsearch jobs instead of `N²`. Both directions use it: the novelty direction seeds
families from the ingroup, the loss direction from the outgroup.

```bash
nextflow run main.nf \
    --config configs/pezizo5.csv \
    --data_dir data \
    --cluster_tool mmseqs \              # ← the only change from the pairwise run
    --run_tool phmmer \                  # still used for the self-search paralog calibration
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

Family-definition knobs (sensible defaults ship; sweep them for a new clade with
`bin/run_param_sweep.sh`): `--family_min_seq_id` (0.3), `--family_cov` (0.8),
`--hmm_presence_evalue` (1e-3), `--hmm_presence_cov` (0.5), and `--family_chunk_size`
(200 — how many families each parallel profile-build task handles).

### Which should I use?

| | `pairwise` (default) | `mmseqs` |
|---|---|---|
| **Best for** | small, well-defined clades | order-scale configs (many proteomes) |
| **Cost** | `|ingroup| × (N−1)` searches | ~`N` hmmsearch jobs |
| **Presence call** | per-gene paralog-calibrated | family profile HMM |
| **Scales to ~200 proteomes** | slow | yes |

They produce the same matrix contract, so you can run both and compare — `novelties.html`
carries a cross-method `support` column that flags candidates found by both methods (higher
confidence) vs one (threshold-sensitive). Add `-profile slurm` to either to run on the
cluster ([below](#running-on-slurm)).

## Quick start

### Run it now (built-in test data)

The repository ships a tiny four-proteome dataset under `tests/data/`, wired up as a
Nextflow `test` profile, so you can run the whole pipeline end-to-end before downloading
any databases:

```bash
pixi install                       # one-time: create the tool environment
nextflow run main.nf -profile test # runs SEARCH → … → REPORT on tests/data/
```

When it finishes, open the collated landing page and follow the links to the three
interactive reports:

```bash
open view/test/report.html   # macOS; Linux: xdg-open
# view/test/report.html links to copies of novelties.html · core.html · losses.html
```

The test data is deliberately small (a handful of genes each), so the reports will be
nearly empty — the point is to confirm the pipeline runs and the reports render on your
machine. For a real analysis, set up the databases and use one of the `configs/*.csv`
files below.

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
    --config configs/pezizo4_asco.csv \
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

## Running with Docker / Singularity

The pipeline ships a containerised build so you can run it without installing
any tools locally.  The image (`ghcr.io/stajichlab/novinvenio`) contains
`phmmer`, `diamond`, `blast`, `mmseqs2`, `famsa`, `hmmsearch`, `hmmbuild`,
and all other bioinformatics tools used by the pipeline.

### Build the image (one time)

```bash
# Clone the repo first
git clone https://github.com/stajichlab/NovInvenio.git
cd NovInvenio

# Build locally
docker build -t ghcr.io/stajichlab/novinvenio:0.5.0 .

# Push to GitHub Container Registry (requires GHCR write access to stajichlab)
docker push ghcr.io/stajichlab/novinvenio:0.5.0
```

Or pull a pre-built image directly:

```bash
# Pull from GHCR (public images — no authentication needed)
docker pull ghcr.io/stajichlab/novinvenio:0.5.0

# Convert to a Singularity/Apptainer SIF for HPC environments
singularity build novenio.sif docker://ghcr.io/stajichlab/novinvenio:0.5.0
```

### Run with Docker

All paths must be **absolute** so they resolve inside the container:

```bash
nextflow run stajichlab/NovInvenio \
    -profile docker \
    --config /absolute/path/to/analysis_config.csv \
    --data_dir /absolute/path/to/proteome_fastas \
    --pfam_hmm /absolute/path/to/Pfam-A.hmm \
    --swissprot_dmnd /absolute/path/to/uniprot_sprot.fasta.dmnd \
    --modelorgs_config /absolute/path/to/modelorgs.yaml
```

Nextflow automatically bind-mounts the launch directory and project directory,
so `--config`, `--data_dir`, and the annotation DB paths all resolve correctly
inside the container as long as they are absolute host paths.

> **UCR HPCC note:** compute nodes do not run a Docker daemon (no root access
> for regular users). Use the Singularity/Apptainer profile below for any
> real SLURM submission — the `docker` profile is only useful on a machine
> where you have Docker installed (e.g. a laptop or a VM).

### Run with Singularity

```bash
nextflow run stajichlab/NovInvenio \
    -profile singularity \
    --config /absolute/path/to/analysis_config.csv \
    --data_dir /absolute/path/to/proteome_fastas \
    --pfam_hmm /absolute/path/to/Pfam-A.hmm \
    --swissprot_dmnd /absolute/path/to/uniprot_sprot.fasta.dmnd
```

### UCR HPCC: Docker + SLURM combined

```bash
nextflow run stajichlab/NovInvenio \
    -profile slurm,docker \
    -c conf/ucr_hpcc_slurm.config \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

### UCR HPCC: Singularity + SLURM combined (recommended)

This is the recommended way to run NovInvenio on the HPCC — no pixi
environment needed, and no Docker daemon required on compute nodes.

```bash
module load nextflow singularity

nextflow run stajichlab/NovInvenio \
    -profile slurm,singularity \
    -c conf/ucr_hpcc_slurm.config \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

`run.sh` sets `NXF_SINGULARITY_CACHEDIR` (if not already set in your
environment) so the SIF is pulled/converted once instead of once per SLURM
array task. `conf/ucr_hpcc_slurm.config` adds `--bind /bigdata` to
`runOptions` — required because `--pfam_hmm`,
`--swissprot_dmnd`, `--modelorgs_config`, and `--data_dir` are passed to
processes as plain string params rather than staged Nextflow `path` inputs
(see `workflows/annotate.nf`), so Singularity's `autoMounts` can't detect
and bind them automatically. If your data lives outside `/bigdata`, add that
path to `runOptions` too.

### UCR HPCC: Docker + SLURM combined

```bash
nextflow run stajichlab/NovInvenio \
    -profile slurm,docker \
    -c conf/ucr_hpcc_slurm.config \
    --config configs/pezizo4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --pfam_hmm db/pfam/38.2/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

`conf/ucr_hpcc_slurm.config` provides the SLURM queue routing, AVX2 node
constraints for `famsa`, and preempt-queue settings specific to the UCR HPCC.
See that file's header for full usage documentation.

### Override the container image tag

```bash
# Use a local build instead of the registry
nextflow run stajichlab/NovInvenio -profile docker \
    --container_version 0.5.0 \
    --config /path/to/config.csv ...
```

### What the container does NOT include

The image ships **only the tool layer**.  User data files are always provided at
runtime via bind-mounts:

| Runtime artefact | How provided |
|---|---|
| Proteome & genome FASTAs | `--data_dir` (bind-mounted by Nextflow) |
| Config CSV | `--config` (bind-mounted by Nextflow) |
| Pfam-A HMM | `--pfam_hmm` (absolute path passed as a param) |
| SwissProt diamond DB | `--swissprot_dmnd` (absolute path passed as a param) |
| Model organism YAML | `--modelorgs_config` (absolute path passed as a param) |

The pipeline source (`bin/`, `lib/`, `modules/`, `workflows/`) is supplied by the
cloned repository and staged into each task's working directory by Nextflow — no
need to bake it into the image.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--config` | *(required)* | Analysis config CSV (see format below) |
| `--data_dir` | *(required)* | Directory containing FASTA files listed in the CSV |
| `--cluster_tool` | `pairwise` | Presence-matrix producer: `pairwise` (exact N² search) or `mmseqs` (scalable family-profile). See [Two analysis approaches](#two-analysis-approaches-pairwise-vs-family-profile-clustering) |
| `--run_tool` | `phmmer` | Search tool: `phmmer`, `diamond`, or `blast` (pairwise search + self-search paralog calibration) |
| `--family_min_seq_id` | `0.3` | mmseqs family-clustering identity threshold (`mmseqs` pathway) |
| `--family_cov` | `0.8` | mmseqs family-clustering coverage (`mmseqs` pathway) |
| `--hmm_presence_evalue` | `1e-3` | Family-HMM full-sequence E-value ceiling for presence (`mmseqs` pathway) |
| `--hmm_presence_cov` | `0.5` | Family-HMM minimum profile coverage for presence (`mmseqs` pathway) |
| `--family_chunk_size` | `200` | Families per parallel profile-build task (`mmseqs` pathway) |
| `--evalue` | `1e-5` | Fallback e-value for proteins with no detectable within-proteome paralog |
| `--parse_evalue` | `0.01` | Noise ceiling applied when parsing raw pairwise hits; the authoritative per-gene cutoff comes from the self-vs-self paralog search |
| `--ingroup_min_frac` | `0.75` | Min fraction of ingroup proteomes that must contain a hit |
| `--outgroup_min_frac` | `0.75` | Min fraction of outgroup proteomes that must contain a hit, for the loss-search direction (`loss_presence_matrix.tsv`) |
| `--loss_ingroup_max_frac` | `0.0` | Max fraction of the ingroup a loss candidate may still be present in. `0.0` = strictly absent from the ingroup; raise it (e.g. `0.1`) to allow genes *nearly*, but not entirely, lost — retained in a few ingroup species |
| `--core_min_frac` | `0.95` | Min presence fraction across *all* proteomes (ingroup + outgroup) for `core.html` |
| `--pfam_hmm` | `null` | Path to Pfam-A.hmm; skips Pfam annotation if unset |
| `--swissprot_dmnd` | `null` | Path to SwissProt `.dmnd` database; skips if unset |
| `--modelorgs_config` | `null` | YAML listing model organisms for gene name lookup (see `configs/modelorgs.yaml`) |
| `--report_sequences` | `novelties` | Which proteins carry a sequence in `novelties.html`: `novelties`, `all`, or `none`. Sequences dominate the file size |
| `--pdf_report` | `true` | Write `view/<project>/summary.pdf` (matplotlib figures). Set `false` to skip the step (gated via the process `when:` directive) |
| `--project` | *(auto)* | Output subdirectory name; defaults to config CSV basename |
| `--outdir` | `results` | Root output directory |
| `--hmm_mpi` | `false` | Run hmmsearch in MPI mode (requires MPI-enabled HMMER) |
| `--hmm_mpi_tasks` | `null` | Number of MPI tasks; defaults to `max_cpus` when `--hmm_mpi true` |
| `--container_version` | `0.5.0` | Docker/Singularity image tag (e.g. `0.5.0`, `latest`). See [Running with Docker / Singularity](#running-with-docker--singularity) |

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

### Generating config files (`bin/make_config.py`)

For larger sample collections, keep one master `configs/samples.csv` (every
available proteome listed once) instead of copy-pasting rows into a new config
per run:

```csv
Species,Strain,Protein,DNA,Short,Lineage
Neurospora crassa,OR74A,Ncra.pep.fa,Ncra.dna.fa,Ncra,Fungi;Ascomycota;Pezizomycotina;Sordariomycetes;Neurospora
```

`Lineage` is a semicolon-delimited list, highest→lowest rank. `bin/make_config.py`
slices a run-specific `GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup` config
out of that pool by matching lineage tokens — it's taxonomy-string matching, not
phylogenetic inference (no species tree is consulted):

```bash
bin/make_config.py \
    --samples configs/samples.csv \
    --ingroup-taxon Pezizomycotina \
    --max-per-outgroup-taxon 2 \
    --output configs/pezio_new.csv
```

This puts every species whose `Lineage` contains `Pezizomycotina` in the
ingroup. With no `--outgroup-taxon` given, the outgroup defaults to "siblings
of the ingroup's parent clade" (everything else under Ascomycota here), capped
at 2 species per sibling lineage so no single outgroup clade dominates.

Other options: `--outgroup-taxon NAME` (repeatable) to name the outgroup pool
explicitly instead of using the sibling default; `--exclude-taxon NAME`
(repeatable) to drop known-bad species; `--ingroup-short SHORT` (repeatable) to
pin exact species into the ingroup for focal-species configs (like
`configs/nirr.csv`, where the ingroup is intentionally one species); `--random`
+ `--seed` for random rather than deterministic stratified sampling; and
`--data-dir PATH` to fail fast if any selected `Protein`/`DNA` file doesn't
actually resolve on disk. If a real species tree becomes available later, a
`--tree` mode for true MRCA/sister-clade selection would be the natural
extension, but isn't implemented yet.

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
| `candidates.pfam.tblout` | Pfam hmmscan tblout (if `--pfam_hmm` set) |
| `candidates.swissprot.tsv` | SwissProt diamond hits (if `--swissprot_dmnd` set) |
| `presence_matrix.function.tsv` | Annotated matrix (gene_name, Best_Swissprot, Pfam_Names) |
| `novelties.<SHORT>.tsv` | Per-species novelty candidates (one file per ingroup species) |
| `loss_presence_matrix.tsv` | Mirror of `presence_matrix.tsv` with outgroup proteomes as the query — candidate lineage-specific losses |
| `loss_candidates.txt` / `loss_candidates.fa` | Outgroup protein IDs / FASTA present in ≥ `outgroup_min_frac` of the outgroup and in ≤ `loss_ingroup_max_frac` of the ingroup |
| `loss_tblastn_summary.tsv` | Loss candidates × ingroup genome TBLASTN hit matrix (the loss-direction validation check) |
| `loss_presence_matrix.function.tsv` | Annotated loss matrix (gene_name, Best_Swissprot, Pfam_Names, from the outgroup side) |
| `novelties.html` | Interactive report — browse the candidates and the presence matrix ([below](#interactive-report-noveltieshtml)) |
| `core.html` | Interactive report — near-universally conserved genes (`--core_min_frac`, default 0.95), for contrast against the novelty candidates ([below](#core-genes-report-corehtml)) |
| `losses.html` | Interactive report — candidate lineage-specific gene losses ([below](#losses-report-losseshtml)) |

A final `COLLATE_REPORTS` step also copies the three reports into a separate, shareable
folder alongside a landing page:

| File | Description |
|---|---|
| `view/<project>/report.html` | Landing page linking to the three reports, with a run summary (ingroup/outgroup proteomes, search tool, thresholds) |
| `view/<project>/novelties.html` · `core.html` · `losses.html` | Copies of the three reports, so the whole result set is one self-contained folder to `scp` or share |

### novelties.\<SHORT\>.tsv columns

Each file contains proteins that are:
- Present in species `<SHORT>`
- Present in ≥ `ingroup_min_frac` of all ingroup species
- Absent from all outgroup proteomes (protein-level search)
- Absent from all outgroup genomes (TBLASTN validation)

Columns include: `protein_id`, `source_proteome`, per-species presence columns,
`gene_name`, `product_description`, `function_source`, `Best_Swissprot`,
`Pfam_Names`.

### Single-ingroup-species configs

With only one `IN` species in the config CSV, `--ingroup_min_frac` is trivially
satisfied (1/1), so the only filter that determines candidates is "absent from every
outgroup proteome." It's normal for this to yield zero candidates — a single species
compared against a handful of outgroups may simply have no protein that fails to hit
any of them. When `candidates.txt`/`candidates.fa` come out empty, `MMSEQS_CLUSTER`
(`modules/mmseqs_cluster.nf`) detects the empty FASTA and skips clustering rather than
invoking mmseqs, which otherwise dies ungracefully (`Cannot seek at the end of file`)
on empty input — it emits empty `clusters_rep_seq.fasta`/`clusters_all_seqs.fasta`/
`clusters_cluster.tsv` and the pipeline continues to completion with empty downstream
outputs.

## Interactive report (novelties.html)

The pipeline writes `results/<project>/novelties.html` — a single file for browsing the
novelty candidates and the presence/absence matrix. It has no external dependencies and
makes no network requests, so copy it anywhere and open it directly:

```bash
# from your laptop
scp cluster:/path/to/NovInvenio/results/pezizo4_asco/novelties.html .
open novelties.html       # macOS  (Linux: xdg-open, or just drag into a browser)
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
    --matrix results/pezizo4_asco/presence_matrix.function.tsv \
    --config configs/pezizo4_asco.csv \
    --tblastn_summary results/pezizo4_asco/tblastn_summary.tsv \
    --novelties results/pezizo4_asco/novelties.*.tsv \
    --candidates_fa results/pezizo4_asco/candidates.fa \
    --output results/pezizo4_asco/novelties.html
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

## Core genes report (core.html)

The pipeline also writes `results/<project>/core.html` — genes present in ≥
`--core_min_frac` (default 0.95) of *every* sampled proteome, ingroup and outgroup alike.
It answers the opposite question from `novelties.html`: not "unique to this lineage" but
"conserved across the whole sampled tree" — useful as a sanity check (a run with almost
no core genes probably has a search-sensitivity or config problem) and as contrast for
the novelty candidates.

It needs no new search or annotation step — it re-reads `presence_matrix.function.tsv`
(or the plain matrix, if annotation was skipped) and `clusters_cluster.tsv`, the same
inputs `novelties.html` uses minus the TBLASTN summary and novelties tables. Sequences
are never embedded (there's no BLASTP-at-NCBI link); instead the detail panel links out
to InterPro (Pfam), UniProt/AlphaFold (SwissProt), FungiDB (model-organism genes), or a
generic NCBI Protein search as a fallback for a core gene with no other hit — worth a
manual spot-check.

Regenerate standalone the same way as `novelties.html`:

```bash
bin/make_core_report.py \
    --matrix results/pezizo4_asco/presence_matrix.function.tsv \
    --config configs/pezizo4_asco.csv \
    --cluster_tsv results/pezizo4_asco/clusters/clusters_cluster.tsv \
    --output results/pezizo4_asco/core.html
```

## Losses report (losses.html)

The pipeline also writes `results/<project>/losses.html` — **gene families conserved
across the outgroups but (nearly) absent from the ingroup**: candidate lineage-specific
gene losses, the kind of signal worth chasing for a missing pathway. It answers a question
`novelties.html` cannot: that report only ever runs searches with ingroup proteomes as the
query, so it has no way to ask "is this outgroup gene present in the ingroup?" `losses.html`
is built from a second, symmetric search direction that queries with the outgroup instead —
the same `--cluster_tool` choice applies: `pairwise` runs `LOSS_SEARCH`
(`workflows/loss_search.nf`), roughly doubling total pairwise search volume; `mmseqs` seeds
gene families from the **outgroup** and mirrors the family-profile pathway (its outputs land
under `loss_families/`, `loss_presence_matrix.tsv`, etc.).

**Which rows appear.** A protein is a loss candidate when it is present in ≥
`--outgroup_min_frac` of the outgroup (conserved there) **and** in ≤
`--loss_ingroup_max_frac` of the ingroup (lost, or nearly so). With the default
`--loss_ingroup_max_frac 0.0` that second condition means "absent from *every* ingroup
proteome"; raise it to surface genes retained in just one or two ingroup species. The
report enforces this predicate itself (the annotated `loss_presence_matrix.function.tsv`
is the full presence table, not a pre-filtered list), so pass the same threshold values
the search used when regenerating it standalone.

**The family is the unit.** A single loss is much more convincing when the gene is
recovered independently across many outgroup species than when it shows up in only one
(which could be that species' own gain, an assembly artefact, or a paralog quirk). Rows
are grouped into gene families (mmseqs clusters) and each row carries two family-level
columns:

- **Outgroup breadth** — how many outgroup species carry *any* member of the family.
- **Ingroup retained** — how many ingroup species still keep a member (`0` is a clean
  loss across the whole ingroup; non-zero only appears when `--loss_ingroup_max_frac > 0`).

The default sort ranks the cleanest, broadest losses first: fewest ingroup species
retaining the gene, then widest outgroup breadth, then no ingroup genomic evidence.
`loss_tblastn_summary.tsv` — TBLASTN of the loss candidates against ingroup *genomic* DNA,
not just the annotated proteins — flags (but does not remove) rows where a hit suggests the
"loss" might just be a missed gene model rather than a real one. Sequences are never
embedded; the detail panel's external links resolve against the *outgroup* protein (e.g.
"this looks like S. cerevisiae's `ERG3`"), since that's where the gene actually is.

Regenerate standalone (pass the same thresholds the run used):

```bash
bin/make_losses_report.py \
    --matrix results/pezizo4_asco/loss_presence_matrix.function.tsv \
    --config configs/pezizo4_asco.csv \
    --tblastn_summary results/pezizo4_asco/loss_tblastn_summary.tsv \
    --cluster_tsv results/pezizo4_asco/clusters/loss_clusters_cluster.tsv \
    --outgroup_min_frac 0.75 \
    --loss_ingroup_max_frac 0.0 \
    --output results/pezizo4_asco/losses.html
```

## Running on SLURM

The SLURM profile works together with the site configuration file
`conf/ucr_hpcc_slurm.config`, which provides queue routing, AVX2 node
constraints, and resource settings specific to the UCR HPCC.

### Pixi environment (no container)

```bash
nextflow run stajichlab/NovInvenio \
    -profile slurm -c conf/ucr_hpcc_slurm.config \
    --config configs/pezizo4_asco.csv \
    --data_dir /bigdata/stajichlab/shared/data/fungi/proteomes \
    --run_tool diamond \
    --ingroup_min_frac 0.8 \
    --pfam_hmm /absolute/path/to/Pfam-A.hmm \
    --swissprot_dmnd /absolute/path/to/uniprot_sprot.fasta.dmnd \
    --modelorgs_config /absolute/path/to/modelorgs.yaml
```

### Docker container on SLURM

```bash
nextflow run stajichlab/NovInvenio \
    -profile slurm,docker -c conf/ucr_hpcc_slurm.config \
    --config /absolute/path/to/config.csv \
    --data_dir /absolute/path/to/proteomes \
    --run_tool diamond \
    --ingroup_min_frac 0.8 \
    --pfam_hmm /absolute/path/to/Pfam-A.hmm \
    --swissprot_dmnd /absolute/path/to/uniprot_sprot.fasta.dmnd
```

### MPI mode for hmmsearch

When running the downstream `HMMSEARCH` process on SLURM with an MPI-enabled
HMMER build (e.g. from conda-forge), pass `--hmm_mpi true` to use
`mpirun -np N hmmsearch --mpi` instead of `hmmsearch --cpu N`.  The SLURM
profile automatically requests `--ntasks N` instead of `--cpus-per-task N`.

```bash
nextflow run stajichlab/NovInvenio \
    -profile slurm,docker -c conf/ucr_hpcc_slurm.config \
    --config configs/pezizo4_asco.csv \
    --data_dir /absolute/path/to/proteomes \
    --run_tool diamond \
    --hmm_mpi true \
    --hmm_mpi_tasks 16
```

If `--hmm_mpi_tasks` is omitted it defaults to `--max_cpus` (32 by default).

### Database paths

The pipeline resolves `--pfam_hmm`, `--swissprot_dmnd`, and
`--modelorgs_config` to absolute paths at startup and creates a `db/` symlink
inside every task work directory.  You can pass relative paths
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
