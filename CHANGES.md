# Changelog

## Unreleased

### New workflows

- **`workflows/annotate.nf`** — `ANNOTATE_MATRIX` process runs Pfam hmmscan
  and/or SwissProt diamond blastp against candidate proteins, then calls
  `annotate_presence_matrix.py` to merge results into a functional presence
  matrix (`presence_matrix.function.tsv`).  Both annotation steps are optional:
  if `--pfam_hmm` or `--swissprot_dmnd` is not set, the corresponding step is
  skipped with a `touch` placeholder.

- **`workflows/validate.nf`** extended — `SUMMARIZE_TBLASTN` process added.
  After per-outgroup TBLASTN runs, results are merged into a protein × genome
  hit matrix (`tblastn_summary.tsv`) via `summarize_tblastn.py`.  The `VALIDATE`
  workflow now emits both raw hits and the summary matrix.

- **`workflows/summarize.nf`** — `MAKE_NOVELTIES` process calls
  `make_novelties.py` to produce per-ingroup-species novelty tables
  (`novelties.<SHORT>.tsv`), filtering to proteins that are ingroup-present,
  outgroup-absent at both proteome and genome levels.

### New scripts

- **`bin/summarize_tblastn.py`** — reads per-outgroup TBLASTN TSVs, expands
  representative protein hits to cluster members via the cluster TSV, and
  writes a protein × outgroup-genome presence/absence matrix.

- **`bin/make_novelties.py`** — reads the annotated presence matrix, the
  TBLASTN summary, and the analysis config to produce per-species novelty
  candidate tables with functional annotation columns.

- **`bin/parse_self_hits.py`** — extracts self-hit bitscores from a proteome's
  self-search result for use as a normalisation reference.

- **`bin/split_fasta.py`** — splits a multi-FASTA by sequence ID (used for
  per-cluster FASTA preparation).

- **`bin/filter_candidates.py`** — post-processes the presence matrix to emit a
  candidates list respecting `--ingroup_min_frac`.

- **`bin/annotate_presence_matrix.py`** rewritten — replaced hardcoded Ncra/Afum
  arguments with `--modelorgs_config` YAML + `--launch_dir`.  Now also accepts
  `--pfam_hits` and `--swissprot_hits` to add `Pfam_Names` and `Best_Swissprot`
  columns to the output matrix.

### Generalised model organism annotation

- **`lib/model_organisms.py`** (new) — `ModelOrgConfig` dataclass and
  `ModelOrgAnnotator` class.  Supports three `id_transform` strategies:
  - `direct` — protein ID used as gene lookup key unchanged
  - `strip` — regex `strip_pattern` removed from protein ID
  - `diamond_fasta` — protein ID → diamond hit table → reference FASTA header
    gene field → gene name CSV
  Relative paths in the YAML are resolved against `--launch_dir` (the Nextflow
  project directory), not the YAML file's location.

- **`configs/modelorgs.yaml`** (new) — example config for Ncra (diamond_fasta
  transform) and Afum (strip transform) using FungiDB gene name tables.

- **`pipeline/02_annotate.sh`** updated — replaced `--ncra_*` / `--afum_*` flags
  with a single `--modelorgs_config` argument.

### MPI support for hmmsearch

- **`modules/hmmsearch.nf`** — new params `--hmm_mpi` and `--hmm_mpi_tasks`
  switch the hmmsearch invocation between `mpirun -np N hmmsearch --mpi` (MPI
  mode) and `hmmsearch --cpu N` (threaded mode).

- **`nextflow.config`** — new `hmmsearch` process label with separate SLURM
  resource rules: MPI mode requests `--ntasks N` via `clusterOptions`; threaded
  mode requests `--cpus-per-task 4` in the normal way.

### Bug fixes and infrastructure

- **`lib/Helpers.groovy`** (new) — Groovy class auto-loaded by Nextflow.
  `Helpers.projectName(params)` derives the output subdirectory from
  `params.project` if set, otherwise from the config CSV basename.  All
  `publishDir` closures now use this helper, eliminating null-value warnings
  when `--project` is not supplied.

- **Absolute DB path resolution** — `main.nf` converts `--pfam_hmm` and
  `--swissprot_dmnd` to absolute paths at startup.  `--modelorgs_config` is
  resolved to an absolute path and passed explicitly through the workflow
  channel to `ANNOTATE_MATRIX`, avoiding DSL2 params-scoping issues.

- **Work directory symlinks** — `beforeScript` in `nextflow.config` creates
  `db -> ${projectDir}/db` and `configs -> ${projectDir}/configs` symlinks
  inside each task's work directory.  This lets relative paths like
  `db/pfam/Pfam-A.hmm` and `configs/modelorgs.yaml` resolve correctly from
  any task without requiring absolute paths.

- **`high_cpu` default CPU cap** — global default for `high_cpu` label lowered
  from `params.max_cpus` (32) to 8 to avoid exhausting local machine CPUs.
  The SLURM profile still overrides to `params.max_cpus`.

- **`pixi.toml`** — added `pyyaml >= 6.0` dependency (required by
  `lib/model_organisms.py`).

- **`workflows/cluster.nf`** — added `candidates_fa` to the `emit:` block so
  the downstream `ANNOTATE` workflow can receive the candidate FASTA.
