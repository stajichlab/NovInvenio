# Changelog

## Unreleased

### New: pangenome cluster-profiling subworkflow (branch `pangenome-profiling-module`)

- **`pangenome.nf`, `workflows/pangenome_profile.nf`, `modules/pangenome/*.nf`,
  `bin/pangenome_*.py`, `lib/pangenome_*.py`, `lib/compressed_io.py`** — a new,
  reusable DSL2 subworkflow for pangenome cluster profiling (tier-1
  clustering, presence/absence, genome-level tblastn rescue pass,
  strain dedup + clade assignment, frequency binning, co-occurrence, and
  physical-linkage pair classification), ported from a one-off study script
  chain (NovInvenio_Investigations' Afumigatus_pangenome study) into a
  generic module any future species-set study can invoke via its own entry
  point (`pangenome.nf`, separate from `main.nf`'s novelty/loss workflow).
  See that study's `NEXTFLOW_MIGRATION_NOTES.md` for the full design record.
  New `params.pangenome_*` config surface in `nextflow.config`; new
  `conf/ucr_hpcc_slurm.config` overrides for the tier-1 clustering (AVX2),
  per-strain tblastn, co-occurrence, and rescue-position-extraction
  processes.
- **Real end-to-end validation (2026-09-15, SLURM job 28428780, AVX2 node)**:
  the full core chain -- CLUSTER_TIER1 -> PRESENCE_MATRIX -> per-strain
  rescue (EXTRACT_ABSENT_QUERIES -> TBLASTN_PER_STRAIN -> RESCUE_PASS ->
  EXTRACT_RESCUE_POSITIONS) -> FREQUENCY_BINS -> COOCCURRENCE ->
  PAIR_CLASSIFICATION, plus MASH-based DEREPLICATE/ASSIGN_CLADES -- ran to
  completion against real data (8 *A. fumigatus* ingroup + 2 outgroup
  strains from the Afumigatus_pangenome study): 57/57 processes succeeded,
  0 failures. Real, non-trivial outputs (20,043 families binned
  5,815 core / 4,539 shell / 9,689 singleton). `cooccurring_pairs.tsv`/
  `pair_classification.tsv` came back empty -- verified via the
  COOCCURRENCE process log this is the statistically correct outcome (0 of
  1,011,439 prefilter-surviving pairs clear FDR 0.05 at n=8 representative
  strains), not a swallowed error. This was the one part of the module an
  independent review had flagged as unconfirmed (most recently added,
  least-tested code path); it is no longer unconfirmed. Still open: a
  real GFF3-parsing gap surfaced testing against a second species
  (Coccidioides funannotate-style GFF3s use CDS `Parent=`, not
  `protein_id=` -- fix proposed, in progress on a separate track) and the
  diamond-backend clustering path remains hard-disabled (unvalidated
  ID-matching, no `restore_mmseqs_cluster_ids.py` equivalent yet).

### Performance: batched DIAMOND_SEARCH (0.6.1)

- **`modules/diamond.nf`** — `DIAMOND_SEARCH` now runs one Nextflow task per
  *query genome* (looping over every target inside the task) instead of one
  task per (query, target) pair. For an N-ingroup × M-other-genome study this
  cuts job count from N×M to N — for an 11-genome study, ~110 individual
  diamond jobs become ~11. Real per-pair diamond compute is small (a few
  seconds even for the largest proteomes measured); job *count*, not diamond
  runtime, was the actual cost driver on a shared SLURM queue with a 50/min
  submission rate limit — the same "trade job count for job size" pattern
  already used for `FAMILY_HMMSEARCH` (see `hmm_search_chunk_size`).
  `workflows/search.nf` and `workflows/loss_search.nf` group `pairs_ch` by
  query genome (`groupTuple`) before calling `DIAMOND_SEARCH`, then
  reconstruct the `[meta_pair, file]` shape `PARSE_HITS` already expects from
  each output file's self-describing name — `PARSE_HITS` itself needed no
  changes. Output file names/format are unchanged.

- Switched `DIAMOND_SEARCH` from `storeDir` to normal resume-based task
  caching + `publishDir`: `storeDir`'s skip-if-exists check isn't well-defined
  against a variable-length glob output (the batched job's file count varies
  per query), so partial-cache correctness wasn't guaranteed under the new
  shape. This also means this step's real cost is now visible in Nextflow's
  own trace files — previously a `storeDir` cache hit skipped the task (and
  its trace row) entirely, so the true all-vs-all search cost was invisible
  to every trace file from any run after the first.

  Validated on a real, cache-cleared, containerized SLURM rerun of
  `pezizo_set1` (11 genomes): job count dropped from ~110 to 11
  (5 `SEARCH:DIAMOND_SEARCH` + 6 `LOSS_SEARCH:DIAMOND_SEARCH`), full
  SEARCH+LOSS_SEARCH stage completed in under 10 minutes wall-clock, and
  every output (`presence_matrix.tsv`, `candidates.txt`,
  `loss_presence_matrix.tsv`, `loss_candidates.txt`) was byte-identical to
  the previously-published pairwise results.

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
