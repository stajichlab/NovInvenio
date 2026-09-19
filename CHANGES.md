# Changelog

## Unreleased

### New: Leiden trans-module detection wired into the pangenome subworkflow

- **`modules/pangenome/trans_modules.nf`** (`LEIDEN_MODULES`, `MODULE_DOMAINS`),
  **`bin/pangenome_detect_trans_modules.py`**, **`bin/pangenome_module_domains.py`**
  — promotes the `coccidioides_pangenome`/`Afumigatus_pangenome` studies'
  manual, study-specific Leiden module-detection scripts into the core
  pipeline as a real Nextflow step. `LEIDEN_MODULES` runs unconditionally
  right after `PAIR_CLASSIFICATION` (cheap; no gating param of its own) and
  collapses `trans`-classified (physically unlinked, statistically
  significant) co-occurrence pairs into gene-family modules via Leiden
  community detection (Traag, Waltman & van Eck 2019). `MODULE_DOMAINS` runs
  inside the existing `--pangenome_island_pfam_hmm` block, cross-referencing
  module assignments against the same Pfam scan `FAMILY_PFAM_SCAN` already
  produces for island enrichment -- one row per module, its member-family
  count, and the Pfam domain names found among them. Deliberately does NOT
  hardcode any specific domain search (e.g. NACHT/HET, PKS/NRPS) -- that's a
  one-line `grep`/`awk` over `module_domains.tsv`, not a pipeline parameter;
  see the `genus_vs_ureesii` real-run analysis this generalizes from
  (`TRANS_MODULE_REPORT.md` in that study) for a worked example.
- **Graceful degradation, not a hard failure, for zero `trans` edges.** The
  original ported script (`sys.exit(1)` on an empty edge list) would crash
  the whole pipeline for any single-species-ingroup study --
  `pangenome_pair_classification.py`'s `min_clades>=2` requirement for a
  `trans` call means a study with only one species commonly produces zero
  `trans` pairs (confirmed for real: both `coccidioides_pangenome`
  species-reciprocal runs, `immitis_in_posadasii_out`/`posadasii_in_immitis_out`,
  had zero). Both new scripts write valid, header-only output and exit 0
  instead. Verified with a real end-to-end `-profile local` smoke test
  against a 5-strain single-species subset (zero `trans` edges, both steps
  degraded gracefully, `[SUCCESS] completed=45 failed=0`) as well as unit
  tests covering the real-clustering path.
- New params: `--pangenome_leiden_resolution` (default 1.0 -- **not**
  validated for every study; a real resolution-stability sweep on the
  `genus_vs_ureesii` run found 1.0 too coarse there, first informative
  structure only appeared at 2.0), `--pangenome_leiden_seed`,
  `--pangenome_module_min_size`.
- New pixi dependencies: `python-igraph`, `leidenalg`.

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
  `protein_id=` -- fix proposed, in progress on a separate track).

### diamond as a tier-1 clustering backend (see docs/adr/0003)

- `--pangenome_cluster_backend diamond` is no longer hard-blocked in
  `pangenome.nf`. `CLUSTER_TIER1` (`modules/pangenome/prefix_and_cluster.nf`)
  now runs the new `bin/verify_diamond_cluster_ids.py` against diamond's raw
  `tier1_cluster.tsv` before anything downstream reads it, failing loud on
  any id that doesn't match a real FASTA header verbatim -- the diamond
  equivalent of what `bin/restore_mmseqs_cluster_ids.py` does for mmseqs
  (diamond is expected to need verification only, not restoration, since it
  doesn't collapse headers the way mmseqs does). Unit-tested against
  adversarial headers, including the previously-untested case of a
  Short-prefixed, multi-pipe UniProt-style id (`Afum|sp|ACC|NAME`).
- `bin/pangenome_cooccurrence.py` and every other downstream consumer needed
  no changes -- they already read `tier1_cluster.tsv`/`presence_matrix.tsv`
  structurally, independent of which clustering tool produced them.
- **2026-09-17: validated end-to-end on real data.** Ran `pangenome.nf
  --pangenome_cluster_backend diamond -profile local` against a real
  5-strain *Coccidioides immitis* subset: 37/37 processes, 0 failures,
  10,467 real gene families with a normal core/shell/singleton frequency
  distribution, `verify_diamond_cluster_ids.py` passing inside the actual
  run. A same-subset mmseqs comparison run for a direct concordance check
  hit `mmseqs easy-cluster` SIGILLing on the test node (no AVX2) instead of
  producing a cluster count -- turned out to be a testing-methodology
  artifact, not a pipeline gap: the `-C ryzen|broadwell|cascade` AVX2
  node-pin for `CLUSTER_TIER1` already exists (line 18 above), it's just
  only applied under `-profile slurm -c conf/ucr_hpcc_slurm.config`, and
  this comparison run used `-profile local` alone. Issue #101, filed over
  this, was closed same-day as not-a-bug. The real-data cluster-quality
  concordance benchmark remains open, unblocked -- it just needs
  `-profile slurm` instead of `-profile local` (tracked in
  `todo/diamond-tier1-cluster-backend.md`).

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
