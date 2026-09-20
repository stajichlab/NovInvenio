# Changelog

## Unreleased

### New: island synteny view

- **`lib/pfam_classes.py`**, **`lib/island_synteny.py`**, **`lib/island_synteny_template.py`**,
  **`bin/pangenome_island_synteny.py`**, **`modules/pangenome/island_synteny.nf`** (issue #116)
  — a self-contained `island_synteny.html` page: for each selected accessory island, a
  presence/absence grid with rows = strains collapsed into distinct haplotypes and columns =
  member families in locus order, so a deletion breakpoint shows up as a vertical edge. Runs
  inside the existing `--pangenome_island_pfam_hmm` conditional block (it consumes
  `report_tables/islands_with_domains.tsv`, which only exists there), and degrades to a valid
  "no islands" page rather than failing when a study has none.
- Two new params: `--pangenome_viz_top_islands` (default `50`) caps how many islands are
  embedded in the page, and reuses the existing `--pangenome_top_islands_min_strains`
  (default `2`) filter that the "Top islands (by size)" report table already applied.
- **Single-strain islands are excluded from the view.** On the real 529-strain
  genus_vs_ureesii study, 62% of located islands (17,218 of 27,836) are present in a single
  strain -- a single-strain island's grid is one filled row with no breakpoint to show, so
  it carries no synteny information this view exists to display.
- **The payload is valid only for the run that produced it.** A family's ID is its mmseqs
  representative, and gene-family IDs are not stable across pipeline re-runs -- 73% ID
  overlap was measured between two runs of the same study. Do not compare or merge
  `island_synteny.html` payloads across separate pipeline invocations.
- **Known scaling limit (not fixed here):** `PresenceMatrix.from_tsv()` materialises the
  whole presence matrix regardless of how few islands are actually drawn -- O(families x
  strains), not O(islands). Measured on the real 529-strain genus_vs_ureesii study (54,421
  families x 530 strains, 10,602 qualifying islands, 50 drawn): peak RSS 5.54 GB, 48s wall.
  `ISLAND_SYNTENY`'s memory is sized from this measurement (12 GB, see below); reducing the
  actual memory footprint needs a real design change (e.g. loading only the rows/columns the
  selected islands touch) and is deferred to a follow-up.

### Fixed: ISLAND_SYNTENY memory sized from a real measurement, truncation now reported

- `modules/pangenome/island_synteny.nf`'s `low_cpu` label default (4 GB) is not enough --
  the process was measured at 5.54 GB peak RSS on the real 529-strain genus_vs_ureesii study
  (see above), which would OOM-kill on exactly the dataset this feature exists to visualise.
  Added a `withName: '.*ISLAND_SYNTENY'` override in `nextflow.config` (matching this repo's
  existing convention for per-process resource overrides) giving it 12 GB, ~2.2x headroom
  over the measurement, while keeping `cpus = 1` from the `low_cpu` label.
- `bin/pangenome_island_synteny.py`'s stderr summary previously reported only "N islands
  drawn, M excluded (< min_strains)", which under-reports: on the same real run, 10,602
  islands qualified but only 50 were drawn, and the summary never said where the other
  ~10,552 went. The summary now also reports `n_islands_truncated` (islands that qualified
  but were cut by `--top_islands`), covered by a new test in
  `tests/test_pangenome_island_synteny.py`.

### Fixed: `lib/compressed_io.py` now fails loudly on a missing or corrupt compressed input

`open_maybe_compressed()` previously returned an **empty stream and exit 0** for a missing
file (`FileNotFoundError`) or a non-zstd/corrupt input, rather than raising. That is a
silent-wrong-answer risk for any of its eight consumers merged in PR #113
(`pangenome_build_islands.py`, `pangenome_report_tables.py`,
`pangenome_build_family_positions.py`, `pangenome_island_synteny.py`, and others) -- a
downstream script would proceed as if the compressed input legitimately had zero rows
instead of failing. It now raises `FileNotFoundError` for a missing input and `RuntimeError`
for a non-zstd/corrupt one. PR #113's `zstd -dc -f` flag is replaced by resolving the
symlink in Python before invoking `zstd`, because `-f` itself was the root cause: it passes
non-zstd bytes straight through with exit 0 (needed only to read through Nextflow's staged
symlinks) -- which is exactly the silent-pass-through this fix closes.

### Changed: the four largest pangenome TSV intermediates are now zstd-compressed

`lib/compressed_io.py` has shipped `open_maybe_compressed_write()` since it was
added, but nothing ever called it. Apart from `RESCUE_PASS`'s tblastn (which
compresses at the shell level in `modules/pangenome/rescue.nf`), every
pangenome output published as plain text -- ~4.1 GB per 529-strain run across
four files, on shared storage. These now publish as `.zst` (issue #111):

| Output | Was (genus_vs_ureesii, 529 strains) |
|---|---|
| `pair_classification.tsv.zst` | 1.99 GB plain |
| `cooccurring_pairs.tsv.zst` | 1.69 GB plain |
| `family_positions.tsv.zst` | 0.24 GB plain |
| `gene_positions.tsv.zst` | 0.22 GB plain |

Producers (`pangenome_cooccurrence.py`, `pangenome_pair_classification.py`,
`pangenome_build_family_positions.py`, `pangenome_build_gene_positions.py`)
write through `open_maybe_compressed_write()`; consumers
(`pangenome_build_islands.py`, `pangenome_report_tables.py`,
`pangenome_build_family_positions.py`) read through `open_maybe_compressed()`.
Compression is selected by the output filename's suffix -- there is no new
parameter, because a flag would duplicate a control the filename already
carries. Plain and `.gz` inputs still work unchanged, so an existing
uncompressed results directory stays readable by these scripts.

`presence_matrix.tsv`/`presence_matrix.rescued.tsv` (0.38 GB) are deliberately
left plain for now: widest consumer fan-out, smallest share of the win. See
issue #111.

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
  this, was closed same-day as not-a-bug.
- **2026-09-17: real-data concordance benchmark, complete.** Re-ran both
  backends via real SLURM submission (`-profile slurm -c
  conf/ucr_hpcc_slurm.config`) against the identical 5-strain subset: 37/37
  processes, 0 failures, for both. diamond: 10,468 families. mmseqs:
  11,830 families. Comparing the two `tier1_cluster.tsv` outputs over the
  same 43,239 proteins: **ARI = 0.9403**, pair recall = 0.8958, pair
  precision = 0.9895 -- substantially stronger concordance than the
  2026-09-08 novelty/loss-pathway benchmark (ARI 0.79, 72%/87%), since this
  uses the pangenome subworkflow's own tight, matched 90%-identity tier-1
  thresholds for both tools. Full writeup in `.living/decisions.md`'s
  2026-09-17 entry; `todo/diamond-tier1-cluster-backend.md` closed.
  `--pangenome_cluster_backend diamond` is now a fully validated
  alternative to mmseqs for this subworkflow.

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
