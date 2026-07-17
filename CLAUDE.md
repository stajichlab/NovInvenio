# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NovInvenio identifies lineage-specific genes: proteins present in a defined ingroup (≥N% of members) but absent from all outgroup proteomes. It uses pairwise protein searches (phmmer/diamond/blast), self-vs-self paralog calibration, mmseqs2 clustering, TBLASTN validation against outgroup genomes, functional annotation (Pfam + SwissProt), and model-organism gene-name lookup to produce per-species novelty candidate tables.

## Commands

```bash
# Set up environment
pixi install

# Run the full pipeline
nextflow run main.nf \
    --config configs/pezio4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool phmmer

# Run tests
pytest tests/ -v

# Run a single test file
pytest tests/test_parse_hits.py -v

# Lint
pixi run lint
```

## Repository Layout

```
NovInvenio/
├── main.nf                        # Nextflow entry point; imports and calls all workflows
├── nextflow.config                # params{} defaults and executor/resource settings
├── pixi.toml                      # [workspace] table — tool dependencies
├── modules/
│   ├── phmmer.nf                  # PHMMER_SEARCH — phmmer pairwise search, storeDir
│   ├── diamond.nf                 # DIAMOND_SEARCH — diamond blastp pairwise, storeDir
│   ├── blast.nf                   # BLAST_SEARCH — blastp pairwise, storeDir
│   ├── self_search.nf             # PHMMER_SELF / DIAMOND_SELF / BLAST_SELF — self-vs-self, storeDir
│   ├── parse_hits.nf              # PARSE_HITS — normalises raw hits → common TSV
│   ├── parse_self_hits.nf         # PARSE_SELF_HITS — extracts per-gene paralog cutoffs
│   ├── build_presence_matrix.nf   # BUILD_PRESENCE_MATRIX — presence/absence matrix + candidates.txt
│   ├── mmseqs_cluster.nf          # MMSEQS_CLUSTER — mmseqs2 easy-cluster of candidates
│   ├── tblastn.nf                 # TBLASTN — translated search vs outgroup genomes
│   ├── hmmbuild.nf                # HMMBUILD (imported but not yet wired into cluster workflow)
│   └── hmmsearch.nf               # HMMSEARCH (imported but not yet wired into cluster workflow)
├── workflows/
│   ├── search.nf                  # SEARCH — pairwise search + self-hits + presence matrix
│   ├── cluster.nf                 # CLUSTER — extract candidates + mmseqs2 clustering
│   ├── validate.nf                # VALIDATE — TBLASTN + SUMMARIZE_TBLASTN
│   ├── annotate.nf                # ANNOTATE — Pfam hmmsearch + SwissProt diamond + matrix merge
│   ├── summarize.nf               # SUMMARIZE — MAKE_NOVELTIES per ingroup species
│   └── report.nf                  # REPORT — MAKE_REPORT, interactive report.html
├── bin/                           # Executable Python helpers; NF adds bin/ to PATH automatically
│   ├── parse_hits.py              # Normalise phmmer/diamond/blast → query_id, target_id, evalue, bitscore, query_proteome, target_proteome TSV
│   ├── parse_self_hits.py         # Self-hit rank-2 → per-gene paralog_cutoffs.tsv
│   ├── build_presence_matrix.py   # Paralog-aware matrix construction + candidates.txt
│   ├── extract_candidates.py      # Pull candidate sequences from ingroup FASTAs
│   ├── summarize_tblastn.py       # Aggregate per-genome TBLASTN TSVs → protein × genome matrix
│   ├── annotate_presence_matrix.py# Add gene_name / Pfam / SwissProt columns
│   ├── make_novelties.py          # Per-species novelties.<SHORT>.tsv (with --skip_tblastn_filter option)
│   ├── make_report.py             # Self-contained interactive report.html
│   ├── make_config.py             # Generate a run config CSV from configs/samples.csv by taxon/lineage matching
│   ├── filter_candidates.py       # (standalone helper — not yet wired into pipeline)
│   ├── split_fasta.py             # (standalone helper)
│   └── summarize_clusters.py      # (standalone helper)
├── lib/                           # Shared Python package (lib/__init__.py present)
│   ├── hits.py                    # PARSERS dict + open_input + Hit namedtuple — used by parse_hits.py and parse_self_hits.py
│   ├── config_parser.py           # parse_config() → list[Sample]; short_to_group()
│   ├── fasta.py                   # FASTA utilities
│   ├── model_organisms.py         # ModelOrgAnnotator — YAML-driven gene name lookup
│   ├── report_data.py             # build_payload() — assembles the report's embedded JSON
│   ├── report_template.py         # HTML_TEMPLATE — the report page (HTML/CSS/JS, no deps)
│   └── Helpers.groovy             # Helpers.projectName(params) — derives results subdirectory name
├── tests/
│   ├── conftest.py
│   ├── test_config_parser.py
│   ├── test_hits.py
│   └── data/                      # test.csv + small FASTAs for -profile test
├── configs/
│   ├── samples.csv                 # Master species pool (Species,Strain,Protein,DNA,Short,Lineage) — source for bin/make_config.py
│   ├── pezio4_asco.csv            # Main analysis config (Pezizomycotina ingroup, Asco outgroup)
│   └── modelorgs.yaml             # Model organism YAML for gene name lookups
└── results/
    └── <config_basename>/
        ├── search_cache/          # storeDir — pairwise + self-hit raw outputs, never re-run
        ├── self_hits/             # Per-species paralog_cutoffs.tsv files
        ├── clusters/              # mmseqs2 cluster FASTA and TSV
        ├── presence_matrix.tsv
        ├── presence_matrix.function.tsv
        ├── candidates.fa
        ├── candidates.txt
        ├── tblastn_summary.tsv
        ├── novelties.<SHORT>.tsv  # One per ingroup species
        └── report.html            # Self-contained interactive report
```

## Architecture: Data Flow

1. **Config CSV** is parsed in `main.nf` using `Channel.fromPath(...).splitCsv(header:true)`. FASTA paths are resolved via `resolve_fa()`, which searches flat layout and subdirs (`pep/`, `dna/`, `genome/`, `scaffolds/`). Three channels are emitted: `ingroup_prot_ch`, `outgroup_prot_ch`, `outgroup_dna_ch`.

2. **SEARCH workflow** (`workflows/search.nf`):
   - Cross-joins every ingroup proteome against all other proteomes (ingroup ∪ outgroup, excluding self-pairs).
   - Runs the selected tool (`PHMMER_SEARCH` / `DIAMOND_SEARCH` / `BLAST_SEARCH`), output cached in `search_cache/` via `storeDir`.
   - `PARSE_HITS` normalises each raw result into a TSV with columns `query_id, target_id, evalue, bitscore, query_proteome, target_proteome`.
   - Also runs self-vs-self searches (`PHMMER_SELF` / `DIAMOND_SELF` / `BLAST_SELF`, `-E 100 --max-target-seqs 2`) to capture the rank-2 (best paralog) hit per protein. Outputs also `storeDir`-cached.
   - `PARSE_SELF_HITS` → `<Short>.paralog_cutoffs.tsv` (columns: `protein_ID, paralog_protein_ID, bitscore, evalue`). Proteins with no within-proteome paralog are omitted; `build_presence_matrix.py` falls back to `--default-evalue` for them.
   - `BUILD_PRESENCE_MATRIX` applies two paralog-aware filters:
     1. **Paralog-cutoff filter**: hit evalue must be < the query's paralog evalue (falls back to `--default-evalue` 1e-5 if no paralog detected).
     2. **Paralog-competition filter**: if the query's paralog hits the same target proteome with a better evalue, the hit is disqualified.
   - Produces `presence_matrix.tsv` (protein × proteome 0/1 matrix with `protein_id` and `source_proteome` columns) and `candidates.txt` (lines of `source_proteome::protein_id`).

3. **CLUSTER workflow** (`workflows/cluster.nf`):
   - `EXTRACT_CANDIDATES` reads `candidates.txt` and pulls matching sequences from all ingroup proteome FASTAs into `candidates.fa`.
   - `MMSEQS_CLUSTER` runs `mmseqs2 easy-cluster` (`--min-seq-id 0.3 -c 0.8 --cov-mode 0`) → `clusters_rep_seq.fasta`, `clusters_all_seqs.fasta`, `clusters_cluster.tsv`.
   - Emits: `candidates_fa`, `representatives`, `cluster_tsv`.

4. **VALIDATE workflow** (`workflows/validate.nf`):
   - `TBLASTN` searches all cluster representatives against each outgroup genome (one job per genome).
   - `SUMMARIZE_TBLASTN` merges per-genome TBLASTN TSVs, expands rep-level hits to cluster members, produces `tblastn_summary.tsv` (protein × outgroup-genome 0/1 matrix).

5. **ANNOTATE workflow** (`workflows/annotate.nf`):
   - Single `ANNOTATE_MATRIX` process runs (conditionally, based on whether paths are set):
     - `hmmsearch` vs Pfam-A → `candidates.pfam.tblout` (MPI-capable via `--hmm_mpi`).
     - `diamond blastp` vs SwissProt → `candidates.swissprot.tsv` (best hit, outfmt 6 with stitle).
   - `annotate_presence_matrix.py` merges hits into the matrix with annotation priority:
     1. Model organism gene names (via `--modelorgs_config` YAML).
     2. Pfam domain names (all unique domains per protein).
     3. SwissProt description (best hit, `sp|ACCN|ID` prefix stripped).
   - Adds columns: `gene_name`, `product_description`, `function_source`, `Best_Swissprot`, `Pfam_Names`.
   - Produces `presence_matrix.function.tsv`.

6. **SUMMARIZE workflow** (`workflows/summarize.nf`):
   - `MAKE_NOVELTIES` calls `make_novelties.py` to produce one `novelties.<SHORT>.tsv` per ingroup species.
   - A protein is a novelty candidate if:
     - `source_proteome == short` (originates from that species).
     - Present in ≥ `ingroup_min_frac` of all ingroup species.
     - Absent from every outgroup proteome.
     - No TBLASTN hit in any outgroup genome (unless `--skip_tblastn_filter` is passed — hits still reported in `tblastn_outgroup_hits` column).
   - Currently wired with `--skip_tblastn_filter` in `workflows/summarize.nf`.

7. **REPORT workflow** (`workflows/report.nf`):
   - `MAKE_REPORT` calls `make_report.py`, which merges the annotated matrix, the
     TBLASTN summary and the novelties tables into `report.html`.
   - `lib/report_data.py` builds the payload; `lib/report_template.py` holds the page.
   - The payload is embedded as JSON in a `<script type="application/json">` block, so
     the page opens from `file://` with no network access — reports get copied off the
     cluster. Any new asset must stay inlined.
   - Annotations come from the **matrix**; the novelties tables supply novelty status
     and protein sequences. When `--novelties` is omitted, novelty status is recomputed
     from the matrix (mirroring `make_novelties.py` with the TBLASTN filter skipped).
   - Sequences dominate file size, so `--sequences` (`params.report_sequences`) defaults
     to `novelties` — roughly 5 MB for a 21k-protein, 4.8k-novelty run.

8. **Final outputs** in `results/<project>/`: `presence_matrix.tsv`, `presence_matrix.function.tsv`, `candidates.fa`, `tblastn_summary.tsv`, `novelties.<SHORT>.tsv` for each ingroup species, and `report.html`.

## Key Parameters (`nextflow.config`)

| Parameter | Default | Notes |
|---|---|---|
| `--config` | (required) | Path to analysis CSV |
| `--data_dir` | (required) | Directory containing all FASTA files; also searched in `pep/`, `dna/`, `genome/`, `scaffolds/` subdirs |
| `--project` | CSV basename | Results subdirectory under `--outdir` |
| `--outdir` | `results` | Root output directory |
| `--run_tool` | `phmmer` | `phmmer`, `diamond`, or `blast` |
| `--evalue` | `1e-5` | Fallback e-value for proteins with no detectable paralog; also TBLASTN significance cutoff |
| `--parse_evalue` | `0.01` | Loose noise ceiling passed to `parse_hits.py`; final filtering uses paralog cutoffs |
| `--ingroup_min_frac` | `0.75` | Fraction of ingroup proteomes that must contain a hit |
| `--use_orthofinder` | `false` | Placeholder — OrthoFinder clustering not yet implemented |
| `--pfam_hmm` | `null` | Path to Pfam-A.hmm; skips Pfam annotation if unset |
| `--swissprot_dmnd` | `null` | Path to SwissProt `.dmnd`; skips SwissProt annotation if unset |
| `--modelorgs_config` | `null` | Absolute path to model organisms YAML; skips gene-name lookup if unset |
| `--report_sequences` | `novelties` | Which proteins carry a sequence in `report.html`: `novelties`, `all`, or `none`. Sequences dominate the file size |
| `--hmm_mpi` | `false` | Run hmmsearch with MPI (`mpirun -np <hmm_mpi_tasks> hmmsearch --mpi`) |
| `--hmm_mpi_tasks` | `null` | MPI task count; defaults to `max_cpus` when `hmm_mpi=true` |
| `--max_cpus` | `32` | Cluster-wide CPU cap |
| `--max_memory` | `64.GB` | Cluster-wide memory cap |
| `--max_time` | `24.h` | Cluster-wide time cap |

## Process Labels and Resources

| Label | Default CPUs | Default Memory | Use for |
|---|---|---|---|
| `high_cpu` | 8 (32 on SLURM) | 16 GB (32 GB on SLURM) | Diamond, Blast, mmseqs2, ANNOTATE_MATRIX |
| `med_cpu` | 4 | 8 GB | Self-searches |
| `hmmsearch` | 4 | 8 GB | hmmsearch (MPI-aware on SLURM) |
| `low_cpu` | 1 | 4 GB | Parse steps, matrix building, make_novelties |

All processes activate the pixi environment via `beforeScript` and symlink `db/` and `configs/` into the work directory.

## Running the Pipeline

### Minimal run (local executor)

```bash
nextflow run main.nf \
    --config configs/pezio4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool phmmer
```

### SLURM cluster

```bash
nextflow run main.nf \
    -profile slurm \
    --config configs/pezio4_asco.csv \
    --data_dir /path/to/fastas \
    --run_tool diamond \
    --ingroup_min_frac 0.8 \
    --pfam_hmm db/pfam/Pfam-A.hmm \
    --swissprot_dmnd db/uniprot/uniprot_sprot.fasta.dmnd \
    --modelorgs_config configs/modelorgs.yaml
```

### Built-in test profile

```bash
nextflow run main.nf -profile test
```

Sets `params.config`, `params.data_dir`, and `params.project` from `tests/data/test.csv` and uses the local executor.

### Resuming after failure

```bash
nextflow run main.nf -resume --config configs/... --data_dir ...
```

`-resume` reuses completed tasks from `work/`. The `search_cache/` directory is additionally protected by `storeDir` — those steps are never re-run even across separate invocations.

### Switching search tools

Pass `--run_tool phmmer|diamond|blast`. Results are keyed by tool name in the cache filenames (`*.phmmer.tblout.gz`, `*.diamond.tsv.gz`, `*.blast.tsv.gz`) so cached results from one tool are not reused when switching.

---

## Development Guide

### Where to add new code

| What you are adding | Where it goes |
|---|---|
| A new external tool | New file `modules/<toolname>.nf` |
| A new analysis stage | New file `workflows/<stage>.nf` |
| A Python helper called from a Nextflow process | `bin/<script>.py` (must be executable) |
| Shared Python logic | `lib/<module>.py` |
| Unit tests | `tests/test_<module>.py` |
| Test fixtures | `tests/data/` |
| A new analysis config | `configs/<name>.csv` |
| New pixi-managed dependencies | `pixi.toml` under `[workspace]` |

### Adding a module

1. Create `modules/<toolname>.nf`. One `process` block per file.
2. Name the process `SCREAMING_SNAKE_CASE` matching the file basename.
3. Add a `label` (`high_cpu`, `med_cpu`, `hmmsearch`, or `low_cpu`) so `nextflow.config` resource rules apply.
4. Use `storeDir` for expensive pairwise/self-search outputs; use `publishDir` only for final user-facing results.
5. Write a minimal test in `tests/` covering the Python helper if one exists.

### Adding a workflow

1. Create `workflows/<stage>.nf`. Import only the modules it directly calls.
2. Declare `take:`, `main:`, and `emit:` blocks — never omit `emit:`.
3. Import and call the new workflow from `main.nf`.

### Modifying `bin/` scripts

- Run `chmod +x bin/<script>.py` after creation.
- All arguments must be `argparse` flags. Positional arguments break when Nextflow reorders tokens.
- Import shared logic from `lib/` via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))`.
- The common parsed-hit TSV columns are: `query_id, target_id, evalue, bitscore, query_proteome, target_proteome`.
- `candidates.txt` format: `source_proteome_short::protein_id` (one per line).

### Modifying `lib/`

- `lib/` is a proper Python package (`lib/__init__.py` exists).
- Current modules:
  - `hits.py` — `PARSERS` dict (keys: `phmmer`, `diamond`, `blast`), `open_input()` (handles `.gz`), `Hit` namedtuple.
  - `config_parser.py` — `parse_config()` → `list[Sample]`; `short_to_group()`.
  - `fasta.py` — FASTA utilities.
  - `model_organisms.py` — `ModelOrgAnnotator`; three ID-transform strategies: `direct`, `strip`, `diamond_fasta`.
  - `report_data.py` — `build_payload()` → the report's embedded JSON; `ROW_FIELDS` names the row array layout.
  - `report_template.py` — `HTML_TEMPLATE`; substitutes `__PROJECT_TITLE__` and `/*__PAYLOAD__*/`.
- Keep `lib/` free of I/O side effects at import time.

---

## The interactive report (`report.html`)

Regenerate standalone from an existing results directory without re-running the pipeline:

```bash
bin/make_report.py \
    --matrix results/pezio4_asco/presence_matrix.function.tsv \
    --config configs/pezio4_asco.csv \
    --tblastn_summary results/pezio4_asco/tblastn_summary.tsv \
    --novelties results/pezio4_asco/novelties.*.tsv \
    --candidates_fa results/pezio4_asco/candidates.fa \
    --output results/pezio4_asco/report.html
```

### Constraints to preserve when editing the page

- **Self-contained.** No CDN scripts, external stylesheets, fonts, or fetches — the page
  must open from `file://` offline. Inline anything new.
- **Untrusted strings.** Protein IDs and product descriptions come from FASTA headers and
  SwissProt/Pfam text. Insert them with `textContent` (the `el()` helper), never
  `innerHTML`. `make_report.py` escapes `</` in the payload so an annotation cannot close
  the `<script>` block early — `tests/test_report_data.py` covers this.
- **Colour encodes evidence type, not group.** series-1 blue = presence from the protein
  search; series-2 green = TBLASTN genome hit. Ingroup vs outgroup is carried by column
  position and the band labels, so no hue does double duty. Both hues are validated for
  colour-blind separation and contrast in light and dark mode — re-validate before
  changing them.
- **Rows are virtualised** on a canvas (~20k rows). `drawGrid()` renders only the visible
  window; the scroll height comes from the spacer div.
- **The table tab is the accessibility twin** of the heatmap and must keep every value
  reachable without hovering.

### Testing template JS changes

`lib/report_template.py` has no test runner of its own — pytest only covers the Python
side (`tests/test_report_data.py` checks the embedded payload, not the JS behaviour). When
editing the JS inside `HTML_TEMPLATE`, verify it in two stages:

1. **Syntax check with `node --check`.** Extract the `<script>` body (the second one — the
   first is the `application/json` payload) and check it parses:

   ```bash
   python3 -c "
   from lib.report_template import HTML_TEMPLATE
   import re
   m = re.search(r'<script>\n(.*)</script>\n</body>', HTML_TEMPLATE, re.S)
   open('/tmp/report_js_check.js', 'w').write(m.group(1))
   "
   node --check /tmp/report_js_check.js
   ```

2. **Drive a real generated report with jsdom** to exercise DOM wiring end-to-end (filters,
   selects, detail-panel buttons, table rendering). This project has no JS dependencies
   checked in, so install jsdom into a scratch dir rather than the repo:

   ```bash
   npm install --prefix /tmp/novinv_check jsdom --no-audit --no-fund --silent
   ```

   Generate a small `report.html` from a hand-written fixture matrix/config (same shape as
   the fixtures in `tests/test_report_data.py`) via `bin/make_report.py`, then load it with
   `JSDOM({ runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true })`.
   jsdom does not implement `<canvas>` or `matchMedia`, and the report page calls both
   unconditionally on load — stub them in the `beforeParse(window)` hook *before* the
   page's inline script runs, or every event handler that calls `drawGrid()` throws:

   ```js
   beforeParse(window) {
     window.matchMedia = () => ({ matches: false, addEventListener() {}, removeEventListener() {} });
     const ctxStub = {
       setTransform() {}, clearRect() {}, fillRect() {}, strokeRect() {},
       measureText: (t) => ({ width: (t || '').length * 6 }),
       fillText() {}, save() {}, restore() {}, translate() {}, rotate() {},
       beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
     };
     window.HTMLCanvasElement.prototype.getContext = () => ctxStub;
   }
   ```

   From there, drive the page like a user: set `select.value` / `input.value`, dispatch a
   plain `new window.Event('change'|'input'|'click', { bubbles: true })` (jsdom doesn't need
   real `MouseEvent`/`InputEvent` for this page's listeners), await a short `setTimeout` for
   the debounced search box (140ms), and assert on `textContent` of `#count`, `#tbl-body`,
   and `#detail` rather than the canvas heatmap (canvas is stubbed to a no-op, so pixel
   output isn't observable — assert through the table tab instead).

   This caught a real behavioural requirement during the gene-family feature: with default
   filters, a family member that doesn't independently clear its own species' novelty
   threshold stays hidden until the detail panel's "Show family members" button is clicked,
   which must override (not stack with) the novelty-only filter. A payload/unit test
   wouldn't have exercised that filter-override interaction — only driving the actual DOM did.

   Clean up the scratch install and generated fixtures afterward (`rm -rf /tmp/novinv_check`);
   none of it belongs in the repo.

---

## Nextflow Linking Rules

### DSL2 import paths

```nextflow
// From main.nf — workflows are siblings of main.nf
include { SEARCH  } from './workflows/search'

// From a workflow — modules are one level up
include { PHMMER_SEARCH } from '../modules/phmmer'
```

Relative paths are resolved from the file containing the `include`. Never use absolute paths.

### Channel contracts

Every channel element is a tuple `[meta, ...files]`. The `meta` map always contains at least:

```groovy
[id: String, group: 'IN'|'OUT', species: String, strain: String, taxon: String]
```

When a process produces a pair result, the meta is `[id: "${meta_q.id}_vs_${meta_t.id}"]`. Downstream processes that only need the file should `.map { meta, file -> file }` before collecting.

### Emit and intake conventions

- Every workflow `emit:` block names its outputs. Callers reference them as `WORKFLOW_NAME.out.<name>`, never positionally.
- Processes use named outputs (`emit: tsv`, `emit: matrix`, etc.).
- When collecting across samples: `.map { meta, tsv -> tsv }.collect()`.

### storeDir vs publishDir

| Directive | Use for | Behaviour |
|---|---|---|
| `storeDir` | Pairwise + self-search results in `search_cache/` | Skips the process if the output file already exists, even across pipeline runs |
| `publishDir` | Final results in `results/<project>/` | Copies/links on every run; does not affect task execution |

Never apply both to the same process output.

### Parameter access in processes

Access global parameters via `params.*` in the `script:` block. Only pass per-sample-varying values through channel tuples.

---

## Nextflow Conventions

- **DSL2** everywhere (`nextflow.enable.dsl=2` at the top of every `.nf` file).
- Meta map: `[id: row.Short, group: row.GROUP, species: row.Species, strain: row.Strain, taxon: row.TaxonGroup]`.
- Cache filename convention: `${meta_query.id}_vs_${meta_target.id}.${params.run_tool}.<ext>.gz`.
- Self-search cache: `${meta.id}_vs_${meta.id}.${params.run_tool}.<ext>.gz`.

## Python Script Conventions

- All scripts in `bin/` must be executable (`chmod +x`) and have `#!/usr/bin/env python3`.
- Shared logic lives in `lib/`; import via `sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))`.
- Use `--input`/`--output` argparse flags; never hardcoded paths.

## Analysis Config CSV Format

```csv
GROUP,Species,Strain,Protein,DNA,Short,TaxonGroup
OUT,Neolecta irregularis,,Neolecta_irregularis.proteins.fa,Neolecta_irregularis.scaffolds.fa,Nirr,Taphrinomycotina
OUT,Schizosaccharomyces pombe,,Schizosaccharomyces_pombe.proteins.fa,Schizosaccharomyces_pombe.scaffolds.fa,Spom,Taphrinomycotina
OUT,Saccharomyces cerevisiae,S288C,Saccharomyces_cerevisiae.aa.fa,Saccharomyces_cerevisiae.dna.fa,Scer,Saccharomycotina
IN,Neurospora crassa,OR74A,Ncrassa.pep.fa,Ncrassa.dna.fa,Ncra,Pezizomycotina
IN,Aspergillus fumigatus,Af239,Afum.pep.fa,Afum.dna.fa,Afum,Pezizomycotina
IN,Pyronema omphalodes,CBS 144459,Pomp.pep.fa,Pomp.dna.fa,Pomp,Pezizomycotina
IN,Coccidioidies immitis,WA_211,Cocci_WA211.pep.fa,Cocci_WA211.dna.fa,Cimm,Pezizomycotina
```

- `Short` is a ≤8-char unique ID used in filenames and output tables throughout.
- `Strain` may be empty.
- The config CSV filename (without `.csv`) becomes the results output subdirectory name.
- `Protein` and `DNA` are basenames resolved relative to `--data_dir` (also checked under `pep/`, `dna/`, `genome/`, `scaffolds/` subdirs).

## Model Organisms YAML Format (`configs/modelorgs.yaml`)

Three `id_transform` strategies:

- `direct` — protein_id used as lookup key unchanged.
- `strip` — regex stripped from protein_id (e.g. `-T-p\d+$` → bare gene ID).
- `diamond_fasta` — protein_id → diamond hit → FungiDB protein_id → FASTA header gene field.

```yaml
model_organisms:
  - short: Ncra
    gene_names_csv: db/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv
    id_transform: diamond_fasta
    diamond_hits: db/modelorgs/Ncra_vs_FungiDB_Ncra.diamond.tsv
    protein_fasta: db/modelorgs/FungiDB-68_NcrassaOR74A_AnnotatedProteins.fasta
    fasta_gene_field: gene
  - short: Afum
    gene_names_csv: db/modelorgs/Aspergillus_fumigatus_gene_names_FungiDB.csv
    id_transform: strip
    strip_pattern: "-T-p\\d+$"
```

## Agent skills

### Issue tracker

Issues live in the repo's GitHub Issues via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Five canonical triage roles mapped to same-named GitHub labels. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context — one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
