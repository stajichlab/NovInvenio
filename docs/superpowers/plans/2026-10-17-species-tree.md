# Species Tree Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Each pairwise run (main.nf) builds a phyling/FastTree tree of all samplesheet species and shows it on `report.html`.

**Architecture:** A `SPECIES_TREE` process (phyling align -> filter -> tree) always emits its files, even on failure (a status file carries the reason). A pure-Python script roots the tree and draws an SVG. `COLLATE_REPORTS` passes the SVG and status to `make_index_report.py`, which embeds them. NII's `sync_reports.sh` keeps the SVG in the committed `report.html` and archives the Newick file.

**Tech Stack:** Nextflow DSL2, Python 3 (pytest), phyling 2.3.1, FastTree, Biopython or a hand-written Newick parser (builder's choice), pixi/Docker.

**Spec:** `docs/superpowers/specs/2026-10-17-species-tree-design.md`

## Global Constraints

- Markerset `fungi_odb12`; `filter -n 10`; `tree -M ft`. Param defaults: `species_tree = true`, `species_tree_markers = 10`, `species_tree_markerset = 'fungi_odb12'`.
- `phyling_markerset_dir` param: runs never download markers. One-time `phyling download fungi_odb12` is documented.
- All GROUPs with a protein file are tips. Tip label = samplesheet `Short`.
- Rooting: best OUT/non-OUT split edge (maximise OUT tips on the OUT side minus non-OUT tips on that side; ties -> longest edge). No OUT tips -> midpoint. The report states the rooting, and `k of n OUT tips outside the OUT side` when the split is impure.
- Failure never fails the run. The report shows `tree not available: <reason>`. Fewer than 4 species -> skipped with a reason.
- Report text names the method: `FastTree, 10 markers (top treeness/RCV), fungi_odb12`.
- No change to pangenome.nf or to novelty/loss calling.
- Report SVG is inline, no network.
- Work on a branch in `/bigdata/stajichlab/jstajich/projects/NovInvenio` (main is 40 commits behind origin; pull first). Code goes by PR.

## Review Focus

- OUT group not monophyletic: report shows the count outside; tips coloured by GROUP.
- No OUT rows in the samplesheet: midpoint root, no crash.
- Fewer than 4 proteomes: skipped, run completes, report says why.
- phyling finds few or no shared markers (distant OUT taxa): status carries the reason, run completes.
- `Short` names with characters that Newick or FastTree can mangle (`.`, `:`, `|`, spaces): labels round-trip to the report.

---

### Task 1: Rooting and SVG drawing (`bin/species_tree.py`)

**Files:**
- Create: `bin/species_tree.py`, `tests/test_species_tree.py`

**Interfaces:**
- Produces:
  - `root_on_outgroup(newick: str, groups: dict[str, str]) -> tuple[str, dict]`. Returns the rooted Newick and `{"method": "outgroup" | "outgroup_best_split" | "midpoint", "n_out": int, "n_out_outside": int}`. `groups` maps tip label -> GROUP.
  - `render_svg(newick: str, groups: dict[str, str]) -> str`. SVG text, tips coloured by GROUP, with a legend.
  - CLI: `species_tree.py --newick IN --config CSV --out_prefix species_tree`. Writes `species_tree.rooted.nwk`, `species_tree.svg` and `species_tree.rooting.json`. Reads GROUP/Short from the config using `lib/config_parser.py` GROUP aliases.

- [ ] **Step 1: Write failing tests** in `tests/test_species_tree.py`:
  - `test_monophyletic_out_roots_on_out_clade`: 6 tips, OUT = {o1, o2}; assert `method == "outgroup"` and `n_out_outside == 0`.
  - `test_non_monophyletic_out_reports_outside_count`: OUT = {o1, o2} split across the tree; assert `method == "outgroup_best_split"` and `n_out_outside >= 1`.
  - `test_no_out_tips_midpoint`: assert `method == "midpoint"`.
  - `test_labels_with_special_characters_roundtrip`: tips `A.fum_1`, `B|2`, `C 3` appear verbatim in the SVG text.
  - `test_render_svg_has_one_label_per_tip`.
- [ ] **Step 2:** Run `pytest tests/test_species_tree.py -v`. Expected: FAIL (module missing).
- [ ] **Step 3:** Implement the two functions and the CLI. Edge scoring as in Global Constraints. Make the file executable (`tests/test_bin_executable.py` checks this).
- [ ] **Step 4:** Run `pytest tests/test_species_tree.py tests/test_bin_executable.py -v`. Expected: PASS.
- [ ] **Step 5:** Commit `feat: species tree rooting and SVG`.

### Task 2: `SPECIES_TREE` process and params

**Files:**
- Create: `modules/species_tree.nf`
- Modify: `nextflow.config` (params), `main.nf` (call it, pass outputs to REPORT)
- Test: `tests/test_species_tree_wiring.py` (static checks, same style as `tests/test_alignment_metrics_wiring.py`)

**Interfaces:**
- Consumes: the `species_tree.py` CLI (Task 1); `samples_ch` in main.nf (`[meta, protein_fa, dna_fa]`).
- Produces: process `SPECIES_TREE(proteomes: path(list), config_csv: path)` emitting `svg` (`species_tree.svg`), `nwk` (`species_tree.rooted.nwk`), `markers` (`species_tree.markers.tsv`) and `status` (`species_tree.status.txt`; first line is `ok`, `skipped: <reason>` or `failed: <reason>`). All four files always exist (empty SVG/Newick on skip or failure). `publishDir "${params.outdir}/${Helpers.projectName(params)}"`.

- [ ] **Step 1:** Write a failing static test: `nextflow.config` defines the three params with the spec defaults; `modules/species_tree.nf` contains `phyling align`, `phyling filter`, `-n ${params.species_tree_markers}`, `phyling tree` and `-M ft`; `main.nf` includes `SPECIES_TREE`.
- [ ] **Step 2:** Run it. Expected: FAIL.
- [ ] **Step 3:** Implement the process. The script counts species first (fewer than 4 -> `skipped: fewer than 4 species`), checks `params.phyling_markerset_dir`, runs phyling with failures captured (never a non-zero exit), and writes `status`. When `params.species_tree` is false, still emit the files with `skipped: disabled`. Use an existing resource label.
- [ ] **Step 4:** Run the static test. Expected: PASS. Run the repo's existing config-parse smoke test (or `nextflow run main.nf --help`) to confirm the scripts parse.
- [ ] **Step 5:** Commit `feat: SPECIES_TREE process`.

### Task 3: Embed in `report.html`

**Files:**
- Modify: `bin/make_index_report.py` (new optional `--species_tree_svg`, `--species_tree_status`), `workflows/report.nf` (`REPORT` and `COLLATE_REPORTS` take the two new inputs), `main.nf` (pass them)
- Test: `tests/test_make_index_report.py` (extend the existing file if present, else create)

**Interfaces:**
- Consumes: Task 2 `svg` and `status` outputs.
- Produces: a `report.html` section `<section id="species-tree">`. When status is `ok`: the inline SVG, the method line and the rooting line. Otherwise: `tree not available: <reason>`. The section is omitted only when neither argument is given (old-run compatibility).

- [ ] **Step 1:** Write failing tests. Status `ok` -> output contains `id="species-tree"`, `<svg` and `FastTree, 10 markers`. Status `failed: no shared markers` -> output contains `tree not available: no shared markers` and no `<svg`. No arguments -> no `species-tree` id.
- [ ] **Step 2:** Run them. Expected: FAIL.
- [ ] **Step 3:** Implement in `make_index_report.py`. Wire the two path inputs through `REPORT` and `COLLATE_REPORTS`. The marker count in the method line comes from the markers/status file, not a hardcoded 10.
- [ ] **Step 4:** Run `pytest tests -k "index_report or species_tree" -v`. Expected: PASS.
- [ ] **Step 5:** Commit `feat: species tree on report landing page`.

### Task 4: Container, dependencies, docs

**Files:**
- Modify: `pixi.toml`, `Dockerfile`, `novinvenio.def`, `nextflow.config` (`container_version`), `README.md`, `CHANGES.md`

- [ ] **Step 1:** Add phyling 2.3.1 and FastTree to pixi, Dockerfile and def. Bump the version in all three and in `container_version`. Check the package names resolve (`pixi install --dry-run` or the repo's equivalent).
- [ ] **Step 2:** Document the params and the one-time `phyling download fungi_odb12` in README. Add a CHANGES entry.
- [ ] **Step 3:** Build or pull the image and run `phyling --version` and `FastTree 2>&1 | head -1` in it. Expected: both print a version. If the image cannot be built here, say so in the PR. Do not claim it works.
- [ ] **Step 4:** Commit `build: phyling and FastTree in image`.

### Task 5: NII side (`NovInvenio_Investigations`)

**Files:**
- Modify: `bin/sync_reports.sh`, `CLAUDE.md` (one line: the tree SVG in `report.html` is class 2; `species_tree.rooted.nwk` goes in `archive/`)

- [ ] **Step 1:** Read how `sync_reports.sh` copies `report.html` and `archive/`. Add copying of `species_tree.rooted.nwk` and `species_tree.markers.tsv` (gzipped) into `archive/`. Confirm `report.html` is not rewritten without the SVG.
- [ ] **Step 2:** Run the script against a completed pairwise study with the new files staged. Expected: `docs/.../report.html` still has `id="species-tree"`, and `archive/` holds both files.
- [ ] **Step 3:** Commit on a branch and open a PR (code and workflow changes go via PR).

### Task 6: Real run and timing

- [ ] **Step 1:** Download `fungi_odb12` once into the chosen `phyling_markerset_dir`.
- [ ] **Step 2:** Run a small pairwise study with the new image tag and `--species_tree true`. Record the real `SPECIES_TREE` wall time and CPU from the trace.
- [ ] **Step 3:** Open `report.html`. Check tip colours, the rooting note and the status line. If a study with a non-monophyletic OUT group exists, run it and check the note.
- [ ] **Step 4:** Put the measured times and shared-marker counts in the PR. Set time limits from these numbers, not a guess.

## Self-review notes

- Spec coverage: process (T2), outputs (T1, T2), report embed (T3), rooting (T1), failure and skip (T2, T3), params, container and markerset (T2, T4), NII side (T5).
- Not covered, by design: taxonomy-based rooting. The samplesheet can already carry `NCBI_TaxID` (read in main.nf's `uniprot_link_ch`), so a later version could use it.
