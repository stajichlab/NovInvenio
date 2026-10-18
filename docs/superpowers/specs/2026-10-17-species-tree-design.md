# Species tree in pairwise runs (main.nf) -- design

Status: implemented on branch `species-tree` (not yet merged). See "Changes during implementation".

## Intent
Each pairwise IN/OUT run shows a sequence-based tree of the species it compared,
on the run's `report.html` landing page. A reader sees how the IN, OUT and context
species relate before reading novelty/loss tables.

## Goals
- `SPECIES_TREE` process in main.nf: phyling align (markerset fungi_odb12) ->
  phyling filter -n 10 -> phyling tree -M ft (FastTree). One tree per run.
- Input: every samplesheet row with a protein file, all GROUPs (IN, OUT,
  NEAR_INGROUP, BROAD_OUTGROUP, DISCOVERY_*). Tip label = `Short`.
- Output: `species_tree.nwk`, `species_tree.markers.tsv` (markers kept, treeness/RCV,
  number of species per marker), `species_tree.svg`.
- `COLLATE_REPORTS` embeds the SVG inline in `report.html` (self-contained, no network).
- Runs in parallel with the search. No novelty/loss call reads the tree.

## Non-goals
- No change to pangenome.nf and its `--pangenome_species_tree` input.
- No change to novelty/loss calling.
- No claim that a 10-marker FastTree is a publication-quality species phylogeny.
  The report says "FastTree, 10 markers (top treeness/RCV), fungi_odb12".

## Anti-goals
- A run must not fail or stall because the tree step fails. On failure the report shows
  "tree not available: <reason>" and the run continues.
- A tree that looks authoritative but is mis-rooted. The report states how it was rooted.

## Rooting (my call; the OUT group may not be monophyletic)
1. Compute the unrooted tree.
2. Pick the edge whose split best separates OUT tips from all other tips
   (maximise OUT tips on one side minus non-OUT tips on that side; ties -> longest edge).
3. If that split is pure (all OUT, only OUT), report "rooted on OUT group (monophyletic)".
   Otherwise report "rooted on best OUT/non-OUT split: k of n OUT tips outside the
   OUT side; m non-OUT tips inside the OUT side", and colour the tips so the mismatch is visible.
4. If there are no OUT tips, midpoint root.
Taxonomy-based rooting is NOT in this version: the samplesheet has no Taxon_ID column.
It can be added later if the samplesheet gains one.

## Parameters (nextflow.config)
- `species_tree = true` (the step also needs `phyling_markerset_dir`; see below)
- `species_tree_markers = 10`
- `species_tree_markerset = 'fungi_odb12'`
- `phyling_markerset_dir` = directory holding the downloaded markerset (my call: runs
  do not download it themselves; compute nodes may lack network). Documented one-time
  `phyling download fungi_odb12`. The step runs only when this is set.
- `species_tree_container_version = '2.5.0'` (tag of the phyling image).

## Constraints
- Container: SPECIES_TREE runs in `ghcr.io/stajichlab/phyling:<species_tree_container_version>`; the main image is unchanged (see "Changes during implementation").
- phyling needs >= 4 species. With fewer, skip with a reason in the report.
- Size every job from measured time, not a guess; record runtime of the first real run.
- Pipeline repo is NovInvenio (code via PR). NII side: `sync_reports.sh` must keep the
  tree SVG in the committed `report.html` and put `species_tree.nwk` in `archive/`
  (class 2: does not grow with candidate count).

## Left to the builder
- SVG drawing method (Biopython/matplotlib is already used in NII's taxonomy_tree.py).
- Process labels/time limits.
- Test fixtures (tiny proteomes, a non-monophyletic OUT case).

## Changes during implementation
Controller rulings made after this draft:
- OUT = `GROUP == OUT` only. NEAR_INGROUP and other groups are not OUT.
- The markerset is passed to phyling as `<phyling_markerset_dir>/<markerset>/hmms`.
- SPECIES_TREE runs in `ghcr.io/stajichlab/phyling:2.5.0` instead of adding phyling to the
  main image. The main image, its version and `container_version` are unchanged.
- The cluster module `phyling/2.3.1` has a broken `filter` subcommand ("No module specified").
  phyling 2.5.0 (pip) was run on real data and passed align, filter and tree.
- A killed job (errorStrategy 'ignore') is covered by an `ifEmpty` stub status
  (`failed: species tree job did not complete`) so the report still renders.
- The step runs only when `--species_tree` is true and `--phyling_markerset_dir` is set.
  Otherwise no job starts and tracked stub statuses (`skipped: disabled`,
  `skipped: phyling_markerset_dir not set`) go to the report.
- The status file records `markerset:` and `n_non_out_inside`; the report shows both.
- Task 5 (NII side) needed no change.
- The container path was not run: no container runtime was available. Only phyling 2.5.0
  installed with pip was run on real data.
