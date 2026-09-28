# Per-group frequency bins and the `singleton` fix (#212, PR 1)

Status: **design agreed 2026-09-28** (grilling session; decisions posted to
[#212](https://github.com/stajichlab/NovInvenio/issues/212#issuecomment-5876428284)).
Nothing here is implemented. PR 2 (overlap figures: class x class heatmap,
UpSet, three-way bar) builds on this and gets its own spec.

## Problem

`bin/pangenome_frequency_bins.py` bins every family over the representative
ingroup strains, and `assign_bin()` returns `singleton` for `strain_count <= 1`.
So `singleton` also holds:

- **outgroup-only families** (ingroup count 0, present in the outgroup), and
- **families present only in non-representative ingroup strains** (dereplicated
  near-duplicates, mash < `pangenome_mash_threshold`).

`per_strain_summary()` then counts these for every matrix column. Measured on the
Coccidioides `rescue_freqpol_immitis_in_posadasii_out` run (pipeline `91e3157`):
29,502 `singleton` families, of which 7,808 are true ingroup singletons, 20,785
outgroup-only (319 carried by >= 90% of the 360 outgroup strains) and 909
non-representative-only. Per-strain `singleton` median: ingroup 37, outgroup 1,582.
Full numbers: #212.

## Goals

1. Every family has a correct class **within each group** (`GROUP` = ingroup,
   outgroup) of one run.
2. Per-genome counts, the per-genome figure and the outlier flags use the
   genome's **own group's** classes.
3. The report's composition section shows the two groups separately.
4. Consumers that read `bin` keep working without edits.

Non-goals: combining bins across runs (family IDs differ between runs: the two
reciprocal cocci runs share 29,755 of ~47.7k family IDs); binning by the
`Species` column; the overlap figures (PR 2).

## 1. Definitions

For a group `G` (ingroup or outgroup) of one run:

- `reps(G)`: strains with `GROUP == G` that are in the matrix and, when
  dereplication is on, have `is_representative == 1` in `strain_inventory.tsv`.
  Without dereplication, `reps(G)` is every `G` strain in the matrix.
- `nonreps(G)`: `G` strains in the matrix that are not in `reps(G)`.
- `count_G(f)`: number of `reps(G)` strains where `f` is present
  (`present` or `genome_only`, as `PresenceMatrix.is_present`).
- `freq_G(f) = count_G(f) / |reps(G)|`.

Class of family `f` in group `G`, in this order:

| Condition | Class |
|---|---|
| `count_G(f) >= 1` | `assign_bin(freq_G, count_G)`: core / soft_core / shell / cloud / singleton, cutoffs unchanged |
| `count_G(f) == 0` and `f` present in any `nonreps(G)` strain | `nonrep_only` |
| `count_G(f) == 0`, present in any strain of the other group | `outgroup_only` when `G` is the ingroup; `ingroup_only` when `G` is the outgroup |
| `count_G(f) == 0`, otherwise (only in strains of some third `GROUP` value) | `absent` |

`singleton` now means exactly one representative strain of the group. Strains
with any other `GROUP` value do not affect either group's classes. The rows are
checked top to bottom: a family in an ingroup non-representative **and** in the
outgroup is `nonrep_only` on the ingroup side, because it is in the ingroup.

## 2. `frequency_table.tsv` schema

Existing columns keep their names and meaning for the ingroup:

| Column | Meaning |
|---|---|
| `family` | unchanged |
| `frequency` | `freq_IN` (unchanged) |
| `strain_count` | `count_IN` (unchanged) |
| `bin` | ingroup class (section 1); values now include `nonrep_only`, `outgroup_only` |
| `frequency_out` | `freq_OUT`, 4 decimals; `-` when the outgroup is not binned |
| `strain_count_out` | `count_OUT`; `-` when not binned |
| `bin_out` | outgroup class (section 1); values include `nonrep_only`, `ingroup_only`; `-` when not binned |

The outgroup is **binned** when `|reps(OUT)| >= pangenome_outgroup_min_bin_strains`
(new param, default 3; the value is a judgement: the smallest size where core
(>= 95%) and cloud (< 15%) can differ). A run with no outgroup strains writes
the three `*_out` columns as `-`.

Readers must accept a table without the `*_out` columns (tables written before
this change) and treat it as "outgroup not binned".

### Consumers of `bin` (checked 2026-09-28)

| Consumer | Uses | Effect |
|---|---|---|
| `pangenome_cooccurrence.py:307` | `shell`, `cloud` | none |
| `pangenome_select_background_reps.py:27` | `shell`, `cloud` | none |
| `pangenome_build_islands.py:90` | `core`, `soft_core` | none |
| `pangenome_assembly_quality_qc.py:422` | `singleton` | set shrinks to true ingroup singletons: intended |
| `pangenome_neighborhood.py:93` (MODULE_NEIGHBORHOOD) | not `core`/`soft_core` = accessory | none (new labels are still accessory) |
| `pangenome_island_loci.py:114` -> `lib/island_locus.py:629` | tooltip text only | shows the new label (e.g. `outgroup_only`) |
| `pangenome_report_tables.py` / `pangenome_report_render.py` | all | changed in this PR |

Implementation must re-grep for any other reader of `bin` / `frequency_table`
before merge.

## 3. `per_strain_summary.tsv`

REPORT_TABLES gains two inputs: the samplesheet and `strain_inventory.tsv`
(stub when dereplication is off), plus `--ingroup_label` / `--outgroup_label`.

One row per matrix column. Columns:

`Short, group, is_representative, bins_from, n_families, core, soft_core, shell,
cloud, singleton, nonrep_only, outgroup_only, singleton_z, is_outlier`

- `group`: `IN`, `OUT` or the strain's other `GROUP` value.
- `is_representative`: `Y`/`N` (`Y` for all when dereplication is off).
- `bins_from`: which class column the counts use:
  - ingroup strains: `bin` (`in`);
  - outgroup strains: `bin_out` (`out`) when the outgroup is binned, else `bin` (`in`);
  - other-group strains: `bin` (`in`).
- Counts: families present in the strain, tallied by the class in `bins_from`.
  `outgroup_only` is non-zero only for strains counted with `bin` that carry
  ingroup-absent families (outgroup strains in the fallback, other-group
  strains). Families whose class is `absent` count only toward `n_families`. For a strain counted with `bin_out`, `ingroup_only` cannot occur
  (the strain is in the outgroup), so there is no `ingroup_only` column.
- `singleton_z` / `is_outlier`: modified z-score of `singleton`, with median and
  MAD computed **within the strain's group over its representative strains
  only**; every strain of the group is scored against them. `-` when the group
  has < 3 representatives or MAD == 0 (current rules).

## 4. Report changes

### Per-genome class figure (`per_genome_class_composition`)

- Segments: core, soft_core, shell, cloud, singleton, `nonrep_only`,
  `outgroup_only` (the last only appears in fallback / other-group bars).
- Each bar uses its row's `bins_from` classes, so the ingroup and outgroup
  blocks are two separate pangenomes. Block order and labels are unchanged
  (GROUP, then Species, then total).
- Non-representative genomes: `†` after the name when names are drawn
  (<= 150 genomes); a thin marker stripe at the left edge when they are not.
  The legend explains both.
- When the outgroup is not binned, the report states that outgroup bars use
  ingroup classes.

### Composition section

- Replace `core_shell_cloud_pie` with `group_composition`: one horizontal stacked
  bar per binned group (ingroup, then outgroup) on a shared family-count axis,
  segments core ... singleton.
- Per-group count lists under the figure, each class with count and percent of
  that group's families (families with a core ... singleton class in that
  group). `nonrep_only` and `outgroup_only` / `ingroup_only` are listed
  separately below each group, not in the bar or the percentages.
- A table without `*_out` columns gives one bar (ingroup).

### Unchanged

Accumulation curve, presence/absence matrix, islands, co-occurrence, enrichment
and their figures.

## 5. Pipeline wiring

- `nextflow.config`: `pangenome_outgroup_min_bin_strains = 3`.
- `FREQUENCY_BINS`: pass `--outgroup_label` and `--outgroup_min_bin_strains`.
- `REPORT_TABLES`: add `path(samplesheet)`, `path(strain_inventory)` inputs and
  the label args; `workflows/pangenome_profile.nf` passes the same
  `effective_samplesheet` / `strain_inventory` channels COOCCURRENCE already gets.
- Output docs that describe `frequency_table.tsv` or `per_strain_summary.tsv`
  are updated in the same PR.

## 6. Testing

TDD, unit level:

- class rules (section 1): each row of the table, including a family present
  only in an outgroup non-representative (`bin_out = nonrep_only`,
  `bin = outgroup_only`);
- outgroup below the minimum and no outgroup: `*_out` = `-`;
- dereplication off: `nonrep_only` never appears;
- `per_strain_summary`: `bins_from` per group, fallback, other-group strain;
  within-group z-score over representatives, scoring non-representatives;
- report render: two-bar composition, one-bar fallback on an old table, `†`
  marks, notes.

Real-data check: rebuild the frequency table offline for
`rescue_freqpol_immitis_in_posadasii_out` and confirm:

- `bin` counts: `singleton` = 7,808; `outgroup_only` + `nonrep_only` +
  `absent` = 21,694; `nonrep_only` >= 909 and `outgroup_only` <= 20,785 (#212
  counted families in both an ingroup non-representative and the outgroup as
  outgroup carriers; section 1 puts them in `nonrep_only`);
- `bin` for every family with `strain_count >= 1` is unchanged from the old table;
- co-occurrence / islands inputs (families with `bin` in shell, cloud, core,
  soft_core) are identical sets to the old table.

## 7. Rollout

1. This PR (fix). 2. PR 2 (overlap figures).
3. After both merge: offline update per run (frequency bins -> REPORT_TABLES ->
   assembly QC -> render), originals backed up, `run.json` records the
   report-step commit. Scope: the 3 published Coccidioides runs, then the 5 v2
   runs (Coccidioides + Afumigatus) when they finish. The known-issue notice on
   the live pages (NII `69a1c48`) is removed by that update.
