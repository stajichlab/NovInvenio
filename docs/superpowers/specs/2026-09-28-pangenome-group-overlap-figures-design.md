# Ingroup-outgroup overlap figures (#212, PR 2)

Status: **draft 2026-09-28, for review.** Decisions from the #212 grilling session
(Q5, Q6, Q13; [#212 design comment](https://github.com/stajichlab/NovInvenio/issues/212#issuecomment-5876428284))
are fixed. Sections marked **Proposed** are this spec's own choices and are open
at review. Builds on PR 1 (`8c55750`,
`docs/superpowers/specs/2026-09-28-pangenome-group-bins-design.md`).

## Goal

Show how the ingroup's and the outgroup's pangenomes overlap, class by class,
within one run. Example: which *C. immitis* core families are absent from, rare
in, or core in *C. posadasii*.

## Input

`frequency_table.tsv` from PR 1: every family has an ingroup class (`bin`) and,
when the outgroup is binned, an outgroup class (`bin_out`).

For the overlap, each side's class is collapsed to 7 values:

| Side value | From `bin` (ingroup) | From `bin_out` (outgroup) |
|---|---|---|
| core, soft_core, shell, cloud, singleton | same label | same label |
| nonrep_only | `nonrep_only` | `nonrep_only` |
| absent | `outgroup_only`, `absent` | `ingroup_only`, `absent` |

Every family falls into exactly one (ingroup value, outgroup value) cell. The
cell (absent, absent) is always 0 for families in the matrix, but it is kept so
the grid is complete.

**Skip rule:** no overlap output (table, figures, report section) when the
table has no `bin_out` column or every `bin_out` is `-` (outgroup not binned,
PR 1 section 2). The report then says why, in one line.

## Real example (Coccidioides `rescue_freqpol_immitis_in_posadasii_out`, 2026-09-28)

Rows: *C. immitis* (IN, 133 representatives). Columns: *C. posadasii* (OUT, 278
representatives). 47,745 families; 48 non-empty cells.

| IN \ OUT | core | soft_core | shell | cloud | singleton | nonrep_only | absent |
|---|---|---|---|---|---|---|---|
| core | 5424 | 281 | 611 | 178 | 45 | 6 | 291 |
| soft_core | 207 | 39 | 174 | 46 | 9 | 1 | 133 |
| shell | 372 | 136 | 1052 | 787 | 204 | 33 | 1428 |
| cloud | 99 | 44 | 666 | 1866 | 523 | 155 | 3433 |
| singleton | 54 | 20 | 299 | 1068 | 647 | 327 | 5393 |
| nonrep_only | 2 | 2 | 25 | 264 | 137 | 59 | 909 |
| absent | 205 | 115 | 1794 | 7657 | 7519 | 3006 | 0 |

The implementation's real-data check must reproduce this table from that run's
`frequency_table.tsv`.

## 1. Table: `report_tables/group_class_overlap.tsv`

Written by REPORT_TABLES (`bin/pangenome_report_tables.py`), next to
`island_size_distribution.tsv`. Long format, all 49 cells, in the row/column
order above:

```
ingroup_class	outgroup_class	n_families
core	core	5424
...
```

When the skip rule applies, the file is written empty (0 bytes), the same
convention as the islands tables on runs without Pfam.

New output in `modules/pangenome/report.nf` REPORT_TABLES
(`path("group_class_overlap.tsv"), emit: group_class_overlap`) and a new
REPORT_RENDER input (`path(group_class_overlap)`), wired in
`workflows/pangenome_profile.nf`. `pangenome_report_render.py` gets an optional
`--group_class_overlap`; absent or empty skips the section.

## 2. Figures (`bin/pangenome_report_render.py`)

All three come from `group_class_overlap.tsv` and are written with
`_savefig_both` (PNG + PDF).

### 2a. Heatmap: `figures/group_class_overlap_heatmap`

- 7 x 7 grid, rows = ingroup value, columns = outgroup value, the order above.
- Each cell shows its count as text.
- Colour: log10(count + 1), sequential colormap; 0 cells blank (white).
- Axis titles name the group and, when every strain of a group has the same
  samplesheet `Species`, that species: e.g. "Coccidioides immitis (IN)". Otherwise
  "Ingroup (IN)" / "Outgroup (OUT)".

### 2b. UpSet plot: `figures/group_class_overlap_upset`

Drawn with matplotlib only (Q6; no `upsetplot` dependency).

- Sets: the 12 group-class sets (5 bands + `nonrep_only`, per group). The
  `absent` value is not a set.
- One intersection per non-empty cell except (absent, absent). A cell (a, b)
  with both non-absent is the intersection of set `IN:a` and set `OUT:b`; a cell
  with one side absent is that single set alone ("only in one group").
- Layout: bar chart of intersection sizes on top, dot matrix below (12 rows, a
  filled dot per member set, a line joining the two dots of a pair), set-size
  bars on the left.
- **Proposed:** order intersections by size, largest first, and show all of
  them, with no cap. At most 48 exist with this collapse (6 x 6 pairs + 6
  ingroup-only + 6 outgroup-only); the example has all 48.

### 2c. Three-way bar: `figures/group_class_overlap_shared`

Q5 decided: shared / ingroup-only / outgroup-only, each split by class.

- **Proposed:** three horizontal bars on one family-count axis:
  - *shared* (both sides non-absent), segments by **ingroup** class;
  - *ingroup only* (outgroup absent), segments by ingroup class;
  - *outgroup only* (ingroup absent), segments by outgroup class.
  The shared bar uses the ingroup's classes only; the heatmap carries the full
  pairing. Colours are the per-genome figure's (`BAND_COLORS`).
- In the example: shared = 15,862 families, ingroup only = 11,587, outgroup
  only = 20,296 (sum 47,745) (outgroup-only includes 3,006 `nonrep_only` on the outgroup side).

## 3. Report section

New section "Ingroup vs outgroup content", after "Per-strain summary":

- one sentence per figure saying what it shows;
- the three figures in the order heatmap, three-way bar, UpSet;
- a line pointing to `report_tables/group_class_overlap.tsv`;
- when skipped, one line with the reason (no outgroup columns / outgroup not
  binned).

## 4. Testing

TDD, unit level:

- collapse rules (table above), including `nonrep_only` on both sides and
  `outgroup_only` / `ingroup_only` -> absent;
- 49-row table in order, counts sum to the family count;
- skip rule: no `bin_out` column; all `-`; empty file written;
- intersections for the UpSet: pairs vs singles, (absent, absent) excluded;
- three-way totals;
- render: section present with a table, skipped with the reason without;
  each figure file written.

Real-data check: reproduce the example table above; look at all three PNGs.

## 5. Rollout

Merge, then the same offline update as PR 1 for the 3 published Coccidioides
runs (REPORT_TABLES + REPORT_RENDER only; frequency bins are unchanged), and
the 5 v2 runs when they finish.
