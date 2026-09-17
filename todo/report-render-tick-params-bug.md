# Fix `plot_classification_counts`'s `fig.tick_params` -> `ax.tick_params` bug

| Field | Value |
|-------|-------|
| **Date** | 2026-09-16 |
| **Author** | Jason Stajich |
| **Priority** | bug |
| **Status** | open |
| **Category** | pangenome report rendering |
| **Related analyses** | pangenome island/Pfam-enrichment report step |
| **Related data** | — |

## Description

`bin/pangenome_report_render.py`'s `plot_classification_counts()` calls
`fig.tick_params(axis="x", rotation=45)`, but `tick_params` is an `Axes`
method, not a `Figure` method — this raises `AttributeError` on any real
run where `classification_counts_dict` is non-empty (i.e. any run where
pair classification actually found physically-linked pairs).

## Motivation

Currently masked because every real smoke test run so far (both the prior
island/Pfam-enrichment plan's own smoke test and this report-enrichment
plan's) happened to produce zero pair-classification counts at the
10-strain test scale, so the buggy code path was never exercised in a
real execution. Will crash the report-render step on the first real study
where classification counts are non-empty (very likely — any pangenome
study with multiple strains and at least one significant island pair).

## Proposed Approach

One-line fix: `fig.tick_params(...)` -> `ax.tick_params(...)` in
`plot_classification_counts()`. Add a test with a non-empty
`classification_counts_dict` (the existing test suite apparently keeps
this fixture empty everywhere, which is exactly why this was never
caught).

## Acceptance Criteria

- [ ] Fix the `fig.tick_params` -> `ax.tick_params` typo
- [ ] Add a test exercising `plot_classification_counts` with a non-empty
      `classification_counts_dict` (would have caught this immediately)
- [ ] Confirm no other `plot_*` function in this file has the same
      `fig.` vs `ax.` mistake

## Notes

Found by a fix-round implementer while writing an unrelated test in the
2026-09-16 pangenome-report-enrichment plan (final-review fix round);
correctly left out of that plan's scope and flagged here instead, since it
predates that plan and isn't part of what it touched.
