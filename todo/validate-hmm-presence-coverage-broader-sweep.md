# Validate hmm_presence_cov / hmm_presence_min_residues against a broader sweep grid and curated controls

| Field | Value |
|-------|-------|
| **Date** | 2026-09-03 |
| **Author** | Jason Stajich |
| **Priority** | idea |
| **Status** | open |
| **Category** | validation |
| **Related analyses** | `.living/decisions.md` entries #12 and #13; `results/sweep_pezizo5_coverage/` |
| **Related data** | `configs/pezizo5.csv`; `busco_pezizo5/` (5 IN-group + 6 OUT-group BUSCO runs) |

## Description

A scoped 4-point sweep (`hmm_presence_cov` {0.5, 0.3} × `hmm_presence_min_residues`
{0, 100}, clustering params held at shipped defaults, `configs/pezizo5.csv`) found:

- `presence_recovery` (curation-free, from real outgroup BUSCO data —
  `bin/busco_presence_recovery.py`): 96.0% at the current shipped default
  (`cov=0.5, residues=0`) → 99.5% at the loosest tested setting (`cov=0.3,
  residues=100`).
- `tblastn_removed` contradiction rate (candidates with a genomic hit despite being
  called protein-absent): 28.2% → 1.9% over the same range.
- `busco_recovery` (clustering quality) was identical (0.753) across all 4 points, as
  expected — it's blind to these parameters.

This is real, encouraging evidence, but it is **only 4 grid points on one clade**, and
`recall`/`fp_rate` (from curated positive/negative gene controls) are still completely
unmeasured — no `configs/controls/pezizo5.controls.csv` exists yet, and it isn't obvious
what the right control genes are for this clade's ingroup/outgroup composition
(Pezizomycotina vs. Basidiomycota/Mucoromycotina/Taphrinomycotina/Saccharomycotina).
`nextflow.config`'s shipped `--hmm_presence_cov` (0.5) was deliberately NOT changed
pending this broader validation.

## Motivation

Before changing a project-wide shipped default that affects every `--cluster_tool
mmseqs` and `novelty_discovery` run, want more confidence that the recommendation
(`--hmm_presence_cov 0.3`) generalizes beyond 4 points and one clade, and a real
precision signal (`fp_rate`) to complement the false-absence-focused
`presence_recovery` metric — a setting that minimizes false absence but inflates false
presence would be a bad trade we currently can't see.

## Proposed Approach

- Broaden the grid, e.g. `hmm_presence_cov` {0.5, 0.4, 0.3, 0.2} ×
  `hmm_presence_min_residues` {0, 50, 100, 150, 200} — `bin/run_param_sweep.sh` already
  supports this via `HMM_COVS_LIST`/`MIN_RESIDUES_LIST` env vars, no code changes needed.
- Run on at least one additional clade (e.g. `sordario`, `Chaetothyriales`) to check the
  recommendation isn't specific to pezizo5's particular composition. BUSCO would need to
  be run for that clade's outgroup species too (fast — the 6 pezizo5 outgroup runs took
  ~1-4 min each on the `short` partition).
- Figure out the right curated controls: draft `configs/controls/pezizo5.controls.csv`
  (see `configs/controls/README.md` and `configs/controls/Chaetothyriales.controls.csv`
  as a template) — positive controls (known lineage-specific genes that SHOULD be
  flagged novel) and negative controls (known core/conserved genes that must NOT be).
  Open question, not yet resolved: what's an appropriate, defensible control set for
  this specific ingroup/outgroup composition? May be worth domain-expert input (known
  Pezizomycotina-specific gene families with independent literature support) rather than
  picking anchors ad hoc.
- Consider whether `bin/busco_presence_recovery.py`'s length-bucket breakdown (already
  computed, just not surfaced in the sweep metrics table) should be reported per-bucket
  in the sweep output, to check the improvement isn't concentrated in one length range.

## Acceptance Criteria

- [ ] Broader grid run completed on at least pezizo5 (and ideally a second clade)
- [ ] `configs/controls/<clade>.controls.csv` populated for at least one clade, with
      `recall`/`fp_rate` genuinely measured (not blank) in the sweep output
- [ ] A `nextflow.config` default change for `--hmm_presence_cov` (or an explicit
      decision not to change it) backed by this broader evidence, logged in
      `.living/decisions.md`

## Notes

Full background, the original bug reports (false novelty calls for real conserved
genes), and the initial 4-point sweep's exact numbers are in `.living/decisions.md`
entries #9–#13.
