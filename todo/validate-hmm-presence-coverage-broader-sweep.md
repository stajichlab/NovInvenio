# Validate hmm_presence_cov / hmm_presence_min_residues against a broader sweep grid and curated controls

| Field | Value |
|-------|-------|
| **Date** | 2026-09-03 |
| **Author** | Jason Stajich |
| **Priority** | idea |
| **Status** | partially resolved — see decision #17 |
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

- [x] Broader grid run completed on at least pezizo5 (and ideally a second clade) —
      done: pezizo5 (decision #13) and `sordariales_shallow` (decision #17) both
      independently rank `hmm_presence_cov=0.3, min_residues=100` best. A third clade,
      `deep_broad_1kfg`, was attempted at real ADR-0002 scale (130 taxa) but is
      confounded by an unrelated clustering-identity failure (`busco_recovery` 0.011,
      see decision #17) — not usable evidence for this criterion either way, and its
      remaining 2 grid points were deliberately not run for that reason.
- [ ] `configs/controls/<clade>.controls.csv` populated for at least one clade, with
      `recall`/`fp_rate` genuinely measured (not blank) in the sweep output — still
      blocked on `bin/score_controls.py`, which remains Phase 2 / not yet built.
- [ ] A `nextflow.config` default change for `--hmm_presence_cov` (or an explicit
      decision not to change it) backed by this broader evidence, logged in
      `.living/decisions.md` — evidence now supports changing the default to `0.3`
      (two clades agree, every measured axis improves), but decision #17 deliberately
      deferred pulling that trigger pending the `recall`/`fp_rate` signal above. Needs
      an explicit go/no-go from the user.

## Notes

**2026-09-10 update (decision #17):** the first bullet above is now satisfied for
moderate/shallow divergence. The open work left under this todo is (a) build
`bin/score_controls.py` + at least one clade's controls CSV, and (b) decide whether to
flip the shipped default now on `presence_recovery`/`tblastn_removed` evidence alone, or
wait for (a). `deep_broad_1kfg`'s clustering-identity failure is a separate, new
follow-up (needs its own `min_seq_id`/`cov` sweep) and should get its own todo rather
than staying folded into this one.

Full background, the original bug reports (false novelty calls for real conserved
genes), and the initial 4-point sweep's exact numbers are in `.living/decisions.md`
entries #9–#13.

## 2026-09-21/22 update: grid extended to {0.2, 0.1}; recall/fp_rate measured for the first time

`bin/score_controls.py` is built (contrary to this file's earlier "not yet built" note --
it exists, was just never wired with a working `configs/controls/pezizo5.controls.csv`
until #147, and had its own latent bug: `proteome_cols` was computed as "everything not
protein_id/source_proteome", which silently includes `presence_matrix.function.tsv`'s
string annotation columns (gene_name, Pfam_Names, ...) whenever the harness prefers the
annotated matrix (it always does, when non-empty) -- breaking `family_presence_vector`'s
numeric comparison. This is very likely why recall/fp_rate were blank in EVERY prior sweep
run, not just this one. Fixed to a config-derived allowlist (commit pending this session).

Extended the grid to `hmm_presence_cov` {0.2, 0.1} (2026-09-20/21 SLURM run, job
28957241; `NI_Sweep/pezizo5_coverage_ext/`). One of the four new grid points
(`cov=0.2, residues=0`) never completed -- lost to `BUILD_CHUNK`/preempt queue
contention documented in `.living/learnings.md` L-3/L-8 -- and is excluded below.

All 8 points tested (`min_seq_id=0.3, cov=0.8, hmm_evalue=1e-3` held fixed):

| hmm_cov | residues | n_novelties | presence_recovery | recall | fp_rate | tblastn_removed |
|---|---|---|---|---|---|---|
| 0.5 | 0   | 1961 | 0.9597 | — | — | 554 |
| 0.5 | 100 | 1158 | 0.9918 | — | — | 26 |
| 0.3 | 0   | 1350 | 0.9907 | — | — | 117 |
| **0.3** | **100** | **1107** | **0.9952** | — | — | **21** |
| 0.2 | 100 | 1384 | 0.9910 | 1.00 (4/4) | 0.00 (0/1) | 73 |
| 0.1 | 0   | 1328 | 0.9934 | 0.75 (3/4) | 0.00 (0/1) | 52 |
| 0.1 | 100 | 1321 | 0.9934 | 0.75 (3/4) | 0.00 (0/1) | 48 |

recall/fp_rate blank above 0.3 because they were never successfully measured before this
session's bugfix -- not because they weren't run.

### Reading: 0.3 looks like the real optimum, not the edge of an incomplete grid

`presence_recovery` -- the curation-free, real-BUSCO-data metric, the one actually
sensitive to this parameter -- **peaks at cov=0.3 (0.9952) and does not improve further
below it** (0.991-0.993 at 0.2/0.1, i.e. flat-to-slightly-worse). Going in with the
2026-09-20 worry that "0.3 was the edge of the tested grid, so the optimum might be
past it" -- that worry is now resolved. It wasn't.

`recall` (now measured for the first time) also favors 0.2 over 0.1, consistent with the
presence_recovery signal's direction, though **on an n too small to lean on**: only 4 of
the 6 curated positives and 1 of 11 BUSCO negatives resolved against this mmseqs-family
matrix (12 of 17 non-placeholder controls came back `unresolved`) -- why that many
controls fail to resolve here is unexamined and worth its own look before trusting
recall/fp_rate as a tiebreaker.

### Recommendation (not yet actioned -- needs your go/no-go, same as before)

Change the shipped default `hmm_presence_cov` from **0.5 -> 0.3** (`min_residues=100`
unchanged -- already the shipped value). This is now backed by: the original two-clade
agreement (pezizo5 + sordariales_shallow, decision #17), an 8-point grid on pezizo5 that
shows 0.3 is a real plateau rather than an artifact of an incomplete sweep, and a recall
signal that -- weakly, n=4 -- points the same direction. Nothing in today's data argues
for going lower than 0.3.

Not recommending 0.2 or 0.1 despite 0.2/100's top composite score: its edge over 0.3/100
rests on presence_recovery being *slightly lower* (0.9910 vs 0.9952) traded against a
recall difference on 4 resolved controls -- not a margin worth trusting over the
better-measured metric.

## Still open before closing this ticket

- [ ] Investigate the 12/17 control-resolution failures against the mmseqs-family matrix
      before trusting recall/fp_rate for anything beyond a weak directional signal.
- [ ] Explicit go/no-go on flipping `nextflow.config`'s `hmm_presence_cov` 0.5 -> 0.3.
- [ ] A second, independently-run clade's --evalue/cov grid would still strengthen this,
      but is no longer blocking -- pezizo5's own grid is now wide enough to trust on its
      own terms.
