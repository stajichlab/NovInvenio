# Handoff: outgroup-side coverage floor + diamond sensitivity — design complete, ready to implement

**Date**: 2026-09-22
**Status**: Design and validation (Phase 0) complete. **Nothing implemented yet.** Written as a
handoff so implementation can continue in a fresh session, possibly with a different model.
**Read this whole document before writing any code.** It supersedes nothing in
`docs/superpowers/specs/2026-09-20-paralog-novelty-disqualifier-analysis.md` (a different,
related filter — see "How this relates to the paralog-rescue floor" below) but extends the same
line of work.

## One-paragraph summary

Issue #135 (closed) found diamond's default search mode misses real, divergent orthologs,
causing false novelty calls. The fix under design is **not** a single flag — it's an
**outgroup-side-only coverage floor** applied after a more sensitive search, because a more
sensitive search also risks false positives from domain cross-reactivity. An independent Opus
design review found the naive design (asymmetric search-mode by ingroup/outgroup) was wrong,
proposed a cheaper filter-side alternative instead, and specified a measurement gate before any
code. That gate (**Phase 0**) is now complete, with real numbers, on real data. The rule is
validated well enough to build, with one concrete correction already identified and confirmed
by the same data.

---

## 1. Background — why this exists

Novelty calling requires: a candidate protein present in `>= --ingroup-min-frac` (default 0.75)
of ingroup proteomes **and** absent from every outgroup proteome
(`bin/build_presence_matrix.py`, `--other-max-frac 0.0`). Presence in a given proteome is
decided by a diamond blastp hit passing an E-value cutoff.

**Issue #135** (closed 2026-09-21) found diamond's **default** search mode misses real,
divergent orthologs: on `NovInvenio_Investigations/results/pezizo_set1` (11 species), 841
candidate x outgroup-protein pairs had a TBLASTN genomic hit but **zero diamond hit** at
default sensitivity, even out to E=0.01. A targeted `diamond blastp --very-sensitive` re-search
recovered 711/841 (84.5%) as significant hits (median E=1.68e-20). `--ultra-sensitive` recovered
734/841 (87.3%), with the 23 additional recoveries showing no obvious contamination on first
pass (later shown to be too small a sample to conclude anything — see §3).

This led to **PR #149** (`--diamond_sensitivity`, merged, default `''` = unchanged behaviour)
and the recorded **DEFER** decision in `todo/diamond-very-sensitive-main-search.md`: don't adopt
a more sensitive default project-wide yet, because ~50% of candidate-count changes under
`--very-sensitive` are unvalidated (see that todo file's A/B benchmark).

**This document is the next step**: designing what would let that DEFER become an ADOPT —
specifically, a filter that lets a more sensitive search be used for outgroup-absence calling
without inheriting its false-positive risk.

## 2. How this relates to the paralog-rescue floor (issue #128/#136 — already shipped)

There are now **two different, unrelated-but-similarly-named "floor" mechanisms** in this
codebase. Do not conflate them:

| | paralog-rescue floor (SHIPPED, #136) | coverage floor (THIS DOCUMENT, not yet built) |
|---|---|---|
| Lives in | `bin/build_presence_matrix.py` filter 2 (paralog competition) | would live in the same file, as a **new**, later filter stage |
| Question it answers | "did the query's own in-genome paralog explain this hit away?" | "is this outgroup hit real orthology or a domain-driven artifact?" |
| Trigger | `--paralog-rescue-evalue` (default `1e-20`, ON) | not yet implemented — this doc proposes `qcov`, OFF by default |
| Motivating case | CorA/A7UWR3 vs its paralog Q7SEN2 | Ribosomal proteins / RAD52 in the BUSCO false-rejection set (§5) |

The paralog-rescue floor is already live and working. This document's coverage floor is
additive on top of it — see §6 "Build order" for exactly where it goes in the filter chain.

## 3. The design conversation and the Opus review (read the full review before building)

The original design conversation proposed three things:
1. A coverage-aware floor to reject domain-cross-reactivity hits.
2. Using `--ultra-sensitive`/`--very-sensitive` for outgroup-absence searching specifically.
3. **Asymmetric search** — less-sensitive search for ingroup-presence, more-sensitive for
   outgroup-absence — motivated by an asymmetric cost model: a false ingroup-presence hit is
   "cheap" (just promotes a candidate for review), a false outgroup-presence hit is "expensive"
   (silently, permanently excludes a real novelty).

An independent Opus-model review (full text preserved in this session's transcript; key points
below) found:

- **The cost-asymmetry framing was right but stated one-sidedly.** The ingroup side has TWO
  error directions, not one: a **missed** ingroup hit (default mode's own known failure) is
  JUST as silently fatal as a false outgroup-presence hit — it drops a candidate below
  `ingroup_min_frac` with nothing downstream to recover it. Two curated positives (`spa-9`,
  `ham-11`) are already known-excluded this way (`configs/controls/pezizo_set1.controls_verification.md`).
- **Therefore: do NOT build asymmetric search.** It would deliberately reinstate a known-defective
  sensitivity level on the ingroup side, using insensitivity as a proxy for "reject cross-reactive
  hits" — a bad proxy, since #135 already showed 84.5% of default-mode misses are genuine
  homology, not cross-reactivity.
- **A cheaper, better alternative exists: filter-side asymmetry.** `build_presence_matrix.py`
  already computes ingroup/outgroup columns from one dataframe in one pass
  (`ingroup_ids`/`outgroup_ids`/`other_ids` already in scope). Applying a coverage filter only to
  `target_proteome ∈ other_ids` rows is a ~3-line mask — versus splitting `modules/diamond.nf`'s
  `DIAMOND_SEARCH` into two invocations per query (it currently batches all targets, ingroup and
  outgroup, in ONE job specifically to keep job count down), which would roughly double job count
  and invalidate the entire `storeDir` cache. **Engineering cost ratio: roughly 1000:1 for the
  same effect.**
- **Do NOT mirror the coverage floor onto the ingroup side.** The error directions invert there:
  the floor firing wrongly on an ingroup hit would DROP a candidate below threshold (the
  expensive error), while failing to fire lets a cheap false positive through. Mirroring the
  rule converts the ingroup's cheap error into the expensive one.
  **Instead: report, don't filter** — a per-candidate sidecar count of ingroup cells carrying
  the coverage-floor signature, surfaced as a report column/sort key (same pattern as
  `--output-evalues`/`--output-targets`), so a reviewer can discount a weak 4/5-ingroup call
  without silently killing it. This also make the ~1400 unresolved "newly novel" candidates
  from the `--diamond_sensitivity` A/B benchmark (the thing blocking that DEFER decision)
  inspectable for the first time.
- **The two pieces of evidence for "ultra-sensitive adds no risk" don't hold up statistically.**
  43/711 and 43/734 share the same numerator — not two measurements. And "0 of 23 additional
  ultra-sensitive recoveries had the cross-reactivity signature" is unsurprising even under the
  null of no real difference (at a 6% base rate, P(0 of 23) ≈ 24%). **Do not carry forward the
  claim "ultra-sensitive is pure gain, no extra risk"** — it was never actually shown.
- **The coverage-floor rule itself (as originally proposed: `qcov<30% AND pident>=40%`) had
  real, unaddressed validation gaps**, most seriously: it was fitted on hits default mode
  *missed*, then proposed for use on ALL outgroup hits — a population-mismatch problem, since
  nobody had measured what fraction of *currently-counted* presence calls it would flip.
- **A curation-free validation set was available and unused**: universally single-copy
  `fungi_odb12` BUSCOs (known-true 1:1 orthologs) x outgroup proteomes. **This became Phase 0's
  gating measurement.**
- **Recommended build order** (adopted below): Phase 0 (measure, zero pipeline change) → Phase 1
  (outgroup-side floor, default off) → Phase 2 (ingroup-side reporting column, no filtering).
  **Do not build**: asymmetric search, an ingroup-side coverage *filter*, Smith-Waterman
  refinement (see §7 for why SW was also dropped).
- **One correction to the user's own cost-model argument, for the record**: the claim that a
  false ingroup-presence candidate gets "downstream scrutiny (TBLASTN validation, annotation,
  human curation)" before being reported isn't actually true in the code —
  `workflows/summarize.nf` wires `--skip_tblastn_filter` (TBLASTN is reporting-only, never
  gates), and `annotate_presence_matrix.py` only appends columns. The asymmetry argument still
  holds, but on "visible-and-reviewable vs silently-invisible" grounds, not on a review step
  that exists in the pipeline.

## 4. Phase 0 — measurement, now complete, with real numbers

All measurements below used real diamond/hmmscan/BUSCO runs against
`NovInvenio_Investigations/results/pezizo_set1` (the same dataset #135 used). **Raw data and
scripts are preserved at**
`NovInvenio_Investigations/results/pezizo_set1/coverage_floor_phase0/` (copied off node-local
`/scratch` before this session ended — that directory is now the durable source, not any
`/scratch` path referenced implicitly by earlier session transcript). BUSCO's own output
(`full_table.tsv` per species) is at
`NovInvenio_Investigations/results/pezizo_set1/busco_phase0/<species>.busco/run_fungi_odb12/`.

### 4.1 Shuffled-decoy null (addresses: does the E-value threshold itself manufacture spurious narrow/high-identity hits by chance?)

Composition-preserving shuffled decoys of the 334 candidate query proteins (position-randomized,
same amino acid composition) searched at `--very-sensitive` against a **full** outgroup
proteome (Ccin, 13,334 sequences — not a reduced target subset, which would understate diamond's
own database-size-scaled E-value).

**Result: 0 significant hits (E<1e-5) out of 30 raw hits total, out of ~4.46M possible
query x target comparisons.** Pure noise essentially never clears the significance bar even
once. The concern that the E-value cutoff manufactures a chance high-identity/low-coverage
population is not supported at this scale.

Files: `query_shuffled.fa`, `decoy_vs_full_Ccin.tsv`.

### 4.2 Domain-overlap mechanism check (does the "suspect" signature actually correlate with shared Pfam domains?)

Two versions run, because the first was uninformative:

- **Naive test** ("does the alignment span touch ANY Pfam domain on the query"):
  **97.7% vs 97.7%** for suspect (`qcov<30 & pident>=40`, n=43) vs a size-matched "broad"
  (`qcov>=50`, n=43) sample. Uninformative — nearly every protein has a domain somewhere, so
  "touches a domain" is true almost universally regardless of whether the hit is real orthology
  or cross-reactivity.
- **Sharper test** ("is the alignment CONFINED TO one domain's envelope, not just touching it"
  — computed as `min(overlap, domain_length) / alignment_span_length`, `>=0.7` = confined):
  **79.1% (34/43) suspect vs 48.8% (21/43) broad.** Real, meaningful discrimination — but not a
  clean separator. A substantial fraction of good "broad" hits (48.8%) are ALSO single-domain
  confined; real orthology in short/simple proteins often IS mediated by one conserved domain.

**Caveat**: Pfam database used was `/bigdata/stajichlab/shared/lib/funannotate2_db/Pfam-A.hmm`
(HMMER3/f, pressed 2025-09-08) — NOT `/bigdata/stajichlab/shared/db/Pfam-24/` (2009-era,
incompatible binary format with the installed HMMER 3.4; `hmmpress`-ing it fresh was not
attempted). If reproducing, use the funannotate2 copy or press a current Pfam-A release; do not
waste time on Pfam-24 as-is.

**A genuine bug was caught and fixed during this analysis**: the first pass indexed diamond's
wide-outfmt columns off-by-one (`qstart`/`qend` read from the `slen`/`qstart` columns instead),
producing biologically impossible spans (`qstart > qend`, which cannot happen in a
protein-protein alignment — no reverse strand). **If you rebuild any of this analysis, verify
column indices against a `head -1 | tr '\t' '\n' | cat -n` check before trusting coordinate-based
logic — this bug produced a plausible-looking but meaningless 14.0%-vs-11.6% result before being
caught.**

Files: `suspects.pfam.domtblout`, `vs_*.coords.tsv` (this run's correct-indexed coordinates),
`domain_overlap.py` (naive, buggy version kept for the record — do not reuse without fixing the
same indexing issue), `domain_overlap2.py` (correct, confinement-fraction version — reuse this
one).

### 4.3 Population-mismatch check (what fraction of ALL currently-counted outgroup presence hits would the rule flip — not just the #135-recovered subset)

Ran fresh: **all 3393 pezizo_set1 candidates** vs all 6 outgroup proteomes, at **default**
sensitivity (i.e. the hits that ALREADY count as presence in production today, before any
sensitivity change) — `E<1e-5`, wide outfmt.

**Result: 1231 significant hits total. 104 (8.4%) would be flagged by `qcov<30 & pident>=40`.**
This is a real, previously-unmeasured, non-trivial effect size — not a rare edge case.

qcov/pident distribution of this population: qcov p50=69.4 (p25=42.5, p75=86.6); pident p50=34.1
(p25=29.1, p75=39.8).

Files: `allcand_vs_*.tsv` (one per outgroup: Ccin, CneoH99, Mcir, Nirr, Scer, Spom).

### 4.4 Named diagnostic cases

Three genes flagged by the Opus review as adjudication cases, using
`configs/controls/pezizo_set1.controls_verification.md`'s own prior investigation:

- **`wsc`** (`sp|U9W802|WSC_NEUCR`): real broad hits, e.g. vs Ccin E=4.1e-48, qcov=59.9,
  pident=47.4. **Correctly NOT flagged** by the rule (qcov too high) — matches the doc's own
  conclusion that this is a real, broader WSC-domain family conservation, not lineage-restricted.
- **`so`/`soft`** (`tr|Q7SDK1|Q7SDK1_NEUCR`): real broad hits across 5 outgroups (qcov 36-44%,
  pident 25-30%). **Correctly NOT flagged** (pident always <40). Matches the doc's conclusion
  that `so/soft`'s clean absence-from-outgroup in the REAL (default-sensitivity) pipeline output
  was a sensitivity artifact, not genuine lineage restriction.
- **`spa-18`** (`tr|Q7SB60|Q7SB60_NEUCR`): its one disqualifying hit (vs Nirr,
  `tr|A0A1U7LH00|`, E=6.76e-08, **qcov=15.1%, pident=27.8%**) is narrow AND low-identity — the
  archetypal shape of a spurious short-motif match against an unrelated 1117aa protein (query is
  only 169aa aligned). **The originally-proposed rule MISSES this case entirely**, because
  `pident=27.8 < 40` fails the AND-gate. This is the concrete evidence that the rule's
  co-requirement of high identity is wrong — see §5's fix.

Files: `named_genes.fa`, `named_vs_*.tsv`.

### 4.5 BUSCO false-rejection rate — the actual gate (per Opus: "if this shows a non-trivial rate, the rule is dead")

Ran fresh BUSCO (v6.0.0, `fungi_odb12` lineage, protein mode) against pezizo_set1's own UniProt
proteomes for the 6 species actually used in #135's analysis (Ncra + Ccin, CneoH99, Mcir, Nirr,
Scer) — **not** the pre-existing `busco_pezizo5/` tables, which are keyed to a **different**
proteome (FungiDB IDs, not the UniProt IDs pezizo_set1 uses) and would have required redoing the
same fragile ID-crosswalk this session already did once for `configs/controls/pezizo5.controls.csv`
— reusing them risked contaminating the validation with mapping error.

**Cluster note for whoever runs this next**: the BUSCO jobs initially failed twice before
succeeding, for two DIFFERENT environment reasons, both worth knowing:
1. `/scratch/$USER/<jobid>` is **node-local**; a `sbatch` job can land on a different node than
   the one whose scratch path was hardcoded in the launch script. Always stage sbatch job I/O
   under `/bigdata`, never a scratch path from a different context (this is the same rule
   already in this repo's CLAUDE.md — it bit this analysis anyway).
2. Even after fixing (1), `module load busco/6.0.0` failed with `conda: error: unrecognized
   arguments: --stack` — a conda-state conflict from inherited environment variables. Fixed by
   `unset CONDA_SHLVL CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE _CE_CONDA _CE_M
   CONDA_PYTHON_EXE` and `module purge` before `module load busco/6.0.0`, and submitting with
   `--export=NONE` instead of `--export=ALL`.

**Universal single-copy BUSCO set for these 6 species: 483** (Complete + exactly 1 copy in
every species) — larger than the 315 found for the full 11-species pezizo5 set, as expected with
fewer species constraining the intersection.

Searched Ncra's copy of each of the 483 BUSCOs against each outgroup's full proteome DB at
`--very-sensitive`, wide outfmt. Of **2415 known-true-ortholog pairs** (483 x 5 outgroups):

- 2367 (98.0%) get a significant hit at all (2.0% miss even at `--very-sensitive` — a separate,
  smaller residual, not examined further here).
- **Of those 2367, the originally-proposed rule (`qcov<30 & pident>=40`) would wrongly reject
  41 (1.73%).**

**This is low, not zero — read as "proceed to build," not "rule is dead."** The 41 failures are
structurally explainable: several are ribosomal proteins and DNA-repair genes (`RM02`, `RAD52`)
— real, ancient, universally single-copy orthologs whose only strongly-alignable region is a
small, highly conserved core domain. That is a real, recognized failure mode of any
coverage-based rule, not evidence of a broken concept.

**Important caveat, not yet resolved**: BUSCO genes are ancient, highly conserved core genes —
systematically different from actual novelty candidates (by construction either genuinely
fast-evolving/lineage-specific, or erroneously flagged). **1.73% is a best-case lower bound**,
not necessarily the true false-rejection rate on the actual candidate population. Nobody has
measured the latter directly (doing so would need independently-verified true orthology among
actual candidates, which is circular without new evidence — e.g. a second, unrelated method of
confirming orthology).

Files: `universal_buscos.tsv` (the 483 x 6-species table), `busco_query.fa`,
`busco_vs_*.tsv` (one per outgroup). BUSCO's own raw output:
`NovInvenio_Investigations/results/pezizo_set1/busco_phase0/`.

## 5. The concrete fix already validated: drop the identity co-condition, use coverage alone

The `spa-18` case (§4.4) showed the AND-gate structure (`qcov<30 AND pident>=40`) misses narrow,
low-identity spurious hits. Re-scored ALL of §4.3 and §4.5's already-collected data (no new
search needed) against `qcov`-alone rules at three thresholds:

| rule | BUSCO false-rejection (§4.5, n=2367) | all-candidate flag rate (§4.3, n=1231) |
|---|---|---|
| `qcov<30 AND pident>=40` (original) | 1.73% | 8.4% |
| `qcov<30` alone | 3.34% (worse) | 17.0% |
| **`qcov<20` alone** | **1.18% (better)** | **11.0% (higher — catches more real signal)** |
| `qcov<15` alone | 0.38% (best false-rej.) | 7.7% (close to original) |

**`qcov<20` alone strictly dominates the original AND-gate rule**: lower false-rejection AND a
higher candidate-flag rate. It also catches `spa-18`'s case (qcov=15.1 < 20), which the original
rule missed entirely (confirmed directly: `qcov<20` alone evaluates `spa-18`'s hit as flagged;
`qcov<30 AND pident>=40` does not, since `pident=27.8 < 40` fails the AND-gate).

**Recommendation: drop `pident` as a co-condition. Use `qcov < 20%` alone as the rejection
rule.** `pident` may still be worth reporting alongside (it's informative for a human reviewing
a flagged case), but it should not gate the decision.

`qcov<15` is more conservative still (best false-rejection rate) at the cost of catching fewer
real spurious hits (7.7% vs 11.0% of the candidate population) — a reasonable alternative if
minimizing false rejection is prioritized over catching more noise. **This exact threshold
choice (15 vs 20) is the one open numeric decision left**, and it should probably be exposed as
a tunable parameter rather than hardcoded, given how directly the trade-off was just measured
and how easily it could be re-measured on another clade.

## 6. Build order (adapted from the Opus review, now with real thresholds)

### Phase 1 — outgroup-side coverage floor, default OFF

Location: `bin/build_presence_matrix.py`, as a new filter stage.

- New parameter: `--outgroup-coverage-floor-qcov` (suggest default `None`/`0` = off, matching
  the existing `--paralog-rescue-evalue` `0`-disables convention already in this file).
- Applied **only** to rows where `target_proteome ∈ other_ids` (the existing outgroup/other-group
  id set already computed in this script) — this is the filter-side-asymmetry mechanism from §3,
  NOT a search-mode change. No changes to `modules/diamond.nf` or the search step at all.
- Rule: reject (treat as absent) a hit where `qcov < <threshold>` (recommend defaulting the
  threshold to `20` when enabled, per §5 — but leave it a parameter, not hardcoded).
- **Requires #129's wide-outfmt fields to be present** (`qcov` column). Per `lib/hits.py`'s own
  documented policy (`Hit`'s docstring): **a run whose hit files lack geometry (narrow 4-column
  cache, or any phmmer `--tblout`-based hit) MUST fail loudly if this floor is requested, not
  silently no-op.** Decide and implement this explicitly — a silent no-op would be invisible and
  would corrupt cross-study comparisons (a study using phmmer would appear to have the floor
  "working" when it never fired).
- **Order it AFTER filter 2 (paralog competition), not before or merged into it.** Log the two
  mechanisms' rejections **separately** (e.g. distinct columns/counts) — otherwise a case like
  `hex-1` (already fully explained by filter 2) becomes unattributable between the two
  mechanisms, and the only clean cross-reactivity control this whole investigation has built up
  becomes useless for future validation.
- Needs the same TDD pattern already used for `--paralog-rescue-evalue`/`--paralog-rescue-delta`
  in this file (see `tests/test_build_presence_matrix.py` for the established fixture style) —
  hand-built hit tables with a known-real broad hit and a known-spurious narrow hit, verifying
  the floor rejects the narrow one and keeps the broad one, at both the default-off and an
  explicit-threshold-set state.

### Phase 2 — ingroup-side reporting column, NO filtering

Location: same file, or a small follow-on to it.

- Per-candidate count of ingroup presence cells whose supporting hit ALSO has `qcov < threshold`
  (same rule, same threshold parameter) — computed but **never used to disqualify or promote**
  a candidate.
- Surface as a new column in the matrix / a sidecar (same pattern as `--output-evalues`,
  `--output-targets`) → then into `novelties.html` as a report column and sort key
  (`lib/report_data.py`, `lib/report_template.py`).
- **Purpose**: makes the ~1400 unresolved "newly novel" candidates from the `--diamond_sensitivity`
  A/B benchmark (`todo/diamond-very-sensitive-main-search.md`) inspectable by a human reviewer
  for the first time — this is explicitly what that todo's DEFER decision is blocked on.

### Explicitly rejected — do not build

- **Asymmetric search** (different diamond sensitivity for ingroup vs outgroup targets). Wrong
  mechanism (§3), and ~1000x more expensive to implement than the filter-side equivalent for no
  additional benefit identified.
- **An ingroup-side coverage FILTER** (as opposed to the Phase 2 reporting column). Would invert
  the ingroup's cheap error into the expensive one (§3).
- **Smith-Waterman refinement for ambiguous cases.** See §7 — sharpens the wrong measurement.

## 7. Why Smith-Waterman was dropped

The original design included an optional SW refinement pass for pairs sitting near the
`qcov`/`pident` threshold boundary, scoped as cheap (compute estimate: negligible at the small
scale discussed — see the session transcript's compute-estimate exchange for the arithmetic, not
reproduced here since it's now moot).

**Dropped per the Opus review**, for reasons that survive independent of the compute cost:

1. SW gives an exact, non-heuristic alignment — but the deficit in this design is not measurement
   *precision* (diamond's heuristic qcov/pident are close enough to exact for this purpose), it's
   *accuracy of the qcov/pident-to-truth mapping itself* (whether narrow+low-coverage actually
   means "cross-reactivity" vs "real but short conserved domain" — exactly what §4.5's BUSCO
   false-rejection cases show it sometimes does NOT mean). Sharpening a noisy proxy doesn't fix
   the proxy.
2. "Ambiguous" (the band of cases that would get SW-refined) is itself an undefined second
   threshold needing its own validation — no one has measured how many cells would fall in it.
3. **A domain-aware tie-breaker that directly tests the actual mechanism already exists and is
   nearly free**: Pfam `hmmscan` both the query and target, and ask whether the aligned span
   sits inside a shared domain envelope — this is exactly what §4.2's `domain_overlap2.py`
   analysis already does, and `ANNOTATE_MATRIX` already runs `hmmscan` against Pfam-A for
   candidates in a real pipeline run. **If a tie-breaker is wanted, extend that, not SW.**

If SW is ever revisited, scope it as a standalone offline analysis script (like the ones in
`coverage_floor_phase0/`), not a pipeline stage.

## 8. Open questions for whoever picks this up

1. **Threshold**: `qcov<20` or `qcov<15`? §5's table gives the trade-off directly; this needs a
   human call (default-off either way, so it's low-stakes to start conservative at 15 and loosen
   later once Phase 2's ingroup reporting column gives more visibility into real-world impact).
2. **Does the false-rejection rate need re-measuring on a second clade** before this leaves
   "default off, opt-in" status? The BUSCO methodology here (fresh BUSCO run against the study's
   own proteomes, not a mismatched pre-existing table) is now a repeatable recipe — see §4.5's
   cluster-environment notes before rerunning it.
3. **The `--diamond_sensitivity` DEFER decision** (`todo/diamond-very-sensitive-main-search.md`)
   should be revisited once Phase 1 + Phase 2 both exist — the floor plus the reporting column
   together may be enough evidence to flip that to ADOPT. Not before.
4. **Issue/PR bookkeeping**: no GitHub issue exists yet for this Phase 1/Phase 2 work — file one
   before starting implementation, linking back to #135, #128/#136, #149, and this document.

## 9. Where everything lives

| what | where |
|---|---|
| This document | `docs/superpowers/specs/2026-09-22-coverage-floor-sensitivity-design-handoff.md` |
| Related, already-shipped paralog floor | `docs/superpowers/specs/2026-09-20-paralog-novelty-disqualifier-analysis.md`, PR #136 |
| The DEFER decision this work feeds into | `todo/diamond-very-sensitive-main-search.md` |
| `--diamond_sensitivity` flag (shipped, default off) | PR #149, `modules/diamond.nf`, `main.nf` |
| Wide-outfmt fields this depends on (shipped) | PR #151 (issue #129), `lib/hits.py` |
| Phase 0 raw data + scripts (persistent) | `NovInvenio_Investigations/results/pezizo_set1/coverage_floor_phase0/` |
| Fresh BUSCO runs used for §4.5 (persistent) | `NovInvenio_Investigations/results/pezizo_set1/busco_phase0/` |
| Named-gene ground truth used in §4.4 | `configs/controls/pezizo_set1.controls_verification.md` |
| Filter chain this extends | `bin/build_presence_matrix.py` |
