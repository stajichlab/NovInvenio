## Priority bumped medium -> high (2026-09-20)

Issue #135's TBLASTN false-positive investigation found this is not a hypothetical
concern: 841/1479 (56.9%) of candidates with only genome-level (TBLASTN) evidence
get literally ZERO diamond hits at default sensitivity, even to E=0.01, at loci
where TBLASTN finds a significant hit and the protein exists in the outgroup
FASTA. Full decomposition: https://github.com/stajichlab/NovInvenio/issues/135#issuecomment-5753915888

---

# Evaluate --very-sensitive for the main pairwise DIAMOND_SEARCH

| Field | Value |
|-------|-------|
| **Date** | 2026-09-10 |
| **Author** | Jason Stajich |
| **Priority** | medium |
| **Status** | complete |
| **Category** | validation |
| **Related analyses** | NovInvenio_Investigations/studies/fungi/pezizo_set1 (HEX1_NEUCR novelty-candidate investigation) |
| **Related data** | results/pezizo_set1/search_cache/Ncra_vs_Ncra.diamond.tsv.gz, ad hoc Ncra-vs-Mcir benchmark (not archived) |

## Description

`DIAMOND_SELF` (`modules/self_search.nf`) was switched from diamond's default
sensitivity mode to `--very-sensitive` (2026-09-10) because default mode failed
to detect the true within-genome HEX-1/eIF-5A paralog pair in *N. crassa* at
all (zero hits even at `-e 100`), which broke the paralog-competition filter
and caused a spurious HEX-1 "present in Mcir" call via a distant cross-hit to
Mcir's eIF-5A — see `.living/decisions.md` entry for that fix. The main
pairwise search (`DIAMOND_SEARCH`, `modules/diamond.nf`) still runs default
mode. Open question: should it also move to `--very-sensitive`?

## Motivation

A real-data benchmark (Ncra vs Mcir proteomes, `e<=0.01`, same as
`params.parse_evalue`) showed:

| | default | `--very-sensitive` |
|---|---|---|
| wall time | 2.4s | 9.1s (~3.9x) |
| Ncra queries with >=1 Mcir hit | 4,078 | 5,379 |

1,302 additional Ncra proteins (~13% of the proteome) get a Mcir hit under
`--very-sensitive` that default mode misses entirely, plus 362 more queries
keep a hit but change which Mcir protein is the best match. This would likely
reduce novelty/loss candidate counts noticeably project-wide (better
outgroup-presence detection = fewer false "absent" calls) but at ~4x the
dominant cost of every study (the pairwise search is O(species^2) and is
already the bulk of runtime for studies with 50-100+ species). It would also
invalidate current results for all six existing studies, same as the
`--paralog-competition-scope`/filter fixes logged in `.living/decisions.md`.

## Proposed Approach

- Benchmark on a full study run (not just one species pair) to get a real
  wall-clock delta at project scale, not just one 9,759x12,143 pair.
- Quantify how novelty/loss candidate counts shift for at least one existing
  study (e.g. re-run `pezizo_set1` with `--very-sensitive` on `DIAMOND_SEARCH`
  and diff `novelties.*.tsv`/`loss_candidates.txt` against the current run).
- Check whether newly-surfaced cross-hits reintroduce paralog-cross-reactivity
  false-absence cases (the same failure class just fixed for self-search) --
  i.e. whether the paralog-competition filter (already running at
  `--very-sensitive` via the self-search fix) is sufficient to keep those in
  check, or whether more spurious "presence" calls slip through.
- Decide per-study or project-wide default; if adopted, needs a
  `.living/decisions.md` entry and a full re-run pass across all studies (this
  repo's own project + every `NovInvenio_Investigations` study), same process
  as prior presence-matrix-affecting fixes.

## Acceptance Criteria

- [ ] Full-study wall-clock benchmark recorded (not just a single species pair)
- [ ] Candidate-count diff for at least one study, characterized qualitatively
      (spot-check a sample of newly-flipped-to-present genes for plausibility)
- [ ] Explicit decision recorded (adopt / defer / reject) with rationale

## Notes (2026-09-10 update)

Concrete real-data example surfaced while curating `configs/controls/pezizo_set1.controls.csv`:
Neurospora's `so`/`soft` gene (NCU02794, hyphal-fusion regulator, literature claims
absent from Basidiomycota/yeast) is called cleanly absent from all 6 `pezizo_set1`
outgroup proteomes by the real (default-sensitivity) `DIAMOND_SEARCH` pathway -- but a
standalone `diamond blastp --very-sensitive` check of the same query against the same
outgroup `.dmnd` databases finds strong hits (e-41 to e-58) in CneoH99, Ccin, Mcir, and
Scer. Default mode misses all four entirely. This is very likely a false "novel" call
caused by main-search insensitivity, not genuine lineage-restriction -- the same failure
class as the `DIAMOND_SELF`/hex-1 fix, but on the presence-calling side rather than the
paralog-competition side. Strengthens the case for benchmarking `--very-sensitive` on
`DIAMOND_SEARCH` project-wide (see `configs/controls/pezizo_set1.controls_verification.md`
for the full comparison).

## Notes

Diamond sensitivity modes tested on the real HEX1 case (Ncra self-search):
`--sensitive`/`--more-sensitive` still missed the HEX-1/eIF-5A pair entirely;
`--very-sensitive` and `--ultra-sensitive` both found it (E=2.8e-06,
`--ultra-sensitive` gave the same result as `--very-sensitive` in that test,
so `--very-sensitive` was chosen as the less extreme option for the
self-search fix).

## Benchmark results (2026-09-20) — A/B, pezizo_set1, 11 species

Two fresh runs from the SAME code (branch `diamond-sensitivity-param`, PR #149); the only
difference is `--diamond_sensitivity very-sensitive`. Compared against each other, NOT against the
old `results/pezizo_set1`, which predates #136/#139 and would confound two changes. Launcher,
comparison script and full report: `NI_Sweep/diamond_sensitivity/{run_benchmark.sh,compare.py,
decompose.py,comparison_report.md}`. Both runs: 211 tasks, 0 failed. The flag was confirmed present
in the executed `DIAMOND_SEARCH` commands (default 0 of 5 tasks, very-sensitive 5 of 5).

### 1. Wall-clock — cheap in absolute terms, ratio NOT reliably measured

| process | sum default (s) | sum very-sens (s) | ratio |
|---|---|---|---|
| DIAMOND_SEARCH (n=11 each) | 294 | 481 | 1.64x |
| DIAMOND_SELF (identical work both arms) | 48 | 33 | 0.69x |
| DIAMOND_MAKEDB (identical work both arms) | 34 | 18 | 0.54x |

The two control rows do identical work in both arms yet differ by 30-46%, so run-to-run/node noise
at this scale is about that size and the 1.64x is barely outside it. Whole-study summed task time
0.23 h vs 0.30 h. **n=11 species cannot speak to the ~3.9x per-pair estimate above or to O(species^2)
scaling at 50-100 species** -- that acceptance criterion is only partly met.

### 2. Candidate diff — the net change hides ~50% gross turnover

| | default | very-sensitive | shared | only-default | only-very |
|---|---|---|---|---|---|
| novelty | 3515 | 3193 (-9.2%) | 1755 | 1760 | 1438 |
| loss | 101 | 141 (+39.6%) | 28 | 73 | 113 |

Two separate mechanisms, each ~98-99% clean:
- **1747 of 1760 dropped** because a real OUTGROUP hit appeared (outgroup 0 -> 1..4), ingroup
  unchanged. Spot-check: MED17 (Mediator; Scer 3e-21, Spom 2e-15), SLX4, CTRA2 -- recognizably
  conserved genes that default mode called "Pezizomycotina-specific". Strongly supports this half.
- **1412 of 1438 newly novel** because INGROUP presence rose (3->4, 2->4, 0->4 ...) with outgroup
  still 0. 290 of them had NO matrix row at all in default mode (no cross-species hit anywhere).
  **UNRESOLVED**: real lineage-specific genes default mode could not see in divergent ingroup
  members, or promiscuous family cross-hits inflating ingroup presence past the 75% threshold?
  The controls cannot say (all 5 positives were already candidates). Distinguishing them needs
  alignment coverage -- i.e. issue #129.

### 3. Controls — no paralog cross-reactivity regression

HEX1, LAH, ADA1, HAM8, SPA1 remain candidates under very-sensitive. CorA/A7UWR3 (hard negative) is
excluded in BOTH arms, so #136's rescue floor is confirmed end-to-end on real data (it IS a
candidate in the pre-#136 `results/pezizo_set1`).

### Caveats
- One study, 11 species; a single run per arm.
- Unexplained: a fresh DEFAULT run differs from the old baseline by 213 only-old (exactly the
  number the 1e-20 floor removes -- good) but ALSO 335 only-new (3515 vs 3393 - 213 + 335). Cause
  not investigated; it does not affect the default-vs-very comparison but should be understood.
- "Ingroup presence rose" is not evidence the newly-novel genes are real.

### Recommendation (decision: DEFER, see checklist below)
**Defer adopting a project-wide default; do not reject.** The outgroup-side gain is large, well
supported, and cheap; the ingroup-side gain is unvalidated. One option worth testing: an
**asymmetric** rule -- very-sensitive search for the "absent from outgroups" call (where a false
absence is the failure), default-mode presence for the "present in >=75% of ingroup" call (where
cross-family false presence is the failure). On these numbers that would drop the 1747 false
novelties while adding none of the unvalidated 1438. Untested. Sequence after #129 so the wider
`outfmt` and this cache invalidation share one re-run pass.

- [x] Full-study wall-clock benchmark recorded (ratio not reliably measured; see caveats)
- [x] Candidate-count diff characterized, spot-checked in both directions
- [x] Explicit decision recorded: **DEFER** (2026-09-20, Jason Stajich) -- do not adopt a project-wide default yet; revisit after #129 (alignment coverage) lets the 1438 newly-novel candidates be validated. Not rejected.


## Resolution (2026-09-23)

Adopted as the default (issue #171). A full-study benchmark covered 3 clades × default / sensitive / very-sensitive, plus a qcov-15 floor and a simulated stricter other-group E-value (NII `notes/diamond-sensitivity/README.md`):
- Share of candidates with an other-group TBLASTN hit: 29-69% (fast mode) → 9-33% (very-sensitive).
- Controls unchanged; search cost ≤1.7×.
- The qcov-15 floor adds little. A stricter other-group E-value is worse at every threshold. Keep `--evalue 1e-5`.
- Existing studies change on rerun.
