# Evaluate --very-sensitive for the main pairwise DIAMOND_SEARCH

| Field | Value |
|-------|-------|
| **Date** | 2026-09-10 |
| **Author** | Jason Stajich |
| **Priority** | medium |
| **Status** | open |
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
