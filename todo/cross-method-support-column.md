# Phase 2 — Cross-method support column + concordance in reports

| Field | Value |
|-------|-------|
| **Date** | 2026-07-20 |
| **Author** | Jason Stajich |
| **Priority** | medium |
| **Status** | open |
| **Category** | feature |
| **Related analyses** | novelty / loss |
| **Related data** | — |

## Description

Phase 2 = validation + integration for the profile pathway: (a) the parameter sweep +
biological-control validation battery (ADR-0002 Q8), and (b) a cross-method `support`
column recording which pathway(s) called each candidate (`pairwise`, `mmseqs`, later
`orthofinder`), with concordance surfaced in the reports.

## Motivation

The sweep + controls turn the family-definition threshold from a guessed knob into an
evidence-backed default. Cross-method concordance is the cheapest convincing validation
once both pathways emit the same matrix: multi-method candidates are high-confidence;
disagreements flag threshold/boundary sensitivity. See `docs/adr/0002` "Resolved design".

## Proposed Approach

**Validation battery (Q8):**
- `bin/score_controls.py` reads `configs/controls/<clade>.controls.csv` (template already
  scaffolded), resolves each anchor (`protein_id` | `fasta` via hmmsearch | `busco`) to a
  family, compares its call to `expected_call` → recall (positives) + FP rate (negatives).
- BUSCO `fungi_odb12` single-copy recovery as the primary curation-free family-quality
  metric (recovered as single family, ~1/species); reuse PHYling/nf_phyling.
- Sweep grid `min-seq-id{0.2,0.3,0.4,0.5} × cov{0.5,0.7} × E{1e-3,1e-5}` scored on
  #families / #novelties / BUSCO recovery / neg-FP / pos-recall → pick the knee = default.
- Report TBLASTN-removal count as the annotation-artifact rate.

**Cross-method support:**
- Merge presence matrices from both pathways on `protein_id` / family.
- Add a `support` field to `ROW_FIELDS` in `lib/report_data.py`; render as filter + column.

## Acceptance Criteria

- [ ] Sweep produces a scored table; a defensible default is selected at the knee.
- [ ] `score_controls.py` reports recall/FP against the controls CSV.
- [ ] Each report row shows which method(s) support it; filter to concordant-only.

## Notes

Depends on Phase 1 (`profile-search-pathway.md`). Future GA/TC/NC per-family thresholds
(ADR-0002 Q5) validated via this same battery. Optionally split the sweep/validation into
its own todo if it grows.
