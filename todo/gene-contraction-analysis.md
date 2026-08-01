# Phase 3 — Gene-contraction analysis + contractions.html

| Field | Value |
|-------|-------|
| **Date** | 2026-07-20 |
| **Author** | Jason Stajich |
| **Priority** | medium |
| **Status** | open |
| **Category** | analysis |
| **Related analyses** | core / contraction |
| **Related data** | Pfam annotations, family MSAs |

## Description

Third analysis dimension: gene families present in **both** ingroup and outgroup where
ingroup members are systematically shorter (fewer AA) while maintaining domain
architecture — i.e., loss of inter-domain length, not domains.

## Motivation

Complements novelty (gene birth) and loss (gene death) with a "shrinking gene" signal.
The family MSAs from Phase 1 make it measurable at homologous positions.

## Proposed Approach

- Family set: matrix rows with high `ingroup_frac` AND high `outgroup_frac` (the "core").
- Stat: ingroup vs outgroup member length (median ratio, Mann–Whitney / Cliff's delta)
  over the family MSA (aligned-column occupancy, not raw length).
- "Domains maintained": require conserved Pfam architecture (reuse ANNOTATE hmmscan);
  rising domain-coverage fraction (domain AA / total AA) in the ingroup is the signal.
- New `build_contractions_payload()` + `contractions.html` following the
  `losses.html`/`core.html` single-table pattern (`lib/report_common.py`).

## Acceptance Criteria

- [ ] `contractions.<taxon>.tsv` + `contractions.html` produced from core families.
- [ ] Confounders controlled (see Notes) with documented null model.
- [ ] Report links out per the existing chain (Pfam/SwissProt/model-organism/NCBI).

## Notes

**Must control (robust-analysis convention):**
1. Genome-wide length shift — normalize against proteome-wide median length.
2. Phylogenetic non-independence — species-level means / comparative method, not pooled proteins.
3. Fragmented gene models — DNA-level length check; require consistency across many ingroup species.

Write a dedicated ADR for the contraction statistics + null model before building.
