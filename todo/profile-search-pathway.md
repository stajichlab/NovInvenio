# Phase 1 — PROFILE_SEARCH pathway (mmseqs → mafft → hmmbuild → hmmsearch)

| Field | Value |
|-------|-------|
| **Date** | 2026-07-20 |
| **Author** | Jason Stajich |
| **Priority** | high |
| **Status** | open |
| **Category** | feature |
| **Related analyses** | novelty / loss (family-based) |
| **Related data** | configs/pezizo5.csv (dev), configs/Chaetothyriales.csv (scale) |

## Description

Add a family-profile search pathway as an alternative producer of the presence-matrix
pivot, selected by `--cluster_tool mmseqs`. See `docs/adr/0002-family-profile-search-pathway.md`.

## Motivation

The pairwise pathway is `|IN| × (N−1)` jobs (~16,848 for Chaetothyriales) and does not
scale, especially for phmmer. Clustering the ingroup collapses the redundancy, reads
ingroup presence from membership, and probes the outgroup with sensitive family HMMs
(~128 hmmsearch jobs). Wires the currently-dormant `hmmbuild`/`hmmsearch` modules.

## Proposed Approach (settled by grilling — see ADR-0002 "Resolved design")

- Add **famsa** to `pixi.toml [workspace]`. New `modules/famsa.nf` (per-family MSA); wire
  `modules/hmmbuild.nf` + rework `modules/hmmsearch.nf` for the family-scan orientation
  (concatenated family-HMM db vs each proteome; target is a proteome, not SwissProt).
- New `workflows/profile_search.nf`:
  1. `mmseqs cluster` the ingroup (cascaded `-s 7`, `--min-seq-id 0.3 -c 0.8 --cov-mode 0
     --cluster-mode 0`); keep **≥2-member** families (min a param), drop singletons.
  2. famsa per family → `hmmbuild` → concatenated family HMM db.
  3. `hmmsearch --domtblout` the HMM db vs **all 209 proteomes** (both groups), fixed `--Z`.
  4. Presence per family per species = per-seq `E < 1e-3` AND profile-coverage ≥ 50%.
  - Emit `matrix`, `candidates`, `cluster_tsv` (same names/shape as `SEARCH`).
- `bin/` helper: cluster membership + hmmsearch presence → `presence_matrix.tsv` +
  `candidates.txt` (per-protein rows via family expansion).
- Branch in `main.nf` on `--cluster_tool` (default `pairwise` = unchanged); profile pathway
  passes the family `cluster_tsv` straight to `VALIDATE` (family-as-cluster; skip
  EXTRACT_CANDIDATES + MMSEQS_CLUSTER). TBLASTN **filters** here (not `--skip_tblastn_filter`).

## Acceptance Criteria

- [ ] `--cluster_tool mmseqs` produces `presence_matrix.tsv` identical in *shape* to the
      pairwise pathway on `pezizo5` (same columns, per-protein rows).
- [ ] Reports (`novelties.html`/`core.html`/`losses.html`) render unchanged from the new matrix.
- [ ] Scales to `Chaetothyriales.csv` within cluster resource limits.
- [ ] Default (`pairwise`) behaviour and outputs are byte-unchanged.
- [ ] Presence params (`--min-seq-id`, `--Z`, E-value, coverage, min family members) are
      exposed as sweepable `nextflow.config` params.

## Notes

Params above are the *shipped default*; the sweep + validation battery is Phase 2
(`cross-method-support-column.md`). Keep the "gene family" (whole-ingroup) vs "Cluster"
(candidate-level, ADR-0001) vocabulary distinct — needs a `/domain-modeling` glossary entry.
Fixed `--Z` value (summed target size vs constant 1e7) settled at implementation.
