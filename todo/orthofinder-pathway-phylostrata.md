# Phase 4 — OrthoFinder pathway + phylostrata

| Field | Value |
|-------|-------|
| **Date** | 2026-07-20 |
| **Author** | Jason Stajich |
| **Priority** | low |
| **Status** | open |
| **Category** | feature |
| **Related analyses** | novelty / loss / phylostrata |
| **Related data** | species tree (PHYling / nf_phyling) |

## Description

Add `--cluster_tool orthofinder`: a symmetric, tree-aware pathway that derives novelty,
loss, and core from one orthogroup matrix, and assigns each family a phylostratum
(gene-age node) using a species tree.

## Motivation

Cross-validation of the mmseqs-profile pathway with an independent clustering method, and
a principled phylogenetic framework: novelty = birth on the ingroup stem, loss = loss on
an ingroup branch (Dollo parsimony), contraction = length trend along ingroup branches.

## Proposed Approach

- `modules/orthofinder.nf` (one storeDir-cached run over all proteomes).
- `bin/orthofinder_to_matrix.py`: `Orthogroups.GeneCount.tsv` (prefer HOGs `N0.tsv`) →
  same `presence_matrix.tsv` + `candidates.txt` + cluster TSV, emitting **both** novelty
  and loss directions (collapses `LOSS_SEARCH` into a matrix re-filter).
- Species tree from OrthoFinder or PHYling; Dollo mapping → `phylostratum` column.

## Acceptance Criteria

- [ ] `--cluster_tool orthofinder` produces the standard pivot artifacts.
- [ ] Loss direction derived from the same run (no separate LOSS_SEARCH).
- [ ] Phylostratum column available to the reports.

## Notes

Prefer Hierarchical Orthogroups (`N0.tsv`) over flat `Orthogroups.tsv` to reduce
over-splitting. Depends on Phase 1's `--cluster_tool` branch and matrix contract.
