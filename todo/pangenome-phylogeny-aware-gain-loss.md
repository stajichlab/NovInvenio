# Gain/loss polarization should use a real strain phylogeny, not outgroup presence/absence counts alone

- **Priority**: high
- **Status**: open
- **Category**: methodology / pangenome co-occurrence
- **Date**: 2026-09-19
- **Author**: Jason Stajich
- **See**: `bin/pangenome_cooccurrence.py::polarize_direction()`, `bin/pangenome_assign_clades.py`, coccidioides_pangenome `genus_vs_ureesii` run

## What

Current gain/loss direction (`direction_a` in `cooccurring_pairs.tsv`/
`pair_classification.tsv`) is called by `polarize_direction()` purely from
raw presence/absence counts in the outgroup: present in every outgroup
strain → "loss" if absent in an ingroup strain; absent from every outgroup
strain → "gain" if present in an ingroup strain. This has no phylogenetic
context at all — it can't distinguish a single gain/loss event from
multiple independent ones, can't place the event on a specific branch, and
treats every ingroup strain as equally informative regardless of how
closely related it is to any other.

`clade_composition`/`bin/pangenome_assign_clades.py` is a mash-distance +
PCoA clustering scheme (confirmed: no real bifurcating tree is built
anywhere in this pipeline) — a reasonable cheap proxy for stratifying the
permutation null test, but not a substitute for a real phylogeny when the
actual question is "on which branch did this gene family appear or
disappear."

## Why

Surfaced during the `genus_vs_ureesii` real run (2026-09-18/19): the
99.3%-gain / 0.7%-loss split, and the Leiden trans-module/island findings
built on top of it (NACHT/HET, PKS islands), are real presence/absence
signal but their gain/loss *direction* and event *count* are only as good
as a single-outgroup-genome comparison. A real phylogeny would let these
claims be made properly: map presence/absence onto tree branches via
ancestral state reconstruction (Dollo parsimony is the standard
assumption for gene families -- genes are lost more easily than gained by
convergent means -- or full ML/Bayesian ancestral state reconstruction),
count actual independent gain/loss events per family instead of a single
binary call, and resolve the multi-outgroup case's "ambiguous" pairs
precisely against the tree instead of a flat presence-count rule.

## How to apply

1. Build a real strain-level phylogeny for a study before drawing gain/loss
   conclusions -- this repo already has `nf_phyling` (BUSCO-marker
   multi-locus phylogenomics) for exactly this; a mash-distance tree is not
   a substitute.
2. Replace or augment `polarize_direction()`'s raw-count logic with
   ancestral state reconstruction on that tree (Dollo parsimony as a
   starting point; tools like `Count`, `BadiRate`, or a direct Fitch/Sankoff
   implementation over presence/absence characters are the standard
   approach in the gene-family gain/loss literature).
3. Re-examine `clade_composition`'s role once a real tree exists --
   whether the permutation-null stratification should key off tree clades
   instead of/in addition to the mash-PCoA grouping.
4. Revisit the `genus_vs_ureesii` gain/loss/module findings once this
   lands, since the current numbers are presence/absence-only, not
   phylogenetically resolved.

## Notes

Not implemented in this session -- explicitly flagged by the user as "for
future" work, recorded here so the methodological gap doesn't get lost.
