# Pin CLUSTER_TIER1's mmseqs branch to AVX2-capable nodes

- **Priority**: medium
- **Status**: open
- **Category**: bug / infra
- **Date**: 2026-09-17
- **Author**: Jason Stajich
- **See**: issue #101, [docs/adr/0003-diamond-tier1-clustering-backend.md](../docs/adr/0003-diamond-tier1-clustering-backend.md)

## What

`modules/pangenome/prefix_and_cluster.nf::CLUSTER_TIER1`'s mmseqs branch has
no SLURM node constraint. Running it on this cluster's AMD "Abu Dhabi" nodes
(no AVX2) crashes `mmseqs easy-cluster` with `SIGILL`, confirmed on a real
`-profile local` smoke test (5-strain *Coccidioides immitis* subset,
2026-09-17). Diamond did not crash on the same node.

## Why

Same root-cause class as the already-documented, already-deferred
2026-07-21 learning (`.living/learnings.md`) about `famsa` SIGILLing on the
same node class — just hitting a different bioconda binary this time.
`modules/mmseqs_cluster.nf` (the novelty/loss pathway's mmseqs module)
already works around this with `clusterOptions = '-C ryzen'`
(`nextflow.config` ~line 163), pinning to AVX2-capable nodes. The pangenome
subworkflow's tier-1 clustering has no equivalent, so a real SLURM run can
silently land on the wrong node class and die, non-deterministically,
depending on scheduler placement.

## How to apply

1. Add the same `-C ryzen` (or an equivalent AVX2-capable constraint) to
   `CLUSTER_TIER1`'s process directives, at minimum when the mmseqs backend
   is selected.
2. Worth checking explicitly whether diamond's binary needs (or benefits
   from) the same pin, or is genuinely more portable across this cluster's
   node fleet as this smoke test suggested — don't assume either way.
3. Once fixed, this unblocks the real mmseqs-vs-diamond cluster-quality
   concordance benchmark tracked in `todo/diamond-tier1-cluster-backend.md`.
