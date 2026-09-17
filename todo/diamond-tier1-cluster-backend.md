# Real-data concordance benchmark for the diamond tier-1 clustering backend

- **Priority**: medium
- **Status**: open
- **Category**: validation
- **Date**: 2026-09-16
- **Author**: Jason Stajich
- **See**: [docs/adr/0003-diamond-tier1-clustering-backend.md](../docs/adr/0003-diamond-tier1-clustering-backend.md), issue #98

## What

`--pangenome_cluster_backend diamond` is no longer hard-blocked — the
ID-fidelity gap is closed by `bin/verify_diamond_cluster_ids.py`
(unit-tested against synthetic adversarial headers, wired into
`CLUSTER_TIER1` in `modules/pangenome/prefix_and_cluster.nf`). What remains
is a real-data concordance benchmark: run both backends' tier-1 clustering
on the same real, Short-prefixed multi-strain proteome set and compare
(ARI, pairwise recovery/precision — same metrics as the existing
2026-09-08 mmseqs-vs-diamond benchmark in `.living/decisions.md`) at the
actual tier-1 identity/coverage regime (`--min_seq_id 0.9`/`--cov 0.8` →
`--approx-id 90`/`--member-cover 80`).

## Why

No `diamond`/`mmseqs` binaries were available in the environment used to
build the fidelity check, so it has never run against a real `diamond
cluster` invocation. The existing benchmark in `.living/decisions.md` used
raw UniProt headers from a different pipeline path and doesn't cover the
pangenome's `Short<id_sep>orig_id` convention or this identity/coverage
regime — it's suggestive, not sufficient, evidence for this use.

## How to apply

1. Pick a real multi-strain proteome set already prefixed via
   `bin/pangenome_prefix_fasta.py` (or run that step fresh).
2. Cluster it with both `pangenome_cluster_backend.py mmseqs-tier1` and
   `diamond-tier1` at the pipeline's real defaults.
3. Confirm `verify_diamond_cluster_ids.py` passes cleanly on the diamond
   output (if it doesn't, that's the actual finding — stop and investigate,
   don't patch around it).
4. Compare cluster membership (ARI, pair recovery/precision) between the two
   backends' output family sets.
5. Record the result in `.living/decisions.md` next to the existing
   2026-09-08 entry, and update this ADR's Status line once done.
