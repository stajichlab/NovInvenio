# Real-data concordance benchmark for the diamond tier-1 clustering backend

**Update 2026-09-16**: a real `diamond cluster` smoke test (pixi env, diamond
v2.2.0.180) against a small adversarial FASTA confirmed header-fidelity
holds for a real invocation, including the multi-pipe Short-prefixed
UniProt-style case (`Afum|sp|O74225|YCF1_SCHPO`) — `verify_diamond_cluster_ids.py`
passed. That FASTA was 5 small near-identical toy sequences, not real
biological data, so it validated ID fidelity only, not clustering quality.

**Update 2026-09-17**: ran the real thing. `pangenome.nf
--pangenome_cluster_backend diamond -profile local` against a real 5-strain
*Coccidioides immitis* subset completed 37/37 processes, 0 failures, 10,467
real gene families, normal core/shell/singleton frequency distribution.
Attempted the same-subset mmseqs comparison for a direct concordance number
and it crashed instead — `mmseqs easy-cluster` SIGILLed on the test node
(no AVX2), unrelated to this work (see issue #101, same root cause as the
2026-07-21 famsa/AVX2 learning). Diamond did not crash on the same node.
The concordance benchmark below is now blocked on #101 (need an
AVX2-capable node to run mmseqs at all for the comparison), not on tooling.

- **Priority**: medium
- **Status**: blocked (on issue #101)
- **Category**: validation
- **Date**: 2026-09-16
- **Author**: Jason Stajich
- **See**: [docs/adr/0003-diamond-tier1-clustering-backend.md](../docs/adr/0003-diamond-tier1-clustering-backend.md), issues #98, #101

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
