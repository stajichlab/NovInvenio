# Pin CLUSTER_TIER1's mmseqs branch to AVX2-capable nodes

- **Priority**: n/a
- **Status**: wont-do (not a bug — already fixed, see correction below)
- **Category**: bug / infra
- **Date**: 2026-09-17
- **Author**: Jason Stajich
- **See**: issue #101 (closed, not planned), [docs/adr/0003-diamond-tier1-clustering-backend.md](../docs/adr/0003-diamond-tier1-clustering-backend.md)

## Correction (2026-09-17, same day)

This was a false alarm. `conf/ucr_hpcc_slurm.config` already has:

```
withName: '.*CLUSTER_TIER1' {
    clusterOptions = '-C ryzen|broadwell|cascade'
}
```

added in commit `3de5252` (the commit that introduced the whole pangenome
subworkflow), well before the smoke test below ever ran. The `SIGILL` this
item was filed over happened because that smoke test used `-profile local`
only, without `-c conf/ucr_hpcc_slurm.config` — `-profile local` just sets
`process.executor = 'local'` and never applies any `clusterOptions`, which
only take effect under `-profile slurm` with that config file included. A
real SLURM run already avoids the non-AVX2 nodes for this process, for
both backends. No code change needed. Issue #101 closed as not-a-bug.

Diamond not crashing under `-profile local` on a non-AVX2 node is still a
genuine, minor point in its favor for local/interactive dev-node smoke
testing (no need to remember `-c conf/ucr_hpcc_slurm.config` just to poke
at the pipeline), but it is not evidence of a production gap.

---

## Original (incorrect) item, kept for the record

`modules/pangenome/prefix_and_cluster.nf::CLUSTER_TIER1`'s mmseqs branch has
no SLURM node constraint. Running it on this cluster's AMD "Abu Dhabi" nodes
(no AVX2) crashes `mmseqs easy-cluster` with `SIGILL`, confirmed on a real
`-profile local` smoke test (5-strain *Coccidioides immitis* subset,
2026-09-17). Diamond did not crash on the same node.
