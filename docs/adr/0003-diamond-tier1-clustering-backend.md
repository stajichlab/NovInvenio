# diamond as a validated tier-1 clustering backend for the pangenome workflow

## Status

Accepted, implemented and validated end-to-end — 2026-09-16 (design),
2026-09-16 (implementation: fidelity check + guard relaxation),
2026-09-17 (real end-to-end pipeline validation). See #98. Extends the
pangenome subworkflow (`pangenome.nf`, `workflows/` under it); does not
touch the novelty/loss pathway (ADR-0001/0002) or `workflows/cluster.nf`.

A real `diamond cluster` invocation (pixi environment, diamond v2.2.0.180)
against a small adversarial FASTA confirmed the header-preservation claim
this ADR was built on: the multi-pipe, Short-prefixed UniProt-style header
`Afum|sp|O74225|YCF1_SCHPO` came through diamond's `*_cluster.tsv` unchanged,
and `bin/verify_diamond_cluster_ids.py` passed against that real output.

**2026-09-17 update — real end-to-end run.** Ran `pangenome.nf
--pangenome_cluster_backend diamond -profile local` against a real 5-strain
*Coccidioides immitis* subset (`NovInvenio_Investigations`' coccidioides
pangenome study data). Completed 37/37 processes, 0 failures. `CLUSTER_TIER1`
produced 10,467 real gene families from 10,549 proteins; the frequency
distribution was exactly what a real pangenome should look like (7,017 core
@ freq 1.0, 1,698 near-singleton @ freq 0.2, a normal shell tail in
between); `verify_diamond_cluster_ids.py` printed its pass confirmation
inside the actual pipeline run, not just a standalone test. This is now
validated against real biological data through real Nextflow orchestration,
not just unit tests and a standalone CLI smoke test.

A same-subset mmseqs comparison run was attempted for a direct concordance
check and **could not complete** — unrelated to this ADR: `mmseqs
easy-cluster` itself crashed with `SIGILL` on the test node, which lacks
AVX2 (confirmed via `/proc/cpuinfo`). **This turned out to be a
testing-methodology artifact, not a pipeline gap**: `conf/ucr_hpcc_slurm.config`
already pins `CLUSTER_TIER1` to AVX2-capable nodes
(`clusterOptions = '-C ryzen|broadwell|cascade'`, added in commit `3de5252`,
predating this smoke test) — the comparison run used `-profile local` only,
which never applies `clusterOptions` (those only take effect under
`-profile slurm` with `-c conf/ucr_hpcc_slurm.config` included). A real
SLURM run already avoids the non-AVX2 nodes for this process, for both
backends. (Filed as issue #101, then closed as not-a-bug the same day once
this was found.) Diamond not crashing under `-profile local` on a non-AVX2
node is still a minor point in its favor for local/interactive dev-node
smoke testing, but it is not evidence of a production gap.

**Still open**: a real mmseqs-vs-diamond cluster-quality comparison (ARI,
pair-recovery) on the *same* real data — unblocked now that the AVX2
question is resolved, it just needs to run via `-profile slurm -c
conf/ucr_hpcc_slurm.config` rather than `-profile local` — tracked in
`todo/diamond-tier1-cluster-backend.md`. The ID-fidelity and pipeline
integration questions this ADR originally raised are now closed; only the
clustering-quality-vs-mmseqs question remains. Treat
`--pangenome_cluster_backend diamond` as validated-but-uncompared: safe to
run, but its cluster granularity relative to mmseqs at the same identity
target is not yet independently confirmed.

## Context

Two independent conversations converged on the same gap.

1. The user asked whether correlated gene gain/loss (a proxy for functional
   interaction) is better detected from raw cluster membership than from the
   `--cluster_tool novelty_discovery` family-HMM pathway, which has repeatedly
   shown itself to be over-sensitive to the outgroup and its own coverage/E-value
   floors (see `lib/family_presence.py`'s module docstring; the 83%-rejection
   finding on HMMs >600aa). Cluster-membership presence/absence is cheaper and
   decouples the correlation question from strict per-family significance
   calibration.
2. That analysis **already exists** for the pangenome workflow:
   `bin/pangenome_cooccurrence.py` computes Fisher/BH-FDR co-occurrence between
   shell/cloud gene families, polarizes direction (gain vs. loss) against the
   outgroup, and corrects for shared ancestry with an exact within-clade
   stratified permutation test (`exact_stratified_pvalue`) keyed on each
   strain's `TaxonGroup`. It consumes `presence_matrix.tsv` /
   `frequency_table.tsv`, which are themselves built from whichever tool
   produced `tier1_cluster.tsv` — **the co-occurrence analysis is already
   cluster-backend-agnostic.** Nothing there needed to change for a second
   backend.

So the actual gap was narrower than "build a new correlation analysis": it
was "make the pangenome workflow's tier-1 clustering step safe to run with a
second backend." `modules/pangenome/prefix_and_cluster.nf::CLUSTER_TIER1`
already dispatched on `params.pangenome_cluster_backend` ('mmseqs' |
'diamond'), and `bin/pangenome_cluster_backend.py` already had a
`diamond-tier1` subcommand (`diamond cluster -d <fasta> -o <out>_cluster.tsv
--approx-id 90 --member-cover 80`). But `pangenome.nf` hard-errored if
`diamond` was requested — it was disabled at the entrypoint, not merely
untested — because the diamond path had no equivalent of
`bin/restore_mmseqs_cluster_ids.py`, the safety net that runs after mmseqs
`easy-cluster` (wired into `modules/mmseqs_cluster.nf` and
`modules/mmseqs_family_cluster.nf`) to undo mmseqs' header-collapsing for
UniProt/NCBI-style deflines (`sp|ACC|NAME` → `ACC`, etc. — see
`lib/fasta.py::mmseqs_id()`).

### Why this was a smaller problem than it looked

The novelty/loss pathway clusters raw proteome FASTAs, so a UniProt-sourced
species' headers (`sp|ACC|NAME`) hit mmseqs' recognized-prefix table directly
and get collapsed to `ACC` — that's what `restore_mmseqs_cluster_ids.py`
exists to fix.

The pangenome subworkflow clusters differently: `bin/pangenome_prefix_fasta.py`
rewrites every header to `Short<id_sep>orig_id` (default `Short|orig_id`,
`params.pangenome_id_sep = '|'`) **before** `CONCAT_PROTEOMES`/`CLUSTER_TIER1`
ever run. `lib/fasta.py::mmseqs_id()` inspects only `fields[0]` (the token
before the first `|`) to decide whether to collapse a header. Once every
header's first field is a `Short` code (`Ncra`, `Afum`, ...), it can never
equal one of mmseqs' recognized prefixes (`sp`, `tr`, `gb`, `ref`, ...) unless
a Short code happens to collide with one verbatim (see Open Questions —
not closed by this change). This means mmseqs' collapsing is **structurally
defused** by the prefix step for the common case. diamond, per
`.living/learnings.md`'s 2026-09-08 entry, already preserves the full
first-token header verbatim and does *not* do UniProt-style accession
extraction the way mmseqs does — but that finding was measured against
**raw** UniProt headers on a different pipeline path (`.living/decisions.md`'s
"Stick with mmseqs2" entry, `NovInvenio_Investigations`' 51,803-protein
`seed_all.faa`), not against `Short|orig_id`-prefixed, potentially
multi-pipe headers (`Short|sp|ACC|NAME` when `orig_id` itself contains
pipes — `pangenome_prefix_fasta.py` does not escape or reject `id_sep`
occurring inside `orig_id`).

## Decision

1. **Verify, don't assume, header fidelity** for both backends under the
   pangenome ID convention — done. `bin/verify_diamond_cluster_ids.py` checks
   every id in a `tier1_cluster.tsv` against the full FASTA header set,
   failing loud (not silently) on an unmatched id, and distinguishing that
   from the more alarming case of an id that only matches after stripping an
   mmseqs-style `sp|`/`tr|`/... wrapper (which would mean diamond had, against
   its documented behavior, collapsed a header the way mmseqs does).
   `tests/test_verify_diamond_cluster_ids.py` exercises adversarial headers,
   including the specific multi-pipe, Short-prefixed case this ADR flagged
   as never having been tested (`Afum|sp|O74225|YCF1_SCHPO`). This is
   verify-only — it never rewrites the file, unlike
   `restore_mmseqs_cluster_ids.py` — because nothing in diamond's documented
   behavior or this test suite gave a reason to build a restoration path;
   if a real run's fidelity check ever fails, that is a stop-and-investigate
   signal, not something to auto-correct.
2. **Wired into `CLUSTER_TIER1`** (`modules/pangenome/prefix_and_cluster.nf`):
   the diamond branch runs `verify_diamond_cluster_ids.py` against the raw
   `tier1_cluster.tsv` immediately after `diamond cluster` produces it, before
   the existing `cut`/`awk` rep-seq reconstruction touches it and before
   anything downstream reads it.
3. **`pangenome.nf`'s hard-error guard relaxed** to a doc comment pointing at
   this ADR and the new safety net, rather than an unconditional `error`.
   `--pangenome_cluster_backend diamond` now runs.
4. **Real-binary header-fidelity smoke test, then a real end-to-end pipeline
   run — both done; mmseqs-vs-diamond cluster-quality concordance — still
   open.** A real `diamond cluster` run against a small adversarial FASTA
   (including the `Afum|sp|O74225|YCF1_SCHPO` multi-pipe case) preserved
   every header verbatim. Then, 2026-09-17, a full `pangenome.nf
   --pangenome_cluster_backend diamond -profile local` run against a real
   5-strain *Coccidioides immitis* subset completed 37/37 processes with
   0 failures, producing 10,467 real gene families with a normal
   core/shell/singleton frequency distribution — this is no longer a toy
   synthetic-data check, it is a real biological dataset through real
   Nextflow orchestration. What's still missing is a same-data mmseqs
   comparison: an attempt at one hit an unrelated AVX2/SIGILL crash on the
   `-profile local` test node, which turned out to be a testing-methodology
   artifact (that profile never applies `conf/ucr_hpcc_slurm.config`'s
   existing `-C ryzen|broadwell|cascade` node pin for `CLUSTER_TIER1`) — not
   a real gap; issue #101 was filed then closed same-day as not-a-bug once
   this was found. The existing mmseqs-vs-diamond benchmark in
   `.living/decisions.md` (ARI 0.79, 72% pair recovery at 87% precision)
   used raw UniProt headers on a different pipeline path and predates this
   ID convention — still suggestive, not load-bearing evidence, for this
   pangenome study. The real-data concordance benchmark is the one item
   left open, and just needs a `-profile slurm` run rather than
   `-profile local`; see `todo/diamond-tier1-cluster-backend.md`.
5. **No changes made to `bin/pangenome_cooccurrence.py`, `lib/pangenome_matrix.py`,
   or the report layer** — confirmed backend-agnostic by inspection (they
   read `tier1_cluster.tsv`/`presence_matrix.tsv` structurally, never caring
   which tool produced them), and no test exercising them needed updating.

## Open Questions

- Does the pangenome path still need `restore_mmseqs_cluster_ids.py` at all,
  given the prefix step already defuses mmseqs' recognized-prefix collapsing
  for the common case? Left wired as-is (cheap, fail-loud, already tested) —
  not re-examined in this change.
- Should a Short code be validated at config-parse time against mmseqs' and
  diamond's recognized-prefix lists (reject `sp`, `tr`, `gb`, ... as a Short),
  closing the one remaining theoretical collision path? Not implemented here
  — no evidence it has ever occurred, and it is a config-validation concern
  orthogonal to the clustering-safety-net fix this ADR scopes.
- Tier-1 identity/coverage thresholds (`--pangenome_tier1_min_id`,
  `--pangenome_tier1_cov`) are shared params across both backends
  (`modules/pangenome/prefix_and_cluster.nf`'s `diamond-tier1` branch maps
  them to diamond's `--approx-id`/`--member-cover` as percentages). Whether
  this mapping gives comparable cluster granularity between the two tools is
  exactly what the still-open real-data benchmark (item 4 above) needs to
  answer — not assumed here.
- Not in scope here: gene-tree-based coevolution (evolutionary rate
  covariation / MirrorTree-style correlated substitution rate across
  orthogroups) is a distinct signal from gain/loss co-occurrence — catches
  functional coupling even where neither family is gained or lost — and
  would be a separate analysis mode built on per-family MSA+tree, not on
  this presence/absence machinery at all.

## Consequences

- Users get a real choice of tier-1 clustering engine for the pangenome
  workflow. The correlated-gain/loss analysis (co-occurrence, direction
  polarization, clade-stratified significance) is available on
  diamond-built families with zero new analysis code, because that analysis
  was already built backend-agnostic.
- The fidelity check is a genuine runtime safety net (fail loud on the first
  bad id, same contract as `restore_mmseqs_cluster_ids.py`), not a
  rubber-stamp, and it has now been exercised both against synthetic
  adversarial data and against a real end-to-end run on real biological
  data (10,467 families, 0 failures). Anyone running
  `--pangenome_cluster_backend diamond` on a real study can trust the
  pipeline mechanics; what remains unverified is only whether diamond's
  cluster granularity at the configured identity/coverage target matches
  mmseqs' closely enough for a given study's purposes — that comparison
  just needs to run under `-profile slurm` rather than `-profile local`.
- A `-profile local` smoke test initially looked like it had surfaced a real
  infra gap (mmseqs SIGILLing with no AVX2 node-pinning). It hadn't — the
  pin already exists in `conf/ucr_hpcc_slurm.config` and simply isn't
  applied by `-profile local`. Worth remembering for future local-profile
  testing on this cluster: a `-profile local` pass exercises pipeline logic
  but not the SLURM-specific node/resource constraints a real run relies on.
