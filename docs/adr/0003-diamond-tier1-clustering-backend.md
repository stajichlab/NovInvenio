# diamond as a validated tier-1 clustering backend for the pangenome workflow

## Status

Accepted, implemented and fully validated — 2026-09-16 (design), 2026-09-16
(implementation: fidelity check + guard relaxation), 2026-09-17 (real
end-to-end pipeline validation), 2026-09-17 (real-data mmseqs-vs-diamond
concordance benchmark). See #98. Extends the pangenome subworkflow
(`pangenome.nf`, `workflows/` under it); does not touch the novelty/loss
pathway (ADR-0001/0002) or `workflows/cluster.nf`.

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

**2026-09-17 update — real concordance benchmark, done.** Ran both backends
via real SLURM submission (`-profile slurm -c conf/ucr_hpcc_slurm.config`)
against the identical 5-strain subset: 37/37 processes, 0 failures, for
both. diamond produced 10,468 families; mmseqs produced 11,830 (mmseqs
splits into more, smaller clusters at these thresholds). Comparing the two
`tier1_cluster.tsv` outputs over the same 43,239 proteins:
**ARI = 0.9403**, pair recall (mmseqs pairs recovered by diamond) =
**0.8958**, pair precision (diamond pairs confirmed by mmseqs) = **0.9895**
— substantially stronger concordance than the 2026-09-08 novelty/loss-pathway
benchmark (ARI 0.79, 72%/87%), because this uses the pangenome subworkflow's
own tight, matched 90%-identity/80%-coverage tier-1 thresholds for both
tools rather than each tool's own defaults at a looser 30%-identity target.
Recorded in `.living/decisions.md`'s 2026-09-17 entry. Every question this
ADR originally raised is now closed: ID-fidelity, pipeline integration, and
cluster-quality-vs-mmseqs. `--pangenome_cluster_backend diamond` is a
validated, real alternative — not merely experimental — though mmseqs
remains the default and the two are not numerically interchangeable
(mmseqs finds ~13% more, smaller families at these settings).

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
4. **Real-binary header-fidelity smoke test, real end-to-end pipeline run,
   and real-data concordance benchmark — all done.** A real `diamond
   cluster` run against a small adversarial FASTA (including the
   `Afum|sp|O74225|YCF1_SCHPO` multi-pipe case) preserved every header
   verbatim. Then, 2026-09-17, a full `pangenome.nf
   --pangenome_cluster_backend diamond -profile local` run against a real
   5-strain *Coccidioides immitis* subset completed 37/37 processes with
   0 failures, producing 10,467 real gene families. An initial attempt at a
   same-data mmseqs comparison hit an unrelated AVX2/SIGILL crash on the
   `-profile local` test node — a testing-methodology artifact (that
   profile never applies `conf/ucr_hpcc_slurm.config`'s existing
   `-C ryzen|broadwell|cascade` node pin for `CLUSTER_TIER1`), not a real
   gap; issue #101 was filed then closed same-day as not-a-bug once this
   was found. Re-run via real SLURM submission (`-profile slurm -c
   conf/ucr_hpcc_slurm.config`), both backends completed 37/37, 0 failures,
   on the identical 5-strain subset: **ARI = 0.9403**, pair recall = 0.8958,
   pair precision = 0.9895 over the same 43,239 proteins — substantially
   stronger concordance than the 2026-09-08 novelty/loss-pathway benchmark
   (ARI 0.79, 72%/87%) at this subworkflow's own tighter, matched tier-1
   thresholds. Full detail in `.living/decisions.md`'s 2026-09-17 entry.
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
  them to diamond's `--approx-id`/`--member-cover` as percentages). Answered
  by the item-4 benchmark above: comparable but not identical (ARI 0.94;
  mmseqs finds ~13% more, smaller families at these settings).
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
  rubber-stamp, and it has now been exercised against synthetic adversarial
  data and two real end-to-end SLURM runs on real biological data (0
  failures each). `--pangenome_cluster_backend diamond` is fully validated:
  pipeline mechanics, and cluster granularity at ARI 0.94 against mmseqs at
  the same configured thresholds. mmseqs remains the default; a study
  sensitive to exact family granularity should still treat the two as
  measurably different (mmseqs: ~13% more, smaller families), not
  interchangeable.
- A `-profile local` smoke test initially looked like it had surfaced a real
  infra gap (mmseqs SIGILLing with no AVX2 node-pinning). It hadn't — the
  pin already exists in `conf/ucr_hpcc_slurm.config` and simply isn't
  applied by `-profile local`. Worth remembering for future local-profile
  testing on this cluster: a `-profile local` pass exercises pipeline logic
  but not the SLURM-specific node/resource constraints a real run relies on.
