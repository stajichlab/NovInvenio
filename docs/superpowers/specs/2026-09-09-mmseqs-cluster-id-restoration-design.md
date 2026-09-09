# mmseqs cluster.tsv ID restoration (design)

Date: 2026-09-09
Status: approved

## Problem

`mmseqs easy-cluster` silently reduces a `sp|ACCESSION|NAME`/`tr|ACCESSION|NAME`
(UniProt-style) FASTA header to its bare accession in its `*_cluster.tsv`
output (rep/member columns), while its own `*_rep_seq.fasta`/`*_all_seqs.fasta`
outputs from the *same run* keep the full header token. Every other artifact
in this pipeline — `candidates.txt`, `presence_matrix.tsv`,
`novelties.<Short>.tsv`, TBLASTN's `qseqid`, the report payload — treats
`protein_id` as that full header token. Any code that joins a `*_cluster.tsv`
row against `protein_id`-keyed data therefore silently matches nothing, for
every UniProt-sourced study.

Confirmed empirically (a dedicated audit, this session) across real on-disk
data:

| Study | Header convention | `*_cluster.tsv` keys | Affected? |
|---|---|---|---|
| pezizo_set1 | `tr\|F9WWV8\|F9WWV8_ZYMTI ...` | `F9WWV8` (bare) | **Yes** |
| zoosporic_dikarya | `tr\|F4NRB2\|...` | bare accession | **Yes** |
| mushrooms_tremella | `tr\|A8N3V4\|...` | bare accession | **Yes** |
| yeast_filamentous | `sp\|P47149\|...` | bare accession | **Yes** |
| UHM_Akkermansia | MAG/prodigal, no `sp\|`/`tr\|` | full header, verbatim | No — mmseqs's UniProt-parser never triggers on this header shape |
| UHM_Koxytoca | MAG/prodigal, no `sp\|`/`tr\|` | full header, verbatim | No |

Confirmed consequences (all silent — no error, no warning, just empty
results):

- `lib/clusters.py`'s `FamilyIndex.index_of()` — every UniProt-sourced
  study's report ("novelties"/"core"/"losses") shows 0 rows assigned to any
  gene family, even though mmseqs found real multi-member families (verified
  on `pezizo_set1`'s live report: 656 families computed, 0/35,597 rows
  actually placed in one).
- `bin/build_alignment_shards.py` — the report's TBLASTN alignment popup has
  no shards to serve (0 genome shards written against real `pezizo_set1`
  data), because its own `member_to_rep`/`candidate_ids` join has the same
  shape.
- `bin/make_novelties.py`'s own `family_id`/`family_size` columns (in the
  *persisted* `novelties.<Short>.tsv`, independent of the report's
  recomputation above) — same root cause, same silent empty result.
- `bin/summarize_tblastn.py` — **traced and confirmed NOT a correctness
  problem**: its buggy bare-accession-keyed rows are only ever the "no hit"
  rows; a fallback branch in its own code (`protein_hit.setdefault(member,
  ...)` when the buggy lookup misses) happens to create correctly
  full-header-keyed entries for every row that *does* have a real hit —
  verified against `pezizo_set1`'s real `tblastn_summary.tsv`: all 1,053
  "hit" rows are full-header-keyed, only the 3,354 "no hit" rows are
  bare-accession-keyed. Since `make_novelties.py`'s TBLASTN exclusion filter
  and the report's TBLASTN column both only ever query "does this protein_id
  have a hit" (defaulting to "no" when absent), the mis-keyed zero-rows are
  inert. **No candidate has ever been wrongly called a novelty because of
  this bug** — it's a feature-loss bug (family grouping, alignment popups),
  not a scientific-correctness one.

A repo-wide grep found **~16-22 files** in nf_NovInvenio read a
`*_cluster.tsv` and build a `protein_id`-keyed map from it — **zero** of them
normalize. NII already solved the equivalent problem on its own side
(`NovInvenio_Investigations/lib/uniprot_ids.py`'s `bare_accession()`, used by
`bin/merge_uniprot_annotations.py`/`go_enrichment.py`/`domain_enrichment.py`),
but nf_NovInvenio has no equivalent, and the two repos don't share it.

Separately audited and confirmed **clean, not in scope for this fix**:
TBLASTN's `sseqid` and genome/scaffold/DNA accession handling (BLAST+ does
not do mmseqs's UniProt-header reparsing; a real genome's FASTA header, GFF3
seqid column, and TBLASTN `sseqid` all agreed end to end in the case
checked). HMM builds (`hmmbuild`/`hmmsearch`, `novelty_discovery`'s
family-profile path) introduce no *new* mismatch — HMMER names a profile
after whatever ID the cluster.tsv rep column already holds and passes it
through domtblout verbatim, so fixing the cluster.tsv fixes the HMM path too,
with no HMM-specific change needed.

## Decision

Fix it once, immediately after mmseqs produces each `*_cluster.tsv`, before
any of the ~16-22 existing consumers ever reads it — not by rewriting FASTA
headers before mmseqs (which would propagate a *new* ID discontinuity into
TBLASTN's query chain, one of the two domains confirmed clean today), and not
by patching every consumer individually (which is exactly the shape of
problem that produced this bug — a 23rd future consumer will forget, with
nothing to catch it).

### 1. Exactly two insertion points

Only two Nextflow processes ever invoke `mmseqs easy-cluster` in this
pipeline:

- `modules/mmseqs_cluster.nf`'s `MMSEQS_CLUSTER` — used only by
  `--cluster_tool pairwise` (`workflows/cluster.nf`'s `CLUSTER`/`LOSS_CLUSTER`,
  both novelty and loss directions share this one process definition).
- `modules/mmseqs_family_cluster.nf`'s `MMSEQS_FAMILY_CLUSTER` — used by
  **both** `--cluster_tool mmseqs` (`workflows/profile_search.nf`) and
  `--cluster_tool novelty_discovery` (`workflows/novelty_discovery.nf` calls
  it directly). Since both call the same shared module, fixing this one file
  covers both `cluster_tool` values in one edit.

`PROFILE_CANDIDATE_CLUSTERS` (`modules/profile_candidate_clusters.nf`, used
by both `mmseqs` and `novelty_discovery` to derive a candidate-scoped
cluster.tsv from `MMSEQS_FAMILY_CLUSTER`'s output) does **not** call mmseqs
itself — it only filters an already-produced `family_cluster_tsv` via
`bin/candidate_families.py`. Once `MMSEQS_FAMILY_CLUSTER`'s own output is
corrected, everything downstream of it — including
`PROFILE_CANDIDATE_CLUSTERS` — inherits the fix automatically, no separate
edit needed there.

So: **2 files touched** (`modules/mmseqs_cluster.nf`,
`modules/mmseqs_family_cluster.nf`), **0 files touched** among the ~16-22
existing consumers, covering all three `--cluster_tool` values
(`pairwise`, `mmseqs`, `novelty_discovery`).

### 2. The corrector script

New script, name chosen to name the actual problem it fixes (mmseqs's own
header-collapsing behavior, not a generic "cluster fixer"):

**`bin/restore_mmseqs_cluster_ids.py`**

```
Usage:
  restore_mmseqs_cluster_ids.py \
      --input-fasta <the FASTA mmseqs clustered> \
      --cluster-tsv <mmseqs easy-cluster's *_cluster.tsv, corrected in place>
```

Algorithm:

1. Parse `--input-fasta`'s headers (the exact FASTA fed to `mmseqs
   easy-cluster` — i.e. the process's own input, always available since
   it's already a process input). For each record, take the full
   first-whitespace-delimited header token (`full_id`) as the ground truth
   `protein_id`, and derive its bare form (`bare_id`) via the *same*
   `sp|`/`tr|`-stripping rule as NII's `lib/uniprot_ids.py::bare_accession()`
   — ported into nf_NovInvenio (new `lib/uniprot_ids.py`, a direct copy, so
   both repos apply the identical rule and never drift apart). A non-UniProt
   header (MAG/prodigal-style) has `bare_id == full_id` — the function is a
   no-op pass-through for those, exactly like NII's version.
2. Build `bare_to_full: {bare_id: full_id}`. **Fail loudly** (exit 1, matching
   this repo's established fail-loud convention) if two distinct `full_id`
   values collapse to the same `bare_id` — a genuine accession collision in
   the source FASTA, not something to silently pick one of.
3. Read `--cluster-tsv` line by line (`rep\tmember` per mmseqs's own format).
   For each column value `x`: if `x in bare_to_full`, replace it with
   `bare_to_full[x]` (the collapsed case — this is the actual fix). If `x`
   already equals some `full_id` exactly (the MAG/prodigal case, or any
   value mmseqs happened to leave untouched), leave it unchanged — the
   replacement is idempotent by construction, so this branch is a no-op that
   still passes the value through. If `x` matches **neither** — an mmseqs
   output ID this script cannot reconcile against the input FASTA at all —
   **fail loudly** rather than silently emit a row nothing downstream can
   use.
4. Write the corrected file. Overwrite the same `--cluster-tsv` path (Step 1
   of "Nextflow wiring" below explains why this needs no output-declaration
   change).

This is a plain, dependency-free script (stdlib only, matching this repo's
existing `bin/*.py` style) — a FASTA header pass plus a two-column TSV
rewrite, no new dependency.

### 3. Nextflow wiring — in-process, no new process

Append the corrector as one more shell command at the end of each of the two
processes' existing `script`/`shell` block, run against the exact filename
`mmseqs easy-cluster` just wrote, **before** the process's `output:` block
takes effect — so `publishDir` copies out the already-corrected file. Example
for `MMSEQS_CLUSTER` (`modules/mmseqs_cluster.nf`):

```groovy
    mmseqs easy-cluster \
        !{candidates_fa} \
        !{out_prefix} \
        ${tmpdir}/mmseqs_$$ \
        --threads !{task.cpus} \
        --min-seq-id 0.3 \
        -c 0.8 \
        --cov-mode 0

    restore_mmseqs_cluster_ids.py \
        --input-fasta !{candidates_fa} \
        --cluster-tsv !{out_prefix}_cluster.tsv
```

Same pattern for `MMSEQS_FAMILY_CLUSTER` against `families_cluster.tsv`
(input FASTA: `seed_fa`). **No change to either process's `output:`
declaration, no change to any of the ~16-22 consumer files** — they keep
reading the exact same filename, now correctly keyed. `*_rep_seq.fasta`/
`*_all_seqs.fasta` need no correction (confirmed unaffected — mmseqs already
keeps full headers there).

### 4. Backfill for already-computed studies

Because the script's inputs (`candidates.fa`/`loss_candidates.fa`/
`seed_all.faa` and the raw `*_cluster.tsv`) are exactly what's already
sitting in `results/<study>/` for any study already run, this can be applied
**retroactively, standalone, with no pipeline rerun**:

```bash
bin/restore_mmseqs_cluster_ids.py \
    --input-fasta results/pezizo_set1/candidates.fa \
    --cluster-tsv results/pezizo_set1/clusters/clusters_cluster.tsv
```

— repeated per study, per direction (`candidates.fa`/`clusters_cluster.tsv`
and `loss_candidates.fa`/`loss_clusters_cluster.tsv`), for every already-run
UniProt-sourced study. This is the same manual/standalone-script pattern
already established this session (Task 6 of the xrefs plan) — no new
convention introduced.

## Documentation

Per explicit instruction: this is a real, non-obvious behavior of a
third-party tool (mmseqs), not a self-explanatory bug — document it where a
future reader would actually look:

- **`METHOD_DESCRIPTION.md`** (repo root) — this file already has an
  algorithm-level description of clustering for both the `pairwise` (line
  ~107: "Candidates are extracted, clustered with `mmseqs easy-cluster`...")
  and `mmseqs` (line ~141: "`mmseqs cluster` (cascaded, ...) clusters the
  whole seed group...") pathways. Add a short paragraph immediately after
  each clustering-step description (or one shared paragraph if a natural
  shared location exists, e.g. right after both pathways are introduced)
  explaining: mmseqs collapses `sp|`/`tr|`-style headers to a bare accession
  in `*_cluster.tsv` specifically (not in `*_rep_seq.fasta`), why that would
  silently break every family/alignment-shard feature without correction,
  and that `restore_mmseqs_cluster_ids.py` runs immediately after every
  `mmseqs easy-cluster` invocation to fix this before anything else reads
  the file.
- **`README.md`**'s "Two analysis approaches: pairwise vs family-profile
  clustering" section (line ~46) — one sentence pointing at
  `METHOD_DESCRIPTION.md`'s fuller explanation, so a reader who never opens
  that file still learns this exists.
- **`bin/restore_mmseqs_cluster_ids.py`**'s own module docstring — the full
  technical explanation (mmseqs's specific behavior, the empirical evidence
  table above in condensed form, the algorithm), following this repo's
  existing convention of putting the detailed "why" in the script that
  implements it (matches `bin/build_alignment_shards.py`,
  `bin/extract_dat_annotations.py`'s existing docstring style).
- **`modules/mmseqs_cluster.nf`** and **`modules/mmseqs_family_cluster.nf`**
  — a one-line comment at the corrector-script call site pointing back at
  the script's docstring, not re-explaining it.

## Testing

- `bin/restore_mmseqs_cluster_ids.py`: unit tests against synthetic fixtures
  covering (a) a UniProt-style FASTA + a bare-accession-collapsed cluster.tsv
  → corrected output matches the full headers; (b) a MAG/prodigal-style
  FASTA + an already-full-header cluster.tsv → output unchanged
  (idempotent no-op); (c) two distinct headers colliding to the same bare
  accession → the script exits non-zero; (d) a cluster.tsv ID matching
  neither a bare nor full form → the script exits non-zero.
- Nextflow-level: extend whichever existing test fixture already exercises
  `MMSEQS_CLUSTER`/`MMSEQS_FAMILY_CLUSTER` (if any) to assert the *published*
  `*_cluster.tsv` is now full-header-keyed for a UniProt-style input fixture.
- Manual/integration check: run the backfill command (Section 4) against
  `pezizo_set1`'s real `clusters_cluster.tsv`, then re-run
  `report_data.py`'s `FamilyIndex` (or just regenerate the report) and
  confirm rows are now actually assigned to families (currently 0/35,597);
  separately re-run `bin/build_alignment_shards.py` against the corrected
  file and confirm it now writes real shards (currently 0).

## Out of scope

- TBLASTN `sseqid` / genome-scaffold-DNA accession handling — confirmed
  clean, no change.
- `bin/summarize_tblastn.py`'s own bare-accession-keyed "no hit" rows —
  confirmed inert (Problem section above); once `*_cluster.tsv` itself is
  corrected at the source, `summarize_tblastn.py`'s *input* is already
  right, so this stops being an issue as a side effect, but it did not
  independently need a targeted fix.
- Migrating existing per-consumer code in the ~16-22 files to use a shared
  normalizer — deliberately not done (Decision section explains why: the
  fix at the source makes this unnecessary, and doing both would be
  redundant defensive layering with no added correctness).
