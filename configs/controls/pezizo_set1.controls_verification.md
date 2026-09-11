# pezizo_set1 positive-control verification (2026-09-10)

Companion to `pezizo_set1.controls.csv`. Records the empirical check behind each row —
real `diamond blastp --very-sensitive` search of each Neurospora crassa (Ncra) anchor
protein against all 11 `pezizo_set1` proteome databases (`results/pezizo_set1/search_cache/*.dmnd`),
not a literature claim taken on faith. Gene list originally supplied by the user
(Jason Stajich) from Neurospora genetics literature; NCU IDs resolved via
`config_support/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv` (NII repo) where the
UniProt `GN=` tag didn't already carry the common name.

Ingroup = Amega, Ncra, Afum, Ztri, Cimm (Pezizomycotina). Outgroup = Nirr, CneoH99, Ccin,
Mcir, Scer, Spom. `ingroup_min_frac` default is 0.75 (>=4/5).

## Result table (best e-value per species, `-` = no hit even at e<=1)

| gene | Amega | Ncra | Afum | Ztri | Cimm | Nirr | CneoH99 | Ccin | Mcir | Scer | Spom |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hex-1 | 4.6e-52 | 4.7e-121 | 1.6e-70 | 8.9e-83 | 6.8e-74 | 1.2e-07 | 3.2e-10 | 1.5e-11 | 3.4e-12 | 1.7e-13 | 7.6e-12 |
| wsc | 3.6e-79 | 5.0e-213 | 3.7e-88 | 1.0e-88 | 4.2e-92 | - | 3.5e-03 | 4.1e-48 | 4.4e-46 | 2.3e-06 | 3.2e-03 |
| lah (NCU02793) | 1.1e-59 | 0.0 | 4.1e-111 | 1.6e-114 | 4.3e-117 | - | - | - | - | - | - |
| so/soft | 1.1e-300 | 0.0 | 0.0 | 0.0 | 0.0 | 6.5e-04 | 4.4e-41 | 3.5e-54 | 2.4e-58 | 1.9e-41 | - |
| ada-1 | 3.7e-83 | 0.0 | 1.9e-73 | 1.6e-62 | 4.9e-77 | - | 1.8e-01 | - | 5.7e-02 | - | - |
| ham-5 | 4.3e-134 | 0.0 | 5.7e-241 | 3.0e-238 | 1.6e-208 | - | - | 4.8e-01 | - | 2.9e-01 | - |
| ham-8 | 6.3e-71 | 0.0 | 8.0e-148 | 1.0e-134 | 2.3e-133 | - | - | - | - | - | - |
| spa-1 | 7.4e-12 | 7.8e-303 | 7.4e-17 | 8.9e-49 | 2.5e-37 | - | - | - | - | - | - |
| spa-9 | - | 0.0 | 2.0e-42 | 5.4e-01 | 1.6e-42 | - | - | - | - | - | - |
| ham-11 | - | 0.0 | 1.3e-06 | - | 3.9e-06 | - | - | - | - | - | - |
| spa-18 | 6.8e-14 | 0.0 | 8.1e-38 | - | 8.5e-29 | 6.8e-08 | - | - | - | - | - |

## Verdicts

| gene | ingroup frac | outgroup contamination | verdict |
|---|---|---|---|
| hex-1 | 5/5 | weak hit in all 6 (eIF-5A paralog cross-reactivity) | in `controls.csv`, flagged: only trust once the paralog-competition filter (post `--very-sensitive` self-search fix, 2026-09-10) is confirmed disqualifying these |
| lah | 5/5 | 0/6 | in `controls.csv`, clean |
| ada-1 | 5/5 | noise-level in 2/6 (e-01, e-02) | in `controls.csv`, clean |
| ham-5 | 5/5 | noise-level in 2/6 (e-01) | in `controls.csv`, clean |
| ham-8 | 5/5 | 0/6 | in `controls.csv`, clean |
| spa-1 | 5/5 | 0/6 | in `controls.csv`, clean |
| spa-18 | 4/5 (missing Ztri) | weak hit in Nirr (6.8e-08) | in `controls.csv`, flagged (same paralog-risk class as hex-1) |
| spa-9 | 3/5 (no hit Amega, marginal Ztri) | 0/6 | **excluded** — fails `ingroup_min_frac=0.75` for this exact species set regardless of broader Pezizomycotina literature support |
| ham-11 | 3/5, and weak (e-06) where present | 0/6 | **excluded**, same reason |
| wsc | 5/5 | **strong** hits in Ccin (4.1e-48) and Mcir (4.4e-46) | **excluded from controls.csv** — contradicts clean lineage restriction; likely a broader WSC-domain family issue, not a paralog the competition filter would catch |
| so/soft | 5/5 | **strong** hits in CneoH99/Ccin/Mcir/Scer (e-41 to e-58) | **excluded from controls.csv** — contradicts literature claim of Basidiomycota/yeast absence; needs its own investigation before reuse anywhere |

## 2026-09-10 update: checked against the real post-fix rerun

After the `DIAMOND_SELF --very-sensitive` fix and a fresh `pezizo_set1` run
(`candidates.txt`, `presence_matrix.tsv`, `self_hits/Ncra.paralog_cutoffs.tsv`):

- **hex-1**: confirmed fixed. Self-search now pairs hex-1 with eIF-5A (2.84e-06); the
  Mcir presence column flipped 1->0; hex-1 is back in `candidates.txt`. Moved to
  "ready to use" in `pezizo_set1.controls.csv`.
- **lah, ada-1, ham-5, ham-8, spa-1**: all clean and all correctly present in
  `candidates.txt` on the real run. Ready to use.
- **spa-18**: re-flagged as **excluded**, not just "watch." Nirr's presence-matrix
  column is `1` on the real run -- the paralog-competition filter did NOT disqualify
  that weak cross-hit, so spa-18 does not make it into `candidates.txt` at all. Needs
  its own root-cause before reuse (real distant ortholog in Nirr vs. an uncaught
  paralog artifact).
- **so/soft**: the real `presence_matrix.tsv` shows it as clean (absent in all 6
  outgroup) -- but this **contradicts** the standalone `--very-sensitive` diamond
  check above, which found strong hits (e-41 to e-58) in 4 outgroup species. The
  pipeline's main pairwise search (`DIAMOND_SEARCH`, `modules/diamond.nf`) still runs
  at diamond's *default* sensitivity, unlike the self-search (fixed to
  `--very-sensitive`). so/soft's clean appearance in real output is most likely a
  **false negative for outgroup presence caused by search insensitivity**, not
  genuine lineage-restriction -- concrete evidence for
  `todo/diamond-very-sensitive-main-search.md`, added there.

## Negative controls (2026-09-10)

Added 10 BUSCO-anchored negative controls (`NEG_BUSCO01`-`10`) to `pezizo_set1.controls.csv`,
`anchor_type: busco`, `expected_call: core` (must never be called novel).

**Used `fungi_odb12`, not `ascomycota_odb12`** (the user asked about ascomycota_odb12
specifically): `ascomycota_odb12`'s single-copy guarantee only holds within Ascomycota.
Three of `pezizo_set1`'s own outgroup proteomes are NOT Ascomycota (Ccin, CneoH99 =
Basidiomycota; Mcir = Mucoromycota) -- an ascomycota_odb12 BUSCO gene has no guaranteed
presence in those, so it wouldn't be a valid "must be universally present" negative for
this specific 11-species panel. `fungi_odb12` (pan-fungal) is the correct choice for a
negative control that must hold across the whole ingroup+outgroup set, and this study
already has real `fungi_odb12` BUSCO runs for all 11 species (`NovInvenio/busco_pezizo5/`,
BUSCO 6.0.0) -- no new BUSCO run was needed.

Computed the intersection of "Complete" (single-copy) BUSCO IDs across all 11 species'
`run_fungi_odb12/full_table.tsv`: **315 of 1122** BUSCOs are universally single-copy in
this exact panel (a real, non-trivial floor below 1122 -- Nirr and Scer recover
noticeably fewer BUSCOs overall, 713 and 780 of 1122, plausibly real reduced-genome
biology for those two lineages, not just annotation completeness). 10 were picked as an
initial representative sample (functional categories: translation, DNA replication/
repair, transcription, proteasome, tRNA synthetase, vacuolar ATPase, general metabolism)
-- easy to extend from the same 315-ID intersection later if a larger negative set is
wanted.

`anchor_type: busco` was used specifically because BUSCO's own per-species `Sequence`
column ID (e.g. `NCU07852-t26_1-p1`) does NOT match this study's UniProt-sourced
`Ncra.pep.fa` protein IDs (`sp|...|..._NEUCR`) -- the `busco_pezizo5/` BUSCO run was
against a different (older, BFD/FGSC-era) Ncra gene-model set. `anchor_type: busco` sidesteps
this ID mismatch entirely (the scorer is expected to resolve BUSCO ID -> family/candidate
independent of which specific protein-ID scheme a run used), which is exactly why the
controls-file template recommends it for negatives.

## First real score_controls.py run (2026-09-10) -- against pezizo_set1_cluster (mmseqs pathway)

`bin/score_controls.py` already existed (committed 2026-09-05, `66be3a2` era) but was
built specifically for `--cluster_tool mmseqs`, not the pairwise pathway `pezizo_set1`
itself uses -- so it was run against `pezizo_set1_cluster` (the mmseqs comparison study,
same species/data) instead. Command + full per-control output:
`results/pezizo_set1_cluster/controls_scoring/pezizo_set1.controls_scored.tsv` (+
`.summary.tsv`, + `Ncra.busco_map.partial.tsv`, the busco_id->protein_id crosswalk built
for this run -- only 5/10 BUSCO negatives resolved, see below).

**Result: recall 0.400 (2/5 resolved), fp_rate 0.000 (0/5 resolved).**

Per-control breakdown:

| control | resolved family (rep) | call | outcome |
|---|---|---|---|
| hex-1 | Amega rep (A0ACF5BXR3) | not_novel | **miss** |
| lah | -- | unresolved | protein not in any profiled family |
| ada-1 | Amega rep (A0ACF5CBF7) | not_novel | **miss** |
| ham-5 | itself (V5IN79) | not_novel | **miss** |
| ham-8 | Cimm rep (J3K3I3) | novel | hit |
| spa-1 | Cimm rep (J3K8R1) | novel | hit |
| NEG_BUSCO02/03/04/05/08 | (resolved) | not_novel | tn (all 5 correct) |
| NEG_BUSCO01/06/07/09/10 | -- | unresolved | NCU locus number didn't match any `GN=` tag in `Ncra.pep.fa` -- needs a real crosswalk (e.g. `config_support/modelorgs/Ncra_self_id_crosswalk.tsv`), not yet chased down |

**Interpretation, not yet root-caused further**:
- **hex-1 miss is expected, not a bug in the scorer**: the mmseqs/family-profile pathway
  has no self-vs-self paralog-competition filter at all (ADR-0002's own stated
  limitation), so it has no mechanism to disqualify the eIF-5A cross-hit the way the
  pairwise pathway's fixed `DIAMOND_SELF` now does. This is real, first-hand evidence of
  the tradeoff flagged earlier in this investigation: HMM-pathway sensitivity helps
  *outgroup absence* detection but has no defense against *paralog* cross-reactivity --
  a structurally different failure mode from pairwise's, not a strictly better one.
- **lah unresolved** likely means it didn't cluster into a >=2-member family at all
  (mmseqs `--cluster_tool mmseqs` drops true singletons per ADR-0002 decision 6) -- not
  confirmed by inspecting the cluster TSV directly, just the most likely explanation.
- **ada-1/ham-5 misses**: not investigated further. Could be real outgroup homology the
  family HMM detects that pairwise diamond (even `--very-sensitive`) didn't, or a
  clustering over-merge pulling in an unrelated family member. Open question.
- Only 5/10 BUSCO negatives resolved because the BUSCO run's own protein IDs
  (`NCU07852-t26_1-p1` style, from an older/different Ncra gene-model set) don't all
  match this UniProt-sourced proteome's `GN=NCU#####` tags 1:1 -- 5 NCU numbers had no
  `GN=` match at all in `Ncra.pep.fa`. Not chased further this session.

## Root-caused: ada-1 / ham-5 misses, and the 5 remaining BUSCO ID mismatches (2026-09-10)

**BUSCO ID mismatches**: chased down 9 of the original 10 via the real NCU locus tag ->
gene-symbol lookup (`Ncra.gff3`) then a `GN=<symbol>` search in `Ncra.pep.fa` (5 matched
directly on the bare NCU number; `rpn-6`/`leu-6`/`vma-4`/`mbf-1` needed the gene-symbol
hop; `mbf-1` specifically was found by matching the BUSCO's own description text,
"multiprotein bridging factor", not the locus tag). The 10th (`100036at4751`, generic
"eukaryotic translation initiation factor 3", locus `NCU03876` has no assigned gene
name and this genome has ~13 different eIF3-subunit proteins) was deliberately left
unresolved rather than guessed -- `results/pezizo_set1_cluster/controls_scoring/Ncra.busco_map.partial.tsv`
now has 9/10.

**ada-1 and ham-5 misses, root-caused via the raw `family_hmmsearch/*.domtblout` files**:
both are real, and both trace to the SAME mechanism -- a **promiscuous shared domain**
fragmenting into two separate weak partial-domain hits whose *merged* span clears
`hmm_presence_min_residues=100` even though neither individual domain match is
biologically meaningful:

- **ada-1**'s family rep's best Mcir hit (`S2JBT0_MUCC1`, "BZIP domain-containing
  protein") has two HMM domain hits at hmm-coords 160-197 (37 aa, not independently
  significant, domain i-Evalue 2.5e+03) and 260-334 (74 aa, i-Evalue 1.3e-06). Neither
  alone reaches 100 residues or 50% of the 647-aa query length, but merged (111 aa
  total) they clear the `min_residues=100` floor. Nearly every one of Mcir's dozen-plus
  weak hits to this family is independently annotated "BZIP domain-containing protein"
  -- ada-1 almost certainly carries a small, common bZIP-like motif shared by many
  unrelated Mcir transcription factors, not a real ortholog.
- **ham-5**'s family rep's best Mcir hit (`S2JCL3_MUCC1`, "Anaphase-promoting complex
  subunit 4 WD40 domain-containing protein") shows the identical pattern: two weak
  partial domain hits (110-230, 603-717 in HMM coordinates) to an unrelated WD40-repeat
  protein, again likely a shared common motif, not orthology.

This is a real, concrete false-positive-enabling side effect of the
`hmm_presence_min_residues=100` alternative-to-coverage floor (added, per
`nextflow.config`'s own comment, specifically to stop penalizing long multi-domain
proteins for having only one conserved region) -- it can also let two *disjoint, weak,
promiscuous-domain* partial hits merge past the same floor with no real full-protein
homology behind them. This is exactly the kind of `fp_rate` evidence
`todo/validate-hmm-presence-coverage-broader-sweep.md` has been missing (it only had a
curation-free false-*absence* proxy via BUSCO recovery, never a real false-*presence*
example) -- worth folding into that sweep's evidence before any default change.

## Extended score_controls.py for the pairwise pathway (2026-09-10)

`bin/score_controls.py` was family-mode-only (`--cluster-tsv`/`--families` required,
hardcoding the mmseqs family-expansion step). Extended it to also support the pairwise
pathway directly: if `--cluster-tsv`/`--families` are both omitted, every protein is
treated as its own one-member "family" (`identity_membership()`), so the exact same
`family_presence_vector()`/`family_call()` scoring path is reused unchanged for both
producers -- consistent with ADR-0002's "same pivot, free UI" design (`fasta` anchors
are the one exception: not yet supported in pairwise mode, no family HMM db exists to
hmmsearch against; reported unresolved).

Ran it against `pezizo_set1` itself (the pairwise run all these controls were originally
verified against by hand): **recall 1.000 (6/6), fp_rate 0.000 (0/9)** -- exactly matches
the manual `candidates.txt`/`presence_matrix.tsv` spot-checks earlier in this file, now
automated and reusable for any future rerun or param-sweep grid point. Output:
`results/pezizo_set1/controls_scoring/pezizo_set1.controls_scored.tsv` (+ `.summary.tsv`).
One BUSCO negative (`100036at4751`, eIF3 subunit, `NCU03876`) is still unresolved --
the one deliberately-not-guessed busco_map entry from the mismatch-chasing above.

## Fix implemented and confirmed: `min_domain_evalue` (2026-09-10)

Implemented in `lib/family_presence.py` (`parse_domtblout`/`family_presence_by_proteome`
gained a `min_domain_evalue` parameter), wired through `bin/profile_to_matrix.py`,
`bin/novelty_presence_matrix.py`, `bin/novelty_screen.py`, `modules/profile_presence_matrix.nf`,
`workflows/novelty_discovery.nf`, `workflows/novelty_screen.nf`, and `nextflow.config`
(new `hmm_presence_domain_evalue`, default `null` = disabled, exact prior behaviour).
15 unit tests added/passing in `tests/test_family_presence.py` (including two that pin
the exact ada-1 domtblout data as a regression case), full suite 403 passed / 3 skipped
(one real Nextflow `val()`-can't-carry-`null` bug caught and fixed by the existing
end-to-end integration test -- see commit).

**Confirmed on real data**: regenerated `pezizo_set1_cluster`'s presence matrix directly
from the already-computed `family_hmmsearch/*.domtblout` files (no full rerun needed)
with `--min-domain-evalue 0.01`, added on top of the shipped
`--evalue 1e-3 --min-coverage 0.5 --min-covered-residues 100`, then re-scored:

**recall improved 0.400 -> 0.600 (3/5 resolved), fp_rate unchanged at 0.000 (0/9)**.

- **ada-1: miss -> hit**, exactly as predicted. Its Mcir cross-hit was two fragments
  (37aa, i-Evalue 2.5e+03 noise; 74aa, i-Evalue 1.3e-06) that only cleared
  `min_residues=100` when merged. The fix excludes the noise fragment, leaving 74aa --
  correctly below both gates now.
- **ham-5: still a miss** -- a genuinely different, harder case. Its Mcir cross-hit
  (`S2JCL3_MUCC1`, "Anaphase-promoting complex subunit 4 WD40 domain-containing
  protein") has domain 1 (120aa, i-Evalue 0.13 -- now correctly excluded) AND domain 2
  (114aa, i-Evalue 6.1e-05 -- INDEPENDENTLY significant and independently clears
  `min_residues=100` alone, no merging involved). No per-domain-significance filter can
  distinguish "genuinely significant 114aa domain match" from "shared promiscuous fold,
  not real orthology" -- that would need a different mechanism entirely (e.g. requiring
  the match also cover a substantial fraction of the TARGET protein, not just the query
  HMM -- a reciprocal-coverage check). Flagged as a follow-up, not solved here.

## Which pathways this touches (asked 2026-09-10)

`lib/family_presence.py` (and therefore this whole fix) is used ONLY by the HMM/
family-profile pathways: `--cluster_tool mmseqs` (`bin/profile_to_matrix.py`) and the
two-phase `--cluster_tool novelty_discovery` (`bin/novelty_presence_matrix.py`,
`bin/novelty_screen.py`) -- both consolidated onto this one shared library by the
2026-09-03 fix (decision #11). **The pairwise pathway
(`bin/build_presence_matrix.py`, `pezizo_set1`'s own default `--cluster_tool pairwise`)
is architecturally immune**: it has no HMM, no per-target domain-coverage/span concept
at all -- a pairwise hit is a flat E-value plus the paralog-competition filter, with
nothing analogous to "merge multiple domain spans" for a promiscuous domain to exploit.
Its own failure mode is different (paralog cross-reactivity, the HEX-1/eIF-5A case
fixed earlier this session), not this one.

## Caveats

- `-` (no hit at e<=1 under `--very-sensitive`) is strong but not absolute evidence of
  true absence — same caveat as the hex-1/eIF-5A investigation this session started from.
- This check used the anchor's *own* best hit per proteome only; it does not run the
  pipeline's actual paralog-competition filter, so "outgroup contamination" here means
  "diamond finds *some* homology," not "the pipeline would call it present" — for hex-1
  and spa-18 specifically, whether they're usable controls depends on that filter
  actually firing correctly on a real run.
- `wsc`/`so` results were a surprise relative to general Pezizomycotina-vs-outgroup
  literature statements; not deeply investigated further (e.g. whether the outgroup hit
  is to a real ortholog, a shared domain, or a same-genome paralog) — flagged, not solved.
