# Pangenome assumptions and cutoffs register

Every `pangenome_*` parameter in `nextflow.config` (the pangenome cluster-profiling
subworkflow, `pangenome.nf` / `workflows/pangenome_profile.nf`), plus the non-parameter
assumptions the subworkflow depends on, with the evidence behind each one.

**Standard this document holds itself to** (from issue #131): every cutoff should be
defensible in a publication as chosen on the basis of (a) statistics, (b) empirical
result, or (c) theory plus a practical constraint. A parameter that cannot name one of
those three is an open liability and is flagged as `unvalidated` below, not given an
invented rationale.

**Evidence classes used**

| class | meaning |
|---|---|
| `validated-empirically` | Measured on real data, with a result that supports the shipped value specifically. |
| `empirical-negative-result` | Measured on real data, but the result only shows the parameter is *not* the lever — it does not show the shipped value is optimal, or better than nearby values. |
| `theory+practical` | Not swept, but follows from a documented mechanical constraint (a hard relationship to another parameter, a resource/wall-clock limit, a data-shape fact) rather than convention alone. |
| `unvalidated-assumption` | Convention, guess, or carried-over default. No measurement exists. Issue #132 tracks which of these are next in line to be swept. |
| `structural-necessity` | Not a scientific cutoff — a required input, a label, a separator, a reproducibility seed, or a value mechanically forced by another parameter. Listed for completeness, not because it needs evidence. |

Every entry below cites its source: the inline `nextflow.config` comment (already
written by whoever set the default), or one of the dated investigation notes in
`NovInvenio_Investigations/notes/pangenome-method-investigations/` (2026-09-20/21,
read-only reference — not duplicated here beyond what's needed to make this entry
stand on its own), or the relevant ADR under `docs/adr/`.

A value's evidence should be treated as stale once the dataset it is applied to differs
materially in scale or taxon from the one it was measured on (see "Calibrated on" per
entry — everything below was measured on the 529-strain *Coccidioides*
`genus_vs_ureesii` study unless stated otherwise).

---

## 1. Tier-1 clustering — defines what a gene family *is*

This is the foundational step: every downstream count (family count, core/accessory
split, islands, co-occurrence pairs, Leiden modules, gain/loss direction) is computed
on top of whatever these two parameters decide a "family" is.

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_tier1_min_id` | `0.9` | `empirical-negative-result` | Issue #132 priority-1 sweep: 12 `mmseqs easy-cluster` runs over `min_id ∈ {0.70,0.80,0.90,0.95} × cov ∈ {0.5,0.8}` plus cov-mode variants. `rho(accessory ~ N50)` stayed in **-0.511 to -0.528** across all eight `min_id×cov` combinations (a 0.017 range) while family count moved 2x (30,406–60,138 families). The identity value does not touch the fragmentation confound in either direction. | 529-strain Coccidioides, `mmseqs easy-cluster --cluster-reassign` | A different taxon (recombining, paralog-rich genus) could show identity sensitivity this clonal dataset does not — see the sweep note's caveats. A re-sweep is needed before extending this conclusion to a new taxon. | Re-run the min_id×cov grid (`tier1_cluster.tsv` per setting) and check whether `rho(accessory~N50)` moves outside ±0.02 of baseline. |
| `pangenome_tier1_cov` | `0.8` | `empirical-negative-result`, with a documented **structural mechanism working against it** | Same sweep as above — cov also does not move the confound. But independently, the 2026-09-20 rescue characterisation found **71.9% of rescue-pass cells** come from a query rep whose length is <80% of the family it hits, with a sharp frequency cliff at exactly 0.8 — mmseqs' own `-c` cutoff. Contig-terminal proteins (truncated by an assembly boundary) are enriched 3.3x in strain-private (singleton) families (31.4% vs 9.38% background). `--cov-mode 1` is the only setting in the grid that reduces the confound (-0.519 → -0.350, a 33% reduction) but was independently shown to be **disqualified**: 76.0% of its newly-created same-strain copy pairs join two *interior, intact* gene models (not two fragments), and 28.9% merge proteins with disjoint Pfam domain architectures (vs 10.2% for the next-safest alternative). | Same as above | Same as above; also: if a future structural pre-filter (see §2) removes contig-terminal singletons before clustering, re-measure whether the confound then responds to `cov`. | Same grid check; additionally cross-tabulate `--cov-mode 1`'s newly-merged pairs against Pfam domain concordance before ever adopting a looser coverage mode. |
| `pangenome_cluster_backend` | `'mmseqs'` | `validated-empirically` | ADR-0003: mmseqs vs diamond concordance benchmarked on a real 5-strain Coccidioides subset via real SLURM submission — ARI = 0.9403, pair recall 0.8958, pair precision 0.9895 over 43,239 proteins at these subworkflow's own 90%/80% thresholds. Both backends produce comparable but **not interchangeable** results (mmseqs finds ~13% more, smaller families). | 5-strain Coccidioides subset (ADR-0003), not the full 529-strain run | A study needing exact family-count reproducibility across backends, or a much larger/more diverse protein set where the 13% granularity gap could widen. | Re-run `verify_diamond_cluster_ids.py` + an ARI/pair-recall comparison on the target study's own data before treating the two backends as interchangeable for it. |

## 2. Genome-level TBLASTN rescue pass

`pangenome_rescue_enable` exists to recover genes present in the genome but missing
from the annotation (a real `ABSENT → GENOME_ONLY` correction). The 2026-09-20
characterisation established that, as configured before issue #133, it mostly did
something else.

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_rescue_enable` | `true` | `theory+practical`, downgraded by evidence to "**needs `structural_filter` on to be safe**" | Exhaustive characterisation of all 12,877,646 rescuable cells: 77.6% overlap a predicted gene of a *different* family already in the matrix (redundant locus, not new content); only ~1.3% (~163,000 cells) look like plausible genuine annotation dropout. Enabling it without the structural filter would move density 15.6% → 60.3%, mostly by re-asserting genes already in the matrix under another family row. | 529-strain Coccidioides, all 530 per-strain TBLASTN outputs (130M HSP rows), exact/exhaustive not sampled | A well-annotated, non-clonal genome set where truncated/fragmented gene models are rare — the 71.9%-from-the-0.8-coverage-cliff mechanism would not apply. | Compare `presence_matrix.tsv` vs `presence_matrix.rescued.tsv` density before/after; if the jump is large and `pangenome_rescue_structural_filter` is off, re-run the redundant-locus overlap check from the characterisation note before trusting it. |
| `pangenome_rescue_evalue` | `1e-10` | `unvalidated-assumption` | No inline rationale beyond being a conventional strict BLAST e-value; not swept in any note. | — | — | Not yet tracked as an open sweep item; would need its own investigation. |
| `pangenome_rescue_max_target_seqs` | `200` | `theory+practical` | Inline comment: "per-strain DB is tiny; generous headroom over any realistic intra-genome copy number." A judgment call with a stated reason, not a measurement. | — | A genome with unusually high intra-genome copy number for some family (e.g. a large expanded gene family) could exceed 200 real hits and get silently truncated. | No current diagnostic; would need a check that no strain/family pair hits exactly 200 target seqs (a truncation signature). |
| `pangenome_rescue_min_pident` | `90.0` | `empirical-negative-result` | Swept 90→99.9 pident × 80→99 qcov: the redundant-locus/intergenic/paralog class mix stayed **flat at 19-23% intergenic across the entire grid**, because the artifact hits are 100%-identical real DNA. Tightening the threshold only discards cells uniformly; it does not change what fraction of the rescue is spurious. | Same 12.88M-cell characterisation | A taxon with genuinely divergent paralogs (not near-clonal) might show real threshold sensitivity — untested here. | Re-run the (pident, qcov) grid table from the rescue-characterisation note on the target dataset; a flat class mix across the grid means thresholds are not the lever there either. |
| `pangenome_rescue_min_qcov` | `80.0` | `empirical-negative-result` | Same evidence as `min_pident` above — swept jointly. | Same as above | Same as above | Same as above |
| `pangenome_rescue_structural_filter` | `true` (issue #133) | `validated-empirically` | Rejecting a rescue hit whose genomic span overlaps a predicted gene already assigned to a different family removes 77.6% of the 12.88M rescuable cells (the redundant-locus class) directly, using `gene_positions.tsv` the pipeline already produces. This is the one lever in the rescue system that the data actually supports. | Same 12.88M-cell characterisation | Would need re-validation if `gene_positions.tsv`'s gene-model boundaries were produced by a substantially different annotation method. | Re-run the three-way overlap breakdown (different-family / intergenic / same-family) on a new study; the different-family share should be the dominant removed class if the mechanism generalizes. |
| `pangenome_rescue_min_rep_length` | `150` (aa) | `validated-empirically` | 83.2% of the remaining intergenic false positives (after the structural filter) come from family reps <150 aa; median real family rep is 201 aa. A conservative genuine-dropout filter (intergenic AND rep ≥150aa AND not in a hotspot) recovers ~163,000 cells = 1.3% of the rescuable set — the estimated real signal. | Same characterisation, post-structural-filter remainder | A gene family type with a genuinely short canonical length (e.g. small secreted peptides) would need this threshold reconsidered per-taxon. | Check the rep-length distribution of rescued cells vs the study's own family-rep median before assuming 150 aa transfers. |
| `pangenome_rescue_hotspot_window` | `2000` (bp) | `validated-empirically` | 47.5% of the remaining intergenic false positives fall in a 2kb window also hit by ≥3 distinct families — a repeat/mobile-element signature, not independent gene losses. | Same characterisation | A genome with a different repeat-element size distribution (e.g. much larger LTR retrotransposons) could need a wider window. | Re-run the hotspot-window multiplicity check on the target genome's repeat annotation, if available. |
| `pangenome_rescue_hotspot_min_families` | `3` | `unvalidated-assumption` (mechanism validated, exact cutoff value not independently swept) | The *existence* of the hotspot effect is validated (above); "3 or more" as the exact threshold was not itself swept against 2, 4, 5, etc. | Same characterisation | — | Would need a small sensitivity check varying this integer against the same 2kb-window statistic. |

## 3. Strain dedup + clade assignment (Mash → PCoA → k-means)

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_dereplicate` | `true` | `structural-necessity` | Gates whether the frequency/co-occurrence denominator restricts to `is_representative==1`; does not gate any compute today (per-strain processes still run on every strain regardless). | — | — | — |
| `pangenome_mash_threshold` | `0.001` | `validated-empirically` | Single-linkage dedup at 0.001 gives 412 groups from 530 strains (22% collapsed), a "safe operating point." Swept 0.0005→0.005: **0.005 collapses the entire 529-strain set to 2 groups** via single-linkage chaining — the whole thing is on a plateau below 0.002 and dangerous above it. | 529-strain Coccidioides | Any single-linkage collapse rule is vulnerable to chaining as N or diversity grows; the plateau boundary (~0.002) was measured on this dataset only. | Plot groups-vs-threshold (as in the phylogrouping note's table) on the target dataset before raising this value; watch for a cliff to a tiny number of groups. |
| `pangenome_assign_clades` | `true` | `structural-necessity`, with a **known silent failure mode** | On a single-species study there is no dominant axis, but the silhouette sweep always returns *some* k, and `FILL_TAXON_GROUP` writes it into `TaxonGroup` with no warning. Measured: within either Coccidioides species alone, every k=2-10 partition has split-half ARI **0.02–0.26** (irreproducible), vs **1.000** for the real species split. | 529-strain Coccidioides (2 real species) | Any dataset with only one real group (a single-species study). This is flagged in issue #131 as the most dangerous non-parameter assumption (see §7, item 3). | Run the label-concordance / structure-diagnostic checklist below (§3 footnote) before trusting `TaxonGroup` on a new study. |
| `pangenome_clades_k_range` | `'2,10'` | `validated-empirically` (the range itself, not the criterion) | "The sweep is cheap and the range is not the problem; the criterion is." Gap statistic and silhouette agree on k=2 for the real species split; both disagree on within-species subsets (a symptom of no real structure there, not of a bad range). | Same dataset | — | Inspect the full per-k curve, not just the winner, especially for a small or single-group study. |
| `pangenome_clades_k_fixed` | `null` | `structural-necessity` / practical override | Recommended whenever the grouping is already known (e.g. a curated 2-species study) — "a known grouping always beats a discovered one." | — | — | — |
| `pangenome_clades_n_pcoa` | `4` | `validated-empirically` | Increasing it does not help: within-species PCoA axes carry 5-15% variance each with no elbow, so extra axes add noise, not signal. | Same dataset | Only raise if PCoA1 explains <~40% variance *and* there is independent evidence of >4 real groups. | Check PCoA1's variance-explained fraction on the target dataset. |
| **(gap, not yet a parameter)** mash sketch size (`-s`) | hardcoded mash default `1000` | `empirical-negative-result` for the current default; a fix is proposed but **not implemented as a param** | At `s=1000`, all 139,656 ingroup pairs reduce to only 323 distinct distance values (quantization step ~2.7e-5). Re-sketching at `s=100000` gives 16,092 distinct values at <1% extra study runtime (~2.5 min of 306 CPU-hours), and the s=1000 within-species clustering is **not reproducible** against s=100000 (ARI 0.30-0.61). It does not, however, reveal real sub-species structure (PCoA1 only rises 12.5%→14.5%). | 529-strain Coccidioides | — | The phylogrouping note recommends adding a `pangenome_mash_sketch_size` param defaulting to 100000; this has not been implemented yet — tracked as an open item here, not a validated shipped default. |

## 4. Frequency binning cutoffs (core/softcore/shell/cloud)

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_core_cutoff` | `0.95` | `unvalidated-assumption` | Conventional pangenomics values (95%/90%/15% core/softcore/shell splits are widely used in the literature) with no dataset-specific measurement behind them in this repo. Issue #132 priority-2 names this explicitly as open: "measure the frequency distribution's actual shape... is there a natural trough to cut at, or is it a continuum where any cutoff is arbitrary? If the latter, say so... rather than implying the values are principled." | **Not yet measured** | — | This is the entire subject of issue #132 priority 2. Do not treat as validated until that sweep runs. |
| `pangenome_softcore_cutoff` | `0.90` | `unvalidated-assumption` | Same as above. | **Not yet measured** | — | Same as above. |
| `pangenome_shell_cutoff` | `0.15` | `unvalidated-assumption` | Same as above. | **Not yet measured** | — | Same as above. |

## 5. Co-occurrence screen (Fisher / BH-FDR)

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_min_strain_count` | `5` | `unvalidated-assumption` | No sweep exists; a plausible minimum-power convention, not measured. | — | — | Not yet tracked as a sweep target. |
| `pangenome_fdr_alpha` | `0.05` | `structural-necessity` / statistical convention | Standard BH-FDR alpha; not itself a biological claim requiring dataset-specific validation. | — | — | — |
| `pangenome_screen_alpha` | `0.2` | `structural-necessity` | Must stay looser than `fdr_alpha`; `pangenome_cooccurrence.py` clamps this automatically so it cannot be misconfigured tighter than `fdr_alpha`. This is a prefilter threshold, not an independent scientific claim. | — | — | — |
| **(the stratification the screen actually relies on)** `pangenome_pair_class_min_clades` used as the co-occurrence null's stratification depth | `2` (see §6) | `validated-empirically` | Measured FPR at α=0.05 on 300 real accessory-family pairs: unstratified Fisher **0.327**; k=2 stratified (species) **0.050**; k=4/10/20/40/80 all **0.033-0.040**. Finer stratification buys ~0.01 and becomes conservative past k=80. All measurable confounding lives on the species axis; within-species alone the unstratified test is already near-nominal (0.048-0.072). **Do not build a phylogenetically-aware or mixed-model null** — measured benefit over the existing test is approximately zero for substantially more code. | 529-strain Coccidioides, 300 real accessory-pair FPR calibration + within-species sub-check | A dataset where confounding structure is *not* concentrated on one clean axis (e.g. no clear species/clade boundary at all) — untested. | Repeat the 300-random-accessory-pair FPR calibration (unstratified vs k=2 vs finer k) on a new dataset before assuming the existing stratification transfers. |

## 6. Pair classification (physical linkage vs trans)

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_pair_class_k` | `10` | `unvalidated-assumption` | The gene-count window separating "physically linked" from "trans." Issue #132 priority-3 names this an open sweep: "Sweep k and measure how the physical/trans/ambiguous split moves. A classification stable across k=5-20 is defensible; one that swings is not." **Not yet run.** | **Not yet measured** | — | Priority-3 in issue #132; still open as of 2026-09-21. |
| `pangenome_pair_class_physical_threshold` | `0.5` | `unvalidated-assumption` | No sweep exists for this threshold specifically. | **Not yet measured** | — | Bundled with the `pair_class_k` sweep (issue #132 priority 3) but not itself named as swept yet. |
| `pangenome_pair_class_trans_threshold` | `0.05` | `unvalidated-assumption` | Same as above. | **Not yet measured** | — | Same as above. |
| `pangenome_pair_class_min_co_carrying` | `5` | `unvalidated-assumption` | No sweep exists. | — | — | — |
| `pangenome_pair_class_perm_alpha` | `0.05` | `structural-necessity` / statistical convention | Standard permutation-test alpha. | — | — | — |
| `pangenome_pair_class_min_clades` | `2` | `validated-empirically`, with a **documented usage caveat** | Same FPR sweep as §5: k=2 stratification/gating reduces FPR from 0.327 to 0.050 and finer k buys almost nothing (0.033-0.040 at k=10-80). **But** the 2026-09-20 design-decision grilling (Q3) resolved to stop using `min_clades` as a `trans` *filter* when the structure diagnostics say the clade partition is noise (i.e., single-species studies — see §3 and §7 item 3), while still using k=2-style stratification for the co-occurrence null itself. | 529-strain Coccidioides | A study where `pangenome_assign_clades` produced a noise partition (see §3) — using `min_clades≥2` as a hard filter there silently discards real `trans` pairs on the basis of meaningless labels. | Run the structure diagnostics (PCoA1 variance, split-half ARI, between/within ratio — §7 footnote) before trusting `min_clades` as a filter, not just as a stratification depth. |

## 7. Optional captain/mobile-element marker and Pfam/GO annotation

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_captain_hmm` / `pangenome_captain_hmm_name` / `pangenome_pfam_hmm` | `null` | `structural-necessity` | Deliberately optional; omitting both simply means every physically-linked pair classifies as `unexplained_physical` rather than `starship_explained`. Not fungal/Starship-specific in its mechanics. | — | — | — |
| `pangenome_captain_hmm_evalue` | `1e-3` | `unvalidated-assumption` | Conventional HMM e-value cutoff; not swept. | — | — | — |
| `pangenome_island_pfam_hmm` | `null` | `structural-necessity` | Gates the whole island+Pfam step when set; distinct from `pangenome_pfam_hmm` above (used only by the captain-by-name branch). | — | — | — |
| `pangenome_island_min_size` | `2` | `unvalidated-assumption` | Issue #131 names this directly: "what counts as an island at all." No sweep exists. | **Not yet measured** | — | Not yet tracked as an open sweep item beyond the issue-#131 framing. |
| `pangenome_pfam_domain_evalue` | `1e-3` | `unvalidated-assumption` | Conventional Pfam domain e-value; not swept. | — | — | — |
| `pangenome_pfam_chunk_size` | `3000` | `theory+practical` | Issue #112: as one job, `FAMILY_PFAM_SCAN` measured 1h50m-3h50m on the 529-strain study (22,061 reps × 30,134 Pfam-A profiles) and **failed twice at 1h59m against a 2-hour wall-clock cap**, losing the whole scan. `3000` gives ~8 chunks at that scale, landing near the HPCC ~1-1.5h/unit sizing guidance. This is a resource-constraint calibration, not a biological one — it should be raised for a smaller study to avoid per-job submission overhead on chunks that finish in seconds. | 529-strain Coccidioides, real wall-clock failures | A much larger or much smaller background/profile set — the chunk count should be re-derived from the HPCC sizing rule (~1-1.5h/unit), not kept fixed. | Check `FAMILY_PFAM_SCAN` chunk wall-clock times in the trace against the 1-1.5h target; re-tune if chunks are far outside that band. |
| `pangenome_marker_names` / `pangenome_marker_hmm_paths` | `''` | `structural-necessity` | Comma-list feature-gate parameters; no default marker is shipped. | — | — | — |
| `pangenome_marker_evalue` | `1e-5` | `unvalidated-assumption` | Conventional; not swept. | — | — | — |
| `pangenome_pfam2go` | `null` | `structural-necessity` | GO-term annotation is skipped entirely when unset; optional feature. | — | — | — |
| `pangenome_accumulation_permutations` | `20` | `unvalidated-assumption` | Not swept; a conventional small permutation count. Measured cheap (<1s for 20 perms × 530 strains), so raising it costs little — but no measurement establishes 20 is enough for stable rarefaction curves specifically. | — | — | — |
| `pangenome_accumulation_seed` | `0` | `structural-necessity` | Reproducibility seed, not a scientific claim. | — | — | — |

## 8. Trans-module detection (Leiden)

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_leiden_resolution` | `1.0` | **known-wrong at genus scale** (`empirical-negative-result`, actively flagged, not merely unvalidated) | Inline `nextflow.config` comment states directly: "No published pangenome-specific precedent exists for the resolution parameter at genome scale... the default below is not validated for every study; expect to tune it empirically per study." Measured: resolution 1.0 gives **4-5 giant modules** on the 529-strain genus-scale run; resolution 2.0 gives **66 modules at AMI 0.79** against the 1.0 partition. The tier-1 sweep (§1) additionally found Leiden module structure at r=1.0 is *uninformative as a discriminator* between clustering settings (NMI 0.68-0.70, ARI 0.77-0.78 across a 2x change in family count) — the resolution problem masks other signal too. Issue #132 priority-4 (open): needs "a principled selection rule rather than a fixed default — resolution-stability sweep with AMI across seeds, pick the most stable plateau." | 529-strain Coccidioides | Any dataset at meaningfully different scale/diversity — the "too coarse" finding is itself scale-dependent (a smaller study might need a different resolution or none of this problem). | Run a resolution-stability sweep (AMI across seeds, at least 2-3 resolution values) before trusting module counts from `pangenome_leiden_resolution=1.0` at genus scale; issue #132 priority 4 is still open. |
| `pangenome_leiden_seed` | `0` | `structural-necessity` | Reproducibility seed. Confirmed exact reproducibility: a full pipeline re-run gave ARI 0.926 and 97.1% identical module labels against the hand-run baseline (byte-identical md5 on a controlled re-run). | 529-strain Coccidioides | — | — |
| `pangenome_module_min_size` | `2` | `unvalidated-assumption` | No sweep exists for the minimum module size worth summarizing. | — | — | — |

## 9. Report/visualization cutoffs

| parameter | default | evidence class | what set it / evidence | calibrated on | what would invalidate it | how to check |
|---|---|---|---|---|---|---|
| `pangenome_top_islands_min_strains` | `2` | `validated-empirically` | Measured on the 529-strain run: 17,218 of 27,836 located islands (62%) are single-strain, and **all 20 largest islands were single-strain** — an unfiltered size ranking shows only strain-private content, often an assembly/annotation artifact (one 51kb, 69-gene "island" on a single contig of a single strain). Filtering to ≥2 strains removes this. | 529-strain Coccidioides | A study specifically interested in strain-private content (e.g. hunting for a single outbreak strain's unique acquisitions) would want this set to 1 to disable the filter — already supported. | Check the size-vs-strain-count distribution of unfiltered islands on a new study before assuming 2 is the right cutoff there too. |
| `pangenome_viz_top_islands` | `50` | `theory+practical` | Payload-size constraint, measured directly: at top-50, the embedded `island_synteny.html` payload is 0.55MB (0.28MB after a proposed `strain_hap` index-array optimization not yet implemented); top-50 references 298 of 54,421 distinct families. Browsers parse a 5-10MB inline JSON block well under a second from `file://`; the current encoding does not cross 10MB below ~1,500 strains at top-500, or ~4,000 strains at top-500 with the `strain_hap` optimization. | 529-strain Coccidioides, real payload rebuild from `lib/island_synteny.build_payload` | A much larger strain count (thousands) would need this revisited alongside the `strain_hap` encoding change; both are quantified in the storage-and-compute note. | Rebuild the payload at the target `top_islands` value and measure its size directly (`lib/island_synteny.build_payload`) before assuming 50 is still small enough. |

## 10. Structural / labeling parameters (no evidence needed)

These are not scientific cutoffs. Listed for completeness per issue #131's structure,
evidence class `structural-necessity` throughout: `pangenome_samplesheet`,
`pangenome_data_dir`, `pangenome_gff3_dir`, `pangenome_project`, `pangenome_ingroup_label`,
`pangenome_outgroup_label`, `pangenome_id_sep`. None of these encode a scientific claim
that could be right or wrong — they are required inputs, path overrides, or labels.

---

## 11. Non-parameter assumptions (per issue #131, "more dangerous for that reason")

These are not knobs in `nextflow.config` at all — they are assumptions baked into how
the pipeline interprets its own output. Carried over verbatim from issue #131 because
they are still open and this document's job is not to invent resolution where none
exists yet.

| # | assumption | evidence class | status |
|---|---|---|---|
| 1 | **Absence means absence** (an `absent` cell in the presence matrix means the gene is not in the genome). | `empirical-negative-result` — **violated on the flagship run** | `RESCUE_PASS` applied zero cells (a `zstd -dc`-on-symlink bug, fixed 2026-09-20/`e086c78`), so every "absent" call included undetected annotation dropout. Coccidioides' resulting gain:loss ratio was 143:1 with 0% ambiguous, against A. fumigatus' 7.1:1 with 11.6% ambiguous where rescue actually ran (5,871,769 cells applied). The flagship study's presence/absence-dependent numbers (islands, Leiden modules, gain/loss, report tables) need re-measuring once #133's structural rescue criterion is applied and the study is re-run — see §2. |
| 2 | **Accessory content is comparable across strains** (a strain's accessory family count reflects biology, not assembly quality). | `empirical-negative-result` — **violated, quantified** | `rho(accessory ~ N50) = -0.53`, `rho(accessory ~ contigs) = +0.53` after controlling for assembly size (issue #130). Within *C. posadasii*, the 15 most fragmented assemblies carry 283 more accessory families (+11%) than the 15 best, on slightly *less* DNA. Mechanism confirmed directly: strain-private-family proteins sit at a contig terminus 31.4% of the time vs 9.38% background (3.3x enrichment) — a gene broken by a contig boundary yields two partial models that fail clustering and become private families. Tracked by issue #130 (assembly-quality QC emission, no auto-exclusion — advisory only per the 2026-09-20 design decisions). |
| 3 | **Clade labels are real.** | `empirical-negative-result` for the single-species case | On a single-species study, `pangenome_assign_clades` produces a confident-looking but meaningless `TaxonGroup` (split-half ARI 0.02-0.26 for any within-species partition, vs 1.000 for a real species split), and `pangenome_pair_class_min_clades >= 2` then silently filters `trans` pairs on that noise. See §3 and §6. The 2026-09-20 design decisions (Q3) resolved to keep k=2-style stratification for the co-occurrence *null* but stop using `min_clades` as a `trans` *filter* when structure diagnostics say the partition is noise, and to emit diagnostics either way — not yet fully implemented as pipeline output (the structure-diagnostic checklist in §3/§7's phylogrouping note is a manual procedure today, not an automated warning). |
| 4 | **`genome_only` is analytically equivalent to `present`.** | `structural-necessity` as currently implemented, flagged as a real risk | `is_present()` collapses the two; no downstream analysis step can distinguish a TBLASTN rescue hit from an annotated gene call. Design-decision Q1 (2026-09-20 grilling) explicitly agreed these are **not** equivalent ("no they aren't equivalent") — this is recorded as agreed but not yet reflected in a code change that lets downstream steps distinguish the two. |
| 5 | **Gene family = mmseqs cluster at the tier-1 threshold.** | `empirical-negative-result` (same evidence as §1) | The same tier-1 clustering question looked at from the "what is a family" angle rather than the "what parameter value" angle — see §1's `pangenome_tier1_min_id`/`pangenome_tier1_cov` entries. A paralog above the threshold merges into one family; a fragmented gene model below it splits into two. |

---

## 12. What issue #132 has and has not settled, as of 2026-09-21

Per the state-of-the-analysis note, of the ~30 scientific cutoffs in the pangenome
subworkflow, four had a real measurement before this document was written
(`pair_class_min_clades`, `mash_threshold`, `top_islands_min_strains`,
`viz_top_islands`); the tier-1 clustering sweep (priority 1) has since been completed as
a **negative result**; core/shell cutoffs (priority 2), `pair_class_k` and its related
thresholds (priority 3), and Leiden resolution (priority 4, though already known-wrong)
remain **open** as of this writing. Do not treat any priority-2/3/4 parameter above as
more validated than "unvalidated-assumption" until its own sweep note exists and is
cited here.

---

## Maintenance

This document is only useful if it is updated when a value changes or new evidence
arrives. Each entry's "calibrated on" column states the dataset the evidence came from;
a value applied to a materially different dataset (different scale, different taxon,
different clonality) should be treated as unvalidated for that dataset until re-checked,
even if the parameter default itself does not change. When issue #132's remaining
priorities (2-4) land, update §4, §6, and §8 in place rather than adding a duplicate
section.
