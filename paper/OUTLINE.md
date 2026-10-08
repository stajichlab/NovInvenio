# Paper outline: lineage-specific gene discovery and accessory-island analysis with NovInvenio

Status: working outline, 2026-10-08. Not a manuscript. Every result below is marked with where
its evidence is and whether it must be rerun before publication.

Status tags: **[have]** result exists in a file named here. **[rerun]** exists, but from an
older pipeline version or with a known flaw, so rerun before use. **[todo]** not done.

## 1. Working title options

1. NovInvenio: finding lineage-specific genes and accessory gene islands across hundreds of fungal genomes.
2. Pairwise and pangenome screens for gene novelty, with evidence of what each call rests on.
3. From presence/absence to islands: a transparent pipeline for lineage-specific and accessory genes.

## 2. One-sentence claim (to be confirmed by the results)

NovInvenio finds genes that are present in a defined ingroup and absent from an outgroup, and
finds accessory gene islands in large strain collections, and it shows for every call which
evidence supports it and which filters shaped it.

## 3. Sections

### Abstract
Problem, two screens (pairwise novelty; pangenome islands), validation, application, availability.
Write last.

### 1. Introduction
- Why lineage-specific genes matter (adaptation, host range, pathogenicity islands, Starships).
- Why presence/absence calls are fragile: annotation gaps, paralogs, assembly fragmentation,
  clonal strains inflating counts.
- What exists (orthology tools, pangenome tools; cite in a literature pass) and what is missing
  (calls that show their evidence; one framework for novelty and for islands).
- Contributions (list, matched to sections 3 and 4).

### 2. Methods (drafts exist: see `methods_pairwise.md`, `methods_pangenome.md`, flow charts)
2.1 Overview and information flow (`flow_pairwise.md`, `flow_pangenome.md`; one combined figure).
2.2 Pairwise novelty pipeline: search, filters 1 to 3, presence matrix, candidate rule,
    candidate clustering, TBLASTN, annotation, outgroup-signal evidence, reports.
2.3 Pangenome pipeline: dereplication, clustering, presence matrix, rescue, bins, co-occurrence
    test, pair classification, islands, locus view, clinker panels, modules, Pfam enrichment,
    openness, QC.
2.4 Validation design: curated controls, BUSCO controls, simulation, known-element recovery.
2.5 Software, data availability, parameter table, compute (see the cost note).

### 3. Results: pairwise screen
| # | Result | Evidence | Status |
|---|---|---|---|
| 3.1 | Paralog competition rescues HEX-1-like orthologs and removes A7UWR3-like cross-hits | controls in `configs/controls/pezizo_set1.controls.csv`; `.living/decisions.md`; `tests/test_build_presence_matrix.py` | [have] needs one clean rerun table |
| 3.2 | Curated controls: pairwise recovered 6/6 resolved pezizo_set1 positives and 0/8 false positives on sordariales_shallow BUSCO negatives | `notes/cluster-vs-pairwise/README.md` (NII) | [have] but positives were screened with the same diamond search (circular), and only about 7 independent genes |
| 3.3 | Search sensitivity changes cost by 1 to 2.8x and changes candidate sets | `notes/diamond-sensitivity/` (NII) | [have] |
| 3.4 | Outgroup evidence explains why a candidate has outgroup hits (filter-removed paralog hits; domain-only coverage) | #208 case A0A090CJI0_PODAN; run `sordariales_shallow_ev208` | [rerun] job was still running when this was written |
| 3.5 | Candidate counts per clade, annotation, TBLASTN support | published reports under `docs/fungi/` | [rerun] older pipeline versions |
| 3.6 | Cost scaling and the taxa bound | `notes/pairwise-cost-model/` (NII) | [have] retrospective only; a designed scaling run is [todo] |

### 4. Results: pangenome screen
| # | Result | Evidence | Status |
|---|---|---|---|
| 4.1 | Species-level grouping agrees with labels (527 of 529 concordant; the two exceptions are mislabels) | `notes/pangenome-method-investigations/2026-09-21-state-of-the-analysis.md` | [have] |
| 4.2 | The co-occurrence null is calibrated: false-positive rate 0.327 unstratified, 0.050 with the shipped 2-stratum test | same note | [have], Coccidioides only; A. fumigatus [todo] |
| 4.3 | Rescue pass: a silent failure was found and fixed; 12.9 million rescuable cells are 74% redundant loci | same note | [have]; rerun on current pipeline [todo] |
| 4.4 | Assembly fragmentation inflates accessory counts (rho -0.53 with N50) | same note | [have]; check singleton batch effect [todo] |
| 4.5 | Island detection needs about 21 to 30 or more representative genomes | `studies/fungi/Afumigatus_pangenome/analysis/island_genome_count/ISLAND_GENOME_COUNT.md` (NII) | [have] subsets reuse 293-genome families |
| 4.6 | Known Starships are recovered: all 7 high-confidence Af293 Starships hit by islands; 21.3% of placed island genes lie in Starships vs 5.2% of non-core genes; 19 of 19 scorable Starships match by presence pattern above a permutation null | same file (full_v070, 11,775 islands) | [have] but with the 8,226 vs 11,775 island gap unexplained until the gap check finishes |
| 4.7 | Island locus view and clinker panels on A. fumigatus and Coccidioides | relaunch runs `full_v071`, `immitis_in_posadasii_out_v3` | [todo] runs in progress |
| 4.8 | Reciprocal immitis vs posadasii comparison | Coccidioides runs | [rerun] |
| 4.9 | Other groups (Zymoseptoria, Neurospora, Cryptococcus, Pyriculariaceae) | BFD counts in project memory | [todo] candidates only |

### 5. Discussion
- What the two screens can and cannot say. Absence means no qualifying hit, not proof of absence.
- Where evidence disagrees (annotation gap vs deletion vs true absence): the DNA presence check.
- Limits: positive-control set is small and partly circular; Cocci-only calibration; one example strain per island
  view; clone handling; mixed clade labels.
- Compute and the taxa bound.

### 6. Availability
Code (GitHub), container, parameters, example data, the NovInvenio_Investigations repository for study configs and published pages.

## 4. Planned figures and tables
1. Information flow of both pipelines (from `flow_*.md`).
2. Pairwise filters: example orthologs kept and removed (HEX-1 vs eIF5A; A7UWR3).
3. Candidate evidence panel: protein hits, filtered hits, TBLASTN coverage (the Outgroup signal).
4. Pangenome frequency classes per group; per-genome class composition.
5. Co-occurrence null calibration (FPR by stratification).
6. Island genome-count sweep (islands and Starship hits vs N).
7. Island synteny and locus view for one island, with the clinker panel and DNA check.
8. Cost vs number of taxa.
Tables: parameter table (from the methods files); control sets; per-clade candidate counts; studies and data provenance.

## 5. What must happen before writing results

See `OPEN_ISSUES_FOR_PAPER.md`. The main blockers are statistical definitions (the `trans` call, the
openness rule, mixed denominators) and a validation set that is independent of the search it tests.
