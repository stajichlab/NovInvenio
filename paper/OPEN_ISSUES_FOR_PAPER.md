# Open issues that affect what the paper can claim

Collected 2026-10-08 from three independent code reviews (pairwise, pangenome, island and clinker
views) and from two methods drafts written from the code. Status key:
**fixed** = merged or in an open PR. **open** = known, not changed. **decision** = needs the
maintainer's choice because it changes results.
"Confirmed" means a second reader or a data check agreed with the reviewer. Items marked
"review only" were not rechecked line by line.

## A. Pairwise pipeline

| ID | Issue | Why it matters for the paper | Status |
|---|---|---|---|
| P1 | `--evalue` did not reach filter 1 (script default 1e-5 always applied) | Methods would state a parameter that did not act | fixed in PR #224; same value at default |
| P2 | Docs said TBLASTN removes candidates; the pipeline passes `--skip_tblastn_filter` | 306 of 568 candidates in sordariales_shallow have an outgroup TBLASTN hit (review figure) | fixed in PR #224 (docs, reports) |
| P3 | Core 95% means 100% for fewer than 20 proteomes | Core claims on small panels | documented on the page (PR #224) |
| P4 | The paralog comes from a self-search at E up to 100; 3,064 of 7,080 Neurcras "paralogs" have E above 1e-5 | A noise paralog can drop a weak true hit | open; test against the HEX-1 and A7UWR3 controls before changing |
| P5 | Outgroup signal mixes diamond best-HSP coverage with TBLASTN genome-wide union coverage | "broad" can be inflated by joining several loci | open; split into two values |
| P6 | A protein with no surviving hit has no matrix row; one ingroup proteome gives zero novelties by construction | Edge case to state in methods | documented (PR #224) |
| P7 | Three copies of filters 2 and 3 with different defaults (`singleton_presence` default scope is `proteome`) | Reproducibility | open; mostly removed by #219 |
| P8 | Positive controls were screened with the same very-sensitive diamond search; about 7 independent genes; `hmm_presence_domain_evalue` was tuned on two of the scoring controls | Validation is partly circular and leaks | open; needs independent or simulated controls |
| P9 | TBLASTN hits and coverage come from the cluster representative at 30% identity and are copied to members | Members can differ from the representative | documented; open as a limitation |
| P10 | Diamond runs with its default of 25 target sequences per query | Filter 2 can be lenient in large families (suspected) | open, unmeasured |

## B. Pangenome pipeline

| ID | Issue | Why it matters for the paper | Status |
|---|---|---|---|
| G1 | The `trans` call uses an uncorrected p-value. `permutation_p` is an exact within-clade hypergeometric test run only on BH-FDR survivors, and `trans` needs `permutation_p < 0.05`. In full_v070 60% of the 564,245 `trans` pairs have 0.01 <= p < 0.05 (review figure). The value is written with 4 decimals, so p below 5e-5 prints as 0.0000 | The main statistical claim for trans modules | **decision**: write full precision, add a BH q over the tested set, and gate `trans` on it. Two reviewers agree |
| G2 | The TaxonGroup column mixes label schemes (252 Barber, 10 DAPC, 31 Mash in full_v070). The null shuffles within these strata and `unknown` counts as a clade | Stratification depth is not what the text implies | open; enforce one scheme or record it |
| G3 | Denominators differ between steps. Bins and co-occurrence use ingroup representatives; pair classification, islands and the accumulation curve use every strain, including near-duplicates and the outgroup | Counts in one table cannot be compared to another | open; see also G13 |
| G4 | The openness rule is `is_open = gamma < 1` for pangenome size P = k N^gamma. Under that model almost any sublinear curve is "open" | The openness verdict says nothing | rule fixed to gamma > 0 in PR #226 (the verdict is still weak: any growing curve has gamma > 0); fitting on ingroup representatives only is open |
| G5 | "Accessory" has three definitions (not core; not core or soft core; shell plus cloud) | Methods must use one | open |
| G6 | Pfam is scanned only on shell plus cloud representatives, so island tables miss domains of singleton, outgroup-only and non-representative members | Domain enrichment | open |
| G7 | Islands are deduplicated on the exact member set, so one locus can appear as several islands; "62% single-strain islands" is partly this | Island counts | open; merge by overlap (the locus view already groups by containment) |
| G8 | The "rescue redundancy" diagnostic fires after the structural filter has already removed the rows, and tells the reader not to enable rescue | Misleading banner | open |
| G9 | Assembly-quality QC tests accessory counts against N50 and contigs, not singletons. In isl45_v1, 6 of 9 singleton outliers are near-complete assemblies (reviewer) | A batch effect would go unnoticed | open |
| G10 | The calibration (false-positive rate 0.327 to 0.050) is from Coccidioides only | Generality | open; repeat on A. fumigatus |
| G11 | The test called `permutation_p` is exact, not a permutation test; the Leiden edge weight uses the unstratified FDR, not the stratified p | Naming and consistency | open |
| G12 | Clinker slices drew zero genes for UniProt-ID proteomes (Af293 and both reference outgroups) | The reference genome was blank in every synteny panel | fixed in PR #222 |
| G13 | The island count depends on gene positions of non-representative strains: with positions of all 295 strains full_v070 gives 11,775 islands, with only the 123 representatives it gives 8,226 (same scripts, same co-occurrence output; checked 2026-10-08, `outputs/gap_check/`) | Island counts depend on clone strains acting as extra evidence | **decision**: restrict `n_co_carrying` and `linkage_fraction` to representatives, or document |
| G14 | Re-running under a new run name rebuilds `storeDir` rescue databases and reruns about 300 task-hours of TBLASTN unless the cache is carried over | Compute reporting | handled by a launch recipe; `ni` support is open |

## C. Island, locus and clinker views

| ID | Issue | Status |
|---|---|---|
| V1 | Columns were labelled by the mmseqs representative (a median of 6 different strains per drawn locus; IDs not stable between runs) | fixed in PR #221 (display strain's own gene IDs; representative kept as the key). Not done: a user-chosen reference strain, and labels in clinker and the island tables |
| V2 | "Model difference" is 56% of locus cells in isl45_v1 (reviewer), because every gene between the block ends in the exemplar, including exemplar-private models, becomes a column | **decision**: down-weight private and rescue-only columns in row classes |
| V3 | 13 of 100 drawn loci use an outgroup genome as exemplar | fixed in PR #221 (ingroup preferred within a tier) |
| V4 | Exemplar rescue columns show no coordinate because `rescue_locations` is not passed | fixed in PR #221 |
| V5 | Sidebar IDs, titles, "families" and "variants" were undefined or used a different strain's coordinates | fixed in PR #221 |
| V6 | No Starship, biosynthetic-cluster or mobile-element markers in the views; no gene name or product; no per-strain table download | open; see the paper plan |

## D. Process lessons (for the supplement)

- A merge into a stacked branch is not a merge into `main` (PR #217 into #216's branch). Retargeted as PR #223.
- A bin-script edit does not change a Nextflow task hash, so a cached result can be reused from an older script version. Check the process command text or add a version argument.

## E. Pull requests (2026-10-08)

| PR | Content |
|---|---|
| #221 | Island and locus views: uniform labels, clearer sidebar, ingroup exemplars, rescue coordinates |
| #222 | Clinker slices draw genes for UniProt-ID proteomes |
| #223 | Losses signal and novelty_discovery singleton evidence (retarget of #217 to `main`) |
| #224 | Pairwise: `--evalue` reaches filter 1; report definitions; stale docs (depends on #223) |
| #225 | This folder |
| #226 | Pangenome report: openness rule, definitions, wording |
