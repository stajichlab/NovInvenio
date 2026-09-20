# Paralog handling in the pairwise novelty call — problem analysis and option space

**Date**: 2026-09-20
**Status**: superseded in part. Sec 3's complementarity argument was measured and **rejected**
the same day -- see **Sec 6**, which is the operative conclusion. Implemented as issue #128.
The original reasoning is left intact below because it records why the delta arm was tried.
**Scope**: `--cluster_tool pairwise`, `bin/build_presence_matrix.py` filters 1 and 2.
**Data**: `NovInvenio_Investigations/results/pezizo_set1/` (diamond, 5 ingroup Pezizomycotina
vs 6 outgroups; 3393 novelty candidates from 35460 matrix rows).

Controls, with trees in
`NovInvenio_Investigations/studies/fungi/pezizo_set1/analysis/`:

- **Positive control — HEX1** (`sp|P87252|HEX1_NEUCR`, Woronin body protein, a
  Pezizomycotina-specific duplicate of eIF5A). Should be called novel.
- **Negative control — CorA** (`tr|A7UWR3|A7UWR3_NEUCR`, NCU11312, Mg2+/Mn2+ transporter,
  in-genome paralog `Q7SEN2`). Should NOT be called novel — the duplication predates the
  Pezizomycotina radiation and real orthologs exist in the outgroups.

---

## 1. How a novelty call is made today

The decision happens at **two levels**. Most of the confusion comes from mixing them up.

**Level 1 — cell level.** For every (ingroup query protein x outgroup proteome) pair, decide
whether the cell gets a `1` or a `0` in `presence_matrix.tsv`. Every raw hit is tested:

```
Filter 1 (absolute):  keep hit if  E < 1e-5
Filter 2 (relative):  drop hit if  the query's in-genome paralog
                                   hits the SAME target protein
                                   with a better E than the query does
```

A cell is `1` if any hit survives both filters. Otherwise `0`.

Filter 2's granularity is set by `--paralog-competition-scope` (`target` is the pipeline
default; `proteome` is the stricter original behaviour).

**Level 2 — protein level.** Once the row is built, a protein becomes a candidate if:

```
present in >= 75% of ingroup columns   AND   present in 0% of outgroup columns
```

(`--ingroup-min-frac`, `--other-max-frac`.)

The second clause is absolute. One surviving outgroup hit anywhere kills novelty.
**This is why a cell-level change is sufficient** — the candidate rule never has to change.

### The two controls through this machinery (Ncra query)

| | raw hits into outgroups | filter 1 | filter 2 | row | candidate? |
|---|---|---|---|---|---|
| **HEX1** P87252 | one, Mcir 3.4e-12 | passes | **drops it** — eIF5A (P38672) hits that same Mcir protein (S2K710) at 9.1e-69 | 0/6 | yes — correct |
| **A7UWR3** CorA | 15 hits across all 6 outgroups, best 3.1e-116 | all pass | **drops all 15** — Q7SEN2 wins every head-to-head | 0/6 | yes — WRONG |

Both are in `candidates.txt` today.

---

## 2. The defect

Filter 2 is purely relative. It asks *"did the paralog win?"* and never asks *"but was the
query's hit good anyway?"*

A7UWR3's best Scer hit is 3.3e-93 to `sp|P43553|ALR2_YEAST`. The tree
(`CorA_family_A7UWR3.aln.clipkit.treefile`) places ALR1/ALR2_YEAST as sister to the A7UWR3
ingroup clade — that hit is a **phylogenetically confirmed ortholog**. Filter 2 erases it
because a sibling gene scored better.

A relative test has no floor, so an arbitrarily strong hit can be erased by an arbitrarily
stronger one. The fix is to add a floor. The question is which floor.

### Filter 2 is load-bearing — do not remove it

HEX1's only outgroup hit (Mcir, 3.4e-12) passes filter 1. Filter 2 is what removes it and
gives HEX1 a clean 0/6 row. Delete filter 2 and HEX1 is present in 1/6 outgroups, fails
`--other-max-frac 0.0`, and stops being a candidate.

The two mechanisms act on **disjoint populations** of the 3393 candidates:

| best raw outgroup hit | n candidates | acted on by |
|---|---|---|
| E < 1e-20 | 213 (6.3%) | a new absolute floor |
| 1e-20 <= E < 1e-5 | ~121 (3.6%) | filter 2 only |
| no hit at all | 3180 (93.7%) | neither |

Filter 2's only remaining legitimate job is the twilight zone — which is exactly where HEX1
lives. **Add, do not retool.**

---

## 3. Option space

Options 1–3 are **cell-level rescues**: they re-admit a hit that filter 2 dropped.
Option 4 is a different axis.

### Option 1 — Absolute floor on the query's own hit

*"If the query's own hit is strong enough on its own, keep it, no matter what the paralog did."*

```
rescue if   query_E <= FLOOR          # e.g. 1e-20
```

**Already implemented** as `--paralog-rescue-evalue` in `bin/build_presence_matrix.py`,
wired through `modules/build_presence_matrix.nf`, exposed as
`params.paralog_rescue_evalue = null` in `nextflow.config`. Currently off.

| control | query E | rescued? | outcome |
|---|---|---|---|
| HEX1 | 3.4e-12 | no | stays novel — correct |
| A7UWR3 | 3.1e-116 | yes | not novel — correct |

75 orders of magnitude separate the controls. Any floor from 1e-15 to 1e-80 gives the same
answer on both.

Population effect at 1e-20: removes **213 of 3393** candidates (6.3%).

**Weakness**: E-value scales with protein length and database size, and differs across
`--run_tool`. Not scale-free.

### Option 2 — Delta / margin rule

*"If the paralog only narrowly beat the query, the query's hit was not explained away — keep it."*

```
delta = log10(query_E) - log10(paralog_E)   # orders of magnitude the paralog won by
rescue if   delta < D
```

Scale-free: immune to protein length and tool choice, which Option 1 is not.

Measured `delta` for all 334 candidates that have at least one filter-2 disqualification
(value shown is for that candidate's strongest suppressed hit):

| query E of suppressed hit | n | delta<10 | 10–30 | 30–50 | 50–100 | >=100 | median |
|---|---|---|---|---|---|---|---|
| < 1e-100 | 55 | 22 | 21 | 4 | 6 | 2 | 16 |
| 1e-100 … 1e-60 | 46 | 10 | 15 | 10 | 6 | 1 | 23 |
| 1e-60 … 1e-40 | 36 | 14 | 14 | 4 | 1 | 2 | 13 |
| 1e-40 … 1e-20 | 76 | 27 | 21 | 13 | 10 | 4 | 18 |
| 1e-20 … 1e-10 | 70 | 27 | 15 | 13 | 11 | 4 | 22 |
| 1e-10 … 1e-5 | 51 | 26 | 11 | 6 | 6 | 2 | 10 |

**Delta does not separate the populations.** The median is 10–23 in every bucket. A protein
with a 1e-100 hit has the same delta distribution as one with a 1e-8 hit.

On the controls it happens to work — HEX1's deltas are 56 and 59; A7UWR3's span 2 to 38 — so
`D ~ 45` separates them. But that threshold is unsupported: 8 of the 55 candidates with a hit
better than 1e-100 have delta >= 50 and would stay "novel" despite a 1e-100 hit.
**Delta alone under-rescues at the strong end.**

### Options 1 and 2 are complementary, not competing  -- **WRONG, see Sec 6**

They fail at opposite ends.

- **Floor alone** misses cells where the paralog barely won. A7UWR3 has a Ccin cell at
  `q=4.6e-13, p=3.9e-15, delta=2` — the paralog won by two orders, so the query's hit is
  plainly real evidence, but 4.6e-13 is above a 1e-20 floor and stays suppressed.
- **Delta alone** misses cells where the paralog crushed the query but the query's hit was
  still enormous.

OR-ing them covers both:

```
rescue if   query_E <= FLOOR    OR    delta < D
```

HEX1 survives both arms (E=3.4e-12 above any sane floor; delta=56 above any sane D), so its
novelty call does not depend on the threshold choice — which is what makes it a usable
positive control. Roughly 53 of the 121 twilight-zone candidates have at least one cell with
delta<10, so the delta arm adds real reach beyond the floor.

### Option 3 — Breadth rule (protein level, not cell level)

*"One strong outgroup hit is an anomaly. Six is an ortholog."*

```
disqualify if strong hits (E < 1e-20) exist in >= K outgroup proteomes
```

Breadth among the 213 affected candidates:

| outgroups hit | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| candidates | 98 | 34 | 21 | 21 | 11 | 28 |

A7UWR3 is 6/6 — maximally systematic. The 98 singletons are a different population: a lone
strong hit is what HGT, contamination, or one promiscuous domain match looks like.

More conservative than Option 1; drops fewer real novelties, but lets a genuine 2-species
ortholog through if K is set high. Composes with the others — a second knob on the *protein*
rule, not a replacement for the *cell* rule.

### Option 4 — Coverage-aware floor

*"A strong E-value over 15% of the protein is a shared domain, not orthology."*

```
rescue if   query_E <= FLOOR   AND   alignment covers >= X% of the query
```

The principled version. Protects against the failure mode Options 1–3 all share: **a
genuinely novel protein carrying one common domain** (kinase, WD40, ankyrin, P-loop NTPase)
will hit outgroup proteins at 1e-30 on that domain alone, and Options 1–3 will all wrongly
disqualify it. Same reasoning as `--hmm_presence_min_residues` in the mmseqs path.

**Cannot be built today.** `search_cache/*.diamond.tsv.gz` has four columns:
`query, target, evalue, bitscore`. No `length`, `pident`, `qcovhsp`. Adding them means
widening the diamond/blast/phmmer `outfmt` and invalidating the entire `storeDir` cache —
every pairwise search re-runs. Treat as a separate ticket that unblocks Option 4 later,
not a prerequisite for this work.

---

## 4. Recommended shape (superseded by Sec 6)

1. **Now**: turn on the OR of Options 1 and 2 at the cell level. Option 1 exists already;
   Option 2 is a small addition in the same code path
   (`build_presence_matrix.py`, filter-2 block).
2. **Keep filter 2.** Load-bearing for HEX1 — see Sec 2.
3. **Separate ticket**: widen pairwise search output to carry `length`, `pident`,
   `qcovhsp`; accept the cache rebuild; then revisit Option 4.
4. **State the caveat in any spec**: neither the floor nor the delta has an elbow in this
   data. Counts ramp smoothly — 55 candidates at 1e-100, 101 at 1e-60, 137 at 1e-40,
   175 at 1e-30, 213 at 1e-20. The only knee is at 1e-5/0.01 (334 -> 622), which is
   filter 1's territory. **1e-20 is a chosen convention, not a discovered threshold.**

---

## 5. Open questions, in priority order

1. **Validation.** Two curated controls set the direction but cannot set a number that will
   be applied to 3393 proteins. Can a larger labelled set be obtained — tree-checking a
   sample of the 213, or a clade where the answer is already known?
2. **Error budget.** Which mistake is worse: a CorA slipping through as novel, or a real
   novelty with a kinase domain being suppressed? These trade directly against each other
   and determine whether to tune strict or permissive.
3. **Presence vs. veto (secondary).** Rescuing makes the matrix cell `1` — honest, and it
   also changes `core.html` and the loss direction. A separate veto would keep the matrix
   unchanged and block only the candidate call. Park until the statistic is settled.

---

## Appendix — reproduction

Scripts used for the measurements above are throwaway; the inputs are:

- `NovInvenio_Investigations/results/pezizo_set1/search_cache/*_vs_*.diamond.tsv.gz`
  (4 cols: query, target, evalue, bitscore)
- `NovInvenio_Investigations/results/pezizo_set1/self_hits/*.paralog_cutoffs.tsv`
- `NovInvenio_Investigations/results/pezizo_set1/candidates.txt`

Method: reload raw hits, apply filter 1 (`E < 1e-5`), reconstruct filter 2 at
`scope=target` from the paralog table, and record every disqualification event for proteins
present in `candidates.txt`.

Control raw numbers (Ncra query, best hit per outgroup):

```
HEX1   sp|P87252|HEX1_NEUCR    Mcir 3.37e-12                          (1/6 outgroups)
A7UWR3 tr|A7UWR3|A7UWR3_NEUCR  Nirr 3.11e-116  Spom 2.93e-105
                               CneoH99 2.70e-96  Ccin 1.67e-95
                               Scer 3.35e-93    Mcir 4.88e-87         (6/6 outgroups)
```

Paralog pair strengths from the self-search (`Ncra.paralog_cutoffs.tsv`):

```
P87252 HEX1   <-> P38672 eIF5A    bitscore  43.9   E 2.84e-06   (ancient, diverged)
A7UWR3        <-> Q7SEN2          bitscore 297.0   E 1.12e-90   (well conserved)
```


---

## 6. Measured outcome (2026-09-20, same day) -- Option 2 rejected

Sec 3's claim that the floor and delta arms are complementary was reasoned from **one cell**
of A7UWR3 (Ccin, `q=4.6e-13, p=3.9e-15, delta=2`) and generalised. Measurement does not
support it.

### The error

The rescue fires **per (protein, proteome) cell**, and any one rescued cell ends a novelty
call. The value that triggers the rule is therefore the **minimum delta across all of a
candidate's suppressed cells**, not the delta at its strongest hit -- which is what Sec 3's
table reported. A minimum over several cells almost always lands on a small number: the
median min-delta across the 334 suppressed candidates is **8.6**, so `delta < 45` fires for
**291 of 334**. It is not a filter; it is close to undoing filter 2 entirely.

### Selectivity, scored against TBLASTN outgroup-genome breadth >= 4

| rule | removed | of those, breadth>=4 | hit-rate | breadth 6/6 |
|---|---|---|---|---|
| **floor 1e-20** | 213 | 122 | **57.3%** | 62 |
| delta<45 | 291 | 137 | 47.1% | 64 |
| delta<10 | 176 | 72 | 40.9% | 33 |
| floor OR delta<45 | 312 | 153 | 49.0% | 73 |
| *remove every suppressed candidate* | 334 | 161 | *48.2%* | 74 |

The last row is the indiscriminate baseline. **The floor is the only rule that beats it.**
Both delta variants fall below it, and OR-ing delta onto the floor drags 57.3% down to
49.0% -- the OR destroys the floor's selectivity rather than extending it.

Direction also runs against the rule: within strong-suppressed candidates, median min-delta
is 10.2 for breadth>=4 versus 6.1 for breadth<4 (AUC 0.64-0.67). Higher delta goes with
*more* genome-level presence, so "rescue if delta is small" preferentially rescues the
candidates least likely to be ancient paralogs.

### The controls could not have distinguished the two rules

| | min-delta | best hit | TBLASTN breadth | floor rescues? | delta<45 rescues? |
|---|---|---|---|---|---|
| A7UWR3 | 2.0 | 3.1e-116 | 6/6 | yes | yes |
| HEX1 | 56.4 | 3.4e-12 | 2/6 | no | no |

Both pass both. This is Sec 4's caveat cashing out: two curated controls set direction but
cannot validate a number applied to 3393 proteins. The caveat earned its place.

### Shipped (issue #128)

- Floor arm **ON by default**, `params.paralog_rescue_evalue = 1e-20`. `null` in config
  disables; the script takes `--paralog-rescue-evalue 0` as the off switch.
- Delta arm implemented but **OFF by default** (`--paralog-rescue-delta`), retained only so
  it can be re-swept once alignment coverage exists (#129). Arms OR when both are on.
- Option 4 unchanged and still blocked on #129.

### Caveat on this evaluation

The breadth>=4 proxy has an **unmeasured false-positive rate**. If TBLASTN is noisier than
assumed, the floor's 57.3%-vs-48.2% margin is softer than it looks. Measuring it is step 1
of issue #135, and is a prerequisite for trusting these numbers rather than a follow-up.

---

## 7. A separate, larger population found while validating (issue #135)

Of the 403 candidates with TBLASTN hits in >= 4 of 6 outgroup genomes, **242 have no
protein-level hit below 1e-5 at all** -- larger than, and disjoint from, the 213 the floor
addresses.

| | n | reading |
|---|---|---|
| no raw hit whatsoever, even out to E=0.01 | **208** (86%) | protein search sees nothing |
| raw hit 1e-5 .. 1e-3 | 27 | threshold near-miss on filter 1 |
| raw hit only 1e-3 .. 0.01 | 7 | noise floor |

Only 27 are near-misses, so this is not filter tuning. Four explanations need opposite
responses: (1) outgroup gene present but unannotated, (2) pseudogene/relic -- arguably a
*genuine* novelty, (3) TBLASTN false positive, (4) fragmented gene model. All four are
currently treated identically, as novelties, because `workflows/summarize.nf` runs
`MAKE_NOVELTIES` with `--skip_tblastn_filter`. See issue #135.

### Also noted

`lib/singleton_presence.py`'s `score_singleton_hits()` applies the same unconditional
filter 2 for the `--cluster_tool novelty_discovery` singleton branch with **no rescue arm**,
so the same protein can be called differently by the two pathways. Not addressed by #128.
