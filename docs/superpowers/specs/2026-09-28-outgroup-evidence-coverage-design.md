# Outgroup evidence record and coverage-based signal class

**Date**: 2026-09-28
**Status**: Design. Nothing implemented.
**Issue**: #208

## 1. Problem

A novelty candidate can have real outgroup hits that the pipeline removes.
Today the report does not show why. The viewer sees "no outgroup hit" and
cannot tell an orphan gene from a lineage-specific member of an old family.

Case: `A0A090CJI0_PODAN` in `NovInvenio_Investigations/results/sordariales_shallow`
(pairwise, DIAMOND).

- Raw DIAMOND finds it in all 8 outgroups (E 1e-7 to 2e-15, 53-79 bits).
- Filter 2 (paralog competition, scope `target`) removes all 33 outgroup hits.
  The in-genome paralog `B2B454_PODAN` hits the same target protein better,
  often at 1e-250 or lower. The rescue floor (1e-20) does not apply.
- TBLASTN (cluster rep `Q2GS92_CHAGB`, 517 aa) hits all 8 outgroup genomes
  (E 2e-7 to 3e-12). Every hit covers only query residues about 62-174,
  that is 19-25% of the protein.
- hmmscan: the shared region is Patatin (PF01734), found in the rep at 8.4e-17.
  The run's `candidates.pfam.tblout` has no row for this protein. The cause
  is not checked.

So the outgroup signal is one shared domain. The "novel" call is defensible
(same class as HEX-1), but the page must show this evidence.

Scale in the same run (TBLASTN, E <= 1e-5, union of HSP query spans per
genome, measured on cluster representatives only):

| | reps |
|---|---|
| at least 1 outgroup TBLASTN hit | 230 |
| best outgroup query coverage < 30% | 60 |
| best outgroup query coverage < 50% | 115 (median 0.50) |
| best outgroup query coverage >= 70% | 72 |

The protein-hit coverage of filter-2-removed cells is not measured yet.

## 2. Goals

1. Record every other-group hit that was judged, including hits that filters
   removed, with the reason and the paralog's E-value. Write it as a
   machine-readable TSV artifact.
2. Compute a per-candidate outgroup signal class from query coverage.
3. Add the class and its inputs as columns in the novelty TSVs and as a
   sortable, filterable column on `novelties.html`.
4. Show the per-gene evidence (including the paralog E-value) in a popup.
5. Do not change any presence call or candidate list. This is evidence only.

## 3. Data

### 3.1 New sidecar: `presence_matrix.other_evidence.tsv.gz`

Long format, one row per (query protein, other-group proteome, evidence type).
Written by `bin/build_presence_matrix.py` (new option `--output-other-evidence`).
The same code path runs for the loss direction (`loss_presence_matrix.other_evidence.tsv.gz`).

| column | meaning |
|---|---|
| `protein_id` | query protein |
| `source_proteome` | query proteome short name |
| `other_proteome` | other-group proteome short name |
| `evidence` | `protein` (this file); TBLASTN rows go in 3.2 |
| `best_target_id` | best-E-value target protein for this cell, before filters |
| `best_evalue`, `best_bitscore` | that hit's values |
| `best_qcov` | that hit's `qcovhsp` (percent). DIAMOND's default max-hsps is 1, so this is one HSP |
| `max_qcov` | highest `qcovhsp` over all hits in the cell with E < `--default-evalue` |
| `n_hits` | hits in the cell with E < `--default-evalue` |
| `status` | `kept`, `paralog_filtered`, `coverage_floor`, `below_evalue` (the cell's best hit's outcome; `kept` if any hit is kept) |
| `paralog_id` | the query's in-genome paralog (from `paralog_cutoffs.tsv`), blank if none |
| `paralog_evalue` | the paralog's E-value against `best_target_id` (scope `target`) or its best E-value in `other_proteome` (scope `proteome`), blank if no paralog hit |
| `paralog_delta` | `log10(best_evalue) - log10(paralog_evalue)`, blank if either is missing |

A cell with no hit at all gets no row. Compressed with gzip, since readers
outside this environment may use it.

The same writer goes into `lib/singleton_presence.py` (the
`--cluster_tool novelty_discovery` singleton branch). Filter 2 has three
copies; `bin/context_presence.py` is report-only and is out of scope here.
Family HMM presence has no per-hit qcov in this form; those rows get
`evidence=family_hmm` with profile coverage in `max_qcov` if it is cheap to
add, else no rows (state which in the implementation PR).

### 3.2 TBLASTN coverage: `tblastn_summary.coverage.tsv.gz`

Written by `bin/summarize_tblastn.py` (new option `--output-coverage`, new
input `--query_fasta` for query lengths). The TBLASTN `-outfmt` does not
change, so cached TBLASTN tasks stay valid on `-resume`.

Long format, one row per (cluster member protein, genome) with a hit at
E <= `--evalue`: `protein_id`, `rep_id`, `genome`, `n_hsps`, `best_evalue`,
`query_span_cov` (union of HSP `qstart..qend`, divided by rep length),
`span_start`, `span_end` (the merged span extremes). Rep values are copied to
members, as the presence summary does now.

### 3.3 Per-candidate columns in `novelties.<short>.tsv` (and loss TSVs)

Added by `bin/make_novelties.py` from 3.1 and 3.2:

| column | meaning |
|---|---|
| `other_protein_max_qcov` | max `max_qcov` over other-group proteomes, blank if no hit |
| `other_protein_n_filtered` | other-group proteomes where hits exist but all were removed |
| `other_protein_min_paralog_evalue` | lowest `paralog_evalue` over the filtered cells |
| `other_tblastn_max_cov` | max `query_span_cov` over other-group genomes, blank if no hit |
| `other_tblastn_n_genomes` | genomes with a TBLASTN hit (same as today's `tblastn_outgroup_hits` count) |
| `other_signal_class` | see 3.4 |

### 3.4 Signal class

Let `C` = max(`other_protein_max_qcov`, `other_tblastn_max_cov` x 100), using
only the evidence that exists. Threshold `T` = new param
`--other_signal_qcov`, default 50 (percent).

- `none`: no other-group protein hit and no TBLASTN hit.
- `domain_only`: some hit exists and `C < T`.
- `broad`: `C >= T`.

The default 50 is the median in section 1. No data yet says it separates
true from false novelty. The page lets the viewer change `T` (section 4), so
the default only sets the TSV column and the initial view.

## 4. Report (`novelties.html`, and `losses.html` with the same code)

- New column "Outgroup signal" with the class as a coloured tag. It sorts and
  filters like other columns.
- A slider or number box for `T`. The page reclassifies rows on the client
  from the per-row `C` values, so the page data carries `C`, not only the class.
- A popup on the tag (the same popup style as the TBLASTN alignment popup)
  that lists, per other-group proteome: best target, E-value, bits, qcov,
  status, paralog ID, paralog E-value, delta; and per genome the TBLASTN span
  and coverage. The page embeds only the candidate rows of 3.1/3.2, so size
  grows with candidate count (class 3 in NII terms, as now).
- `lib/report_data.py` reads the two sidecars. Missing or empty sidecars mean
  "not measured": the column shows blank, not `none`.

## 5. Tests

- `build_presence_matrix.py`: a fixture where the paralog wins on the same
  target writes one `paralog_filtered` row with the paralog E-value and delta;
  a kept hit writes `kept`; presence matrix and candidates are byte-identical
  with and without `--output-other-evidence`.
- `summarize_tblastn.py`: two overlapping HSPs and one separate HSP give the
  correct union coverage.
- `make_novelties.py`: class boundaries at `C = T - 1`, `T`, and no evidence.
- Report: page data carries `C` per row; missing sidecars give blank.
- Real data check: rerun `sordariales_shallow` with `-resume` and confirm
  `A0A090CJI0_PODAN` gets `domain_only`, `other_protein_n_filtered = 8`, and
  the paralog `B2B454_PODAN` E-values from section 1.

## 6. Out of scope

- Changing any presence call or default filter.
- `bin/context_presence.py` evidence.
- The pangenome pipeline.
- Choosing `T` from data. That needs a labelled control set (issue #135 notes).
