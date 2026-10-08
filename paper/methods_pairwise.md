# Methods: pairwise novelty and loss pipeline (`--cluster_tool pairwise`)

Scope: the default `--cluster_tool pairwise` pathway of `main.nf`, in the novelty
direction and the loss direction. The pipeline also has `--cluster_tool mmseqs` and
`--cluster_tool novelty_discovery` pathways (`main.nf` lines 176-290). These are being
removed and are not described here.

Source of truth: the code at this checkout. Every statement comes from the cited file
unless it carries the tag [UNVERIFIED]. Tool defaults that the code does not set
(for example DIAMOND or BLAST built-in limits) are tagged [UNVERIFIED], because I did
not read the tool source or manual in this task.

---

## 1. Overview

The pipeline asks which genes are specific to a set of genomes. It also asks which
genes a set of genomes has lost.

The user puts each proteome into one of two groups. The ingroup (`GROUP` = `IN`) is
the lineage of interest. The outgroup (`GROUP` = `OUT`) is the comparison set.

The novelty direction looks for ingroup proteins that are present in most ingroup
proteomes and absent from every outgroup proteome (`workflows/search.nf`).

The loss direction looks for outgroup proteins that are present in most outgroup
proteomes and absent from the ingroup (`workflows/loss_search.nf`). It is the same
algorithm with the query and target groups swapped. The thresholds are not symmetric
(Section 4.4).

Both directions run in every pairwise analysis (`main.nf` lines 292-390).

The unit of the call is one protein in one source proteome. The pipeline calls
presence or absence of that protein's homolog in each other proteome. It then
clusters the candidates, checks the candidates against genome DNA with TBLASTN,
annotates them, and writes reports.

TBLASTN, the outgroup evidence record, and the context search do not remove
candidates. They add evidence to the reports (`workflows/summarize.nf`,
`--skip_tblastn_filter`; `lib/other_evidence.py` line 6; `bin/context_presence.py`
lines 5-8).

---

## 2. Inputs

### 2.1 Config CSV (`--config`, required)

`lib/config_parser.py` `parse_config()` reads these columns:

| Column | Required by parser | Use |
|---|---|---|
| `GROUP` | yes | Group role. Pairwise uses `IN` and `OUT`. `NEAR_INGROUP` and `BROAD_OUTGROUP` feed only the report-only context search. Aliases `TARGET`, `DISC_OUT`, `NEAR_IN`, `BROAD_OUT` are normalized (`main.nf` `normalizeGroup()`; `lib/config_parser.py` `GROUP_ALIASES`). |
| `Species` | yes | Label. |
| `Strain` | no | Label. |
| `Protein` | yes | Proteome FASTA basename. |
| `DNA` | no (per row) | Genome FASTA basename, for TBLASTN. |
| `Short` | yes | Proteome key used everywhere. The code comment says "unique ≤8-char ID" (`lib/config_parser.py` line 16). I found no code that enforces the length or uniqueness [UNVERIFIED]. |
| `TaxonGroup` | yes | Label. |
| `GFF3` | no | Report-only chromosome and start columns. |
| `SourceDB` | no | Report-only gene link-out. |
| `NCBI_TaxID` | no | Report-only taxonomy link; also passed to UniProt linking. |
| `UniProtDatGz` | no | Deprecated. `main.nf` lines 129-132 read it only as a `--restrict-proteome` value for UniProt linking. |

`build_presence_matrix.py` uses only `IN` and `OUT` rows (lines 291-301). The comment
there says this is deliberate.

### 2.2 Sequence files (`--data_dir`, required)

`main.nf` `resolve_fa()` looks for each basename in `data_dir`, then in the
subdirectories `pep/` and `proteins/` (proteins) or `dna/`, `genome/`, `scaffolds/`
(DNA). A missing file stops the run.

Proteomes: one protein FASTA per config row.

Genomes: optional per row. At launch, `main.nf` lines 82-88 require at least one `OUT`
row with DNA (novelty TBLASTN) and at least one `IN` row with DNA (loss TBLASTN). The
run stops otherwise.

GFF3: optional per row. It is not staged as a channel. The report scripts resolve it
from `data_dir` (`main.nf` lines 97-101). It only adds chromosome and start columns.

### 2.3 Optional databases

- `--pfam_hmm` (default `null`): Pfam-A HMM file for annotation.
- `--swissprot_dmnd` (default `null`): SwissProt DIAMOND database for annotation.
- `--modelorgs_config` (default `null`): YAML of model-organism gene names.
- `--uniprot_index` (default `null`): UniProt library index for linking.

When a database is not set, that annotation step writes an empty file
(`workflows/annotate.nf`).

---

## 3. Steps

Step numbers match `paper/flow_pairwise.md`.

### Step 1. Pairwise protein search (SEARCH / LOSS_SEARCH)

Purpose: find candidate homologs of every query-group protein in every other proteome.

Algorithm:
1. Take every query-group proteome. This is the ingroup in the novelty direction and
   the outgroup in the loss direction.
2. Pair it with every proteome in `IN ∪ OUT` except itself
   (`filter { meta_q.id != meta_t.id }`, `workflows/search.nf` and
   `workflows/loss_search.nf`).
3. Search all proteins of the query proteome against the target proteome.

The number of query-target pairs is |Q| × (N − 1), where Q is the query group and N
is |IN| + |OUT|.

Tool (`--run_tool`, default `phmmer`, `nextflow.config` line 13):

| Tool | Process | Command details read in code |
|---|---|---|
| phmmer | `PHMMER_SEARCH` (`modules/phmmer.nf`) | `phmmer --tblout --noali`. No `-E` is set, so phmmer's own reporting threshold applies [UNVERIFIED: tool default]. |
| diamond | `DIAMOND_MAKEDB`, `DIAMOND_SEARCH` (`modules/diamond.nf`) | `diamond blastp --outfmt 6 qseqid sseqid evalue bitscore length pident qcovhsp scovhsp qlen slen --evalue ${parse_evalue}` plus `--${diamond_sensitivity}` (default `very-sensitive`). No `--max-target-seqs` is set, so DIAMOND's own default applies [UNVERIFIED: tool default, believed to be 25]. One job per query proteome loops over all targets. |
| blast | `BLAST_MAKEDB`, `BLAST_SEARCH` (`modules/blast.nf`) | `blastp -outfmt "6 qseqid sseqid evalue bitscore length pident qcovhsp qlen slen" -evalue ${parse_evalue}`. No `-max_target_seqs` is set [UNVERIFIED: tool default]. |

Outputs: `<Q>_vs_<T>.phmmer.tblout.gz`, `<Q>_vs_<T>.diamond.tsv.gz` or
`<Q>_vs_<T>.blast.tsv.gz` in `search_cache/`.

### Step 2. Hit parsing (PARSE_HITS)

Purpose: convert raw hits to one common table and drop noise.

Algorithm (`bin/parse_hits.py`, `lib/hits.py`):
1. Parse each hit. phmmer: per-sequence E-value (column 5) and score (column 6) from
   `--tblout`. diamond and blast: tabular columns.
2. Keep a hit when `evalue <= parse_evalue` (`lib/hits.py` `filter_hits()`; the hit is
   dropped when `hit.evalue > cutoff`).
3. Write the alignment metrics `length pident qcov scov qlen slen` when the tool
   reports them. phmmer hits have blank metric cells. blastp has no `scov`.

Parameter: `--parse_evalue`, default `0.01`. The script's own default is `1e-15`,
but `modules/parse_hits.nf` always passes `params.parse_evalue`.

Output: `<Q>_vs_<T>.parsed.tsv` with columns `query_id target_id evalue bitscore
query_proteome target_proteome length pident qcov scov qlen slen`.

### Step 3. Self-search and paralog assignment (PHMMER_SELF / DIAMOND_SELF / BLAST_SELF, PARSE_SELF_HITS)

Purpose: find each query protein's best within-proteome paralog, for filter 2.

Algorithm:
1. Search each query-group proteome against itself (`modules/self_search.nf`).
   This runs for the ingroup in the novelty direction and for the outgroup in the
   loss direction.
2. For each query protein, ignore the self hit (`target_id == query_id`).
3. Keep the remaining hit with the lowest E-value as `paralog_protein_ID`
   (`bin/parse_self_hits.py`).
4. A protein with no non-self hit has no paralog and is omitted.

Tool parameters (hard-coded, not exposed as params):
- phmmer: `-E 100`. No hit-count limit.
- diamond: `--max-target-seqs 2 --evalue 100 --very-sensitive`. This is independent of
  `--diamond_sensitivity`.
- blast: `-max_target_seqs 2 -evalue 100`.

Output: `self_hits/<SHORT>.paralog_cutoffs.tsv` with columns `protein_ID
paralog_protein_ID bitscore evalue`. Only `paralog_protein_ID` is used. The E-value
column is informational (`bin/parse_self_hits.py` docstring;
`bin/build_presence_matrix.py` `load_paralog_info()`).

Each protein has at most one paralog. The pipeline does not use second-best or
further paralogs.

### Step 4. Presence matrix and candidate call (BUILD_PRESENCE_MATRIX)

Script: `bin/build_presence_matrix.py`. Process: `modules/build_presence_matrix.nf`.

Inputs: all `*.parsed.tsv`, all `*.paralog_cutoffs.tsv`, the config CSV, the query
group (`IN` or `OUT`), the query-group presence threshold, and the other-group
maximum.

Notation: q is a query protein in query proteome Q. t is a target protein in target
proteome T. E(q,t) is the E-value of the hit q→t. p = paralog(q).

Only hits whose `query_proteome` is in the query group are scored (line 323).

#### Filter 1: flat E-value

A hit survives when `E(q,t) < 1e-5` (strict less-than, line 327).

The value comes from the script constant `DEFAULT_EVALUE = 1e-5` (line 152).
`modules/build_presence_matrix.nf` does not pass `--default-evalue`. Therefore
`--evalue` does not change filter 1 (see Section 7).

#### Filter 2: paralog competition

Purpose: drop a hit that the query's own paralog explains better.

Lookup table (lines 313-320): `best_ev` is the minimum E-value over all parsed hits
(E ≤ `parse_evalue`), keyed by:
- scope `target`: (query proteome, query protein, target protein);
- scope `proteome`: (query proteome, query protein, target proteome).

Rule (lines 333-363). For each hit q→t that passed filter 1:
1. If q has no paralog p, keep the hit.
2. E_p = `best_ev[(Q, p, t)]` for scope `target`, or `best_ev[(Q, p, T)]` for scope
   `proteome`. If no such entry exists, keep the hit.
3. The hit is disqualified when `E_p < E(q,t)`.
4. A disqualified hit is rescued, and kept, when either arm holds:
   - floor arm: `E(q,t) <= paralog_rescue_evalue` (default `1e-20`). A value of 0 or
     null disables it (`modules/build_presence_matrix.nf` passes `?: 0`).
   - delta arm (off by default): `log10(E(q,t)) − log10(E_p) < paralog_rescue_delta`.
     When `E_p = 0` the delta is +inf and the arm never rescues.
5. Hits that remain disqualified get status `paralog_filtered` and are removed.

Parameters: `--paralog_competition_scope` (default `target`), `--paralog_rescue_evalue`
(default `1e-20`), `--paralog_rescue_delta` (default `null`).

#### Filter 3: other-group coverage floor (opt-in)

Purpose: drop narrow hits to the absence side, which look like shared-domain hits.

Rule (lines 369-392). When `--other_coverage_floor_qcov` is set to F > 0:
1. Judge only hits whose target proteome is in the other group (outgroup in the
   novelty direction, ingroup in the loss direction).
2. Remove a judged hit when `qcov < F` (query coverage in percent). Status becomes
   `coverage_floor`.
3. If any judged hit has no `qcov`, stop with an error.

Filter 3 runs after filter 2. Query-group hits are never judged.

Parameter: `--other_coverage_floor_qcov`, default `null` (off). The config comment
says "15 is the chosen value when enabled". `main.nf` line 73 stops the run when the
floor is set with `--run_tool phmmer`.

Side outputs:
- `<matrix>.coverage_floor_rejections.tsv`: rejected hits (only when the floor is set).
- `<matrix>.query_lowcov.tsv`: report-only. Per row, the number of query-group cells
  whose best `qcov` is below F (lines 403-428). Header-only when the floor is off.

#### Presence rule

A cell (q, T) is present (1) when at least one hit q→t, with t in T, survives
filters 1-3 (lines 430-433).

The self cell (q, Q) is always set to 1 (`all_present = hit_proteomes | {qp}`,
line 465). The search never compares a proteome with itself, so the self cell is
not based on a hit.

A protein gets a matrix row only when it has at least one surviving hit
(rows are built from the `presence` dictionary, line 464). A protein with no
surviving hit to any other proteome has no row.

The matrix has one column per `IN` and per `OUT` proteome, in sorted order.

#### Candidate rule

From lines 499-510, for each row:

```
query_frac = (number of query-group columns equal to 1) / N_query
             (this count includes the self cell)
other_frac = (number of other-group columns equal to 1) / N_other
             (other_frac = 0 when N_other = 0)

candidate  ⇔  query_frac >= min_frac  AND  other_frac <= other_max_frac
```

| Direction | Query group | `min_frac` | `other_max_frac` |
|---|---|---|---|
| Novelty | `IN` | `--ingroup_min_frac` (0.75) | `0.0`, hard-coded in `workflows/search.nf` |
| Loss | `OUT` | `--outgroup_min_frac` (0.75) | `--loss_ingroup_max_frac` (0.0) |

Candidate IDs are written as `<source_proteome>::<protein_id>`.

#### Outputs

Novelty / loss file names:
- `presence_matrix.tsv` / `loss_presence_matrix.tsv`: 0/1 matrix, all rows.
- `candidates.txt` / `loss_candidates.txt`.
- `presence_matrix.evalues.tsv` / `loss_...`: best surviving E-value per cell; blank
  for the self cell and absent cells. Report-only.
- `presence_matrix.targets.tsv` / `loss_...`: target protein of that best hit.
  Report-only.
- `presence_matrix.other_evidence.tsv.gz` / `loss_...`: see below.
- `presence_matrix.query_lowcov.tsv` / `loss_...`.

#### Other-group evidence record (issue #208)

`lib/other_evidence.py` `build_other_evidence()` writes one row per
(query protein, other-group proteome) cell that has at least one hit after filter 1.
It includes hits that filters 2 and 3 removed.

Columns: `protein_id source_proteome other_proteome evidence best_target_id
best_evalue best_bitscore best_qcov max_qcov n_hits status paralog_id
paralog_evalue paralog_delta`.

- `status` is `kept` when any hit in the cell was kept. Otherwise it is the status of
  the cell's lowest-E hit (`paralog_filtered` or `coverage_floor`).
- `paralog_delta = log10(best_evalue) − log10(paralog_evalue)`; +inf when the paralog
  E-value is 0.

This file does not change any call. Hits with E ≥ 1e-5 are not in it.

In `main.nf`, only the novelty-direction file is consumed (by SUMMARIZE and REPORT).
The loss-direction file is produced and published but not consumed.

### Step 5. Candidate extraction and clustering (CLUSTER / LOSS_CLUSTER)

Purpose: group near-identical candidates so TBLASTN runs once per group.
ADR: `docs/adr/0001-cluster-candidates-before-tblastn.md`.

Algorithm:
1. `bin/extract_candidates.py` reads `candidates.txt` and writes the sequences of the
   listed proteins from the query-group proteomes to `candidates.fa`
   (`loss_candidates.fa`). The FASTA ID is the original protein ID, without the
   source prefix.
2. `modules/mmseqs_cluster.nf` runs
   `mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0` on all candidates of
   all query-group proteomes together. These values are hard-coded. They are not the
   `family_*` params, which belong to the `mmseqs` pathway.
3. `bin/restore_mmseqs_cluster_ids.py` rewrites the IDs in `*_cluster.tsv` to the full
   FASTA header token. It stops on an ID collision (from its docstring; I did not
   read its body [UNVERIFIED in detail]).
4. An empty candidate FASTA gives empty outputs and no error.

Outputs (`clusters/`): `clusters_rep_seq.fasta`, `clusters_all_seqs.fasta`,
`clusters_cluster.tsv` (loss: `loss_clusters_*`). The cluster TSV maps
representative → member.

Clustering does not make a presence call. It defines (a) which sequences TBLASTN
searches and (b) the `family_id` shown in reports.

### Step 6. Genome validation with TBLASTN (VALIDATE / LOSS_VALIDATE)

Purpose: look for the candidate in the other group's genome DNA, including
unannotated genes.

Algorithm (`workflows/validate.nf`, `modules/tblastn.nf`):
1. Build one nucleotide BLAST database per genome (`makeblastdb -dbtype nucl`).
   Novelty: outgroup genomes. Loss: ingroup genomes.
2. Run one TBLASTN job per genome with all cluster representatives as the query:
   `tblastn -outfmt "6 qseqid sseqid evalue bitscore pident length qstart qend sstart
   send sframe qseq sseq" -evalue ${evalue}`. No `-max_target_seqs` or filtering flag
   is set [UNVERIFIED: tool defaults apply].
3. `bin/summarize_tblastn.py`: a representative "hits" a genome when it has at least
   one HSP with `evalue <= evalue` (less-or-equal; line 54).
4. Expansion to members: every member of that representative's cluster gets 1 for
   that genome (lines 169-182). A member never gets its own TBLASTN search.
5. Coverage (issue #208): for each (representative, genome) with hits, take the union
   of all HSP query spans that pass the E-value cutoff. Divide by the representative
   length. Cap at 1.0. Write the same value for every member (lines 184-199).

Parameter: `--evalue`, default `1e-5`. It sets both the TBLASTN `-evalue` and the
summary cutoff.

Outputs: `tblastn/<SHORT>.tblastn.tsv`; `tblastn_summary.tsv` (protein × genome 0/1);
`tblastn_summary.coverage.tsv.gz` (columns `protein_id rep_id genome n_hsps
best_evalue query_span_cov span_start span_end`); loss: `loss_tblastn_summary.tsv`
and `loss_tblastn_summary.coverage.tsv.gz`. `BUILD_ALIGNMENT_SHARDS` writes
`alignments/` or `loss_alignments/` for the report popup (script body not read
[UNVERIFIED in detail]).

TBLASTN does not remove candidates in either direction (Steps 8 and 9).

The span union does not check that HSPs come from one locus, contig or strand
(`parse_tblastn_spans()` collects all HSPs per representative per genome).

### Step 7. Functional annotation (ANNOTATE / LOSS_ANNOTATE)

Purpose: add names and domains.

Algorithm (`workflows/annotate.nf`, `bin/annotate_presence_matrix.py`):
1. Pfam (when `--pfam_hmm` is set): `hmmsearch --tblout --noali <Pfam-A.hmm>
   candidates.fa`. No `-E`, `--cut_ga` or other threshold is set, so hmmsearch's own
   reporting threshold applies [UNVERIFIED: tool default]. Every reported domain name
   is listed, once per name, in `Pfam_Names`, with `Pfam_Accessions` and
   `Pfam_Evalues`. MPI mode is on by default (`--hmm_mpi true`).
2. SwissProt (when `--swissprot_dmnd` is set): `diamond blastp --evalue 1e-5
   --max-target-seqs 1 --outfmt 6 qseqid sseqid stitle evalue bitscore`. No
   sensitivity flag is set. The best hit title becomes `Best_Swissprot`.
3. Model organisms (when `--modelorgs_config` is set): gene name lookup
   (`lib/model_organisms.py`, not read in detail [UNVERIFIED in detail]).
4. UniProt (when `--uniprot_index` is set): `UNIPROT_LINK` runs once per config
   proteome. It matches by accession, RefSeq ID, or exact sequence
   (`bin/uniprot_link.py` docstring). Columns `uniprot_*` are joined by
   (source proteome, protein ID).
5. `product_description` priority: model-organism name, then Pfam names, then the
   SwissProt title (lines 189-218).

Scope: the search steps (Pfam, SwissProt) use `candidates.fa` only. The output
`presence_matrix.function.tsv` (`loss_presence_matrix.function.tsv`) keeps every row
of the presence matrix. Non-candidate rows therefore have empty Pfam and SwissProt
columns. Model-organism and UniProt columns can be filled for any row.

Outputs: `presence_matrix.function.tsv`, `candidates.pfam.tblout`,
`candidates.swissprot.tsv` (loss: `loss_` prefix).

### Step 8. Per-species novelty tables (SUMMARIZE / MAKE_NOVELTIES), novelty only

Script: `bin/make_novelties.py`. It runs with `--skip_tblastn_filter`
(`workflows/summarize.nf`).

Algorithm, for each ingroup proteome S:
1. Keep rows with `source_proteome == S`.
2. Keep rows where `(ingroup columns equal to 1) / N_in >= ingroup_min_frac`.
3. Keep rows where every outgroup column is 0.
4. Do not drop rows with TBLASTN hits. List the genomes in `tblastn_outgroup_hits`.
5. Add `family_id`, `family_size`, `family_members` from the cluster TSV. A protein in
   a one-member cluster gets `family_size = 1` and an empty `family_id`.
6. Add the outgroup signal columns (below) when evidence files exist.

With the default `other_max_frac = 0.0`, steps 1-3 select the same proteins as
`candidates.txt` (same predicate on the same matrix).

Output: `novelties.<SHORT>.tsv`, one per ingroup proteome.

#### Outgroup signal class (issue #208)

`lib/other_evidence.py` `candidate_signal()` computes, per candidate:

- `other_protein_max_qcov`: max `max_qcov` over the candidate's evidence rows.
- `other_protein_n_filtered`: number of evidence rows with status other than `kept`.
- `other_protein_min_paralog_evalue`: min paralog E-value over those rows.
- `other_tblastn_max_cov`: max TBLASTN `query_span_cov` over genomes.
- `other_tblastn_n_genomes`: number of genomes with a TBLASTN hit.
- `other_signal_class`, from

```
C = max( { protein max_qcov over evidence rows }
         ∪ { 100 × TBLASTN query_span_cov over genomes } )

class = none         if there is no evidence row and no TBLASTN hit
        ''           if a hit exists but no coverage value exists (e.g. phmmer and no TBLASTN hit)
        broad        if C >= T
        domain_only  if C <  T
```

T = `--other_signal_qcov`, default `50`. The config comment says "observed median,
not validated".

The protein evidence includes only hits with E < 1e-5 (after filter 1). For a
novelty candidate every such cell was removed by filter 2 or filter 3, because any
kept cell would make the protein present in the outgroup.

`novelties.html` recomputes the same class in the browser from C and lets the viewer
change T (`lib/report_template.py` `osigClass()`). It is computed only for rows
flagged as novel (`lib/report_data.py` line 670).

The class never changes the novelty call.

### Step 9. Reports (REPORT)

Processes: `workflows/report.nf`.

- `novelties.html` (`bin/make_report.py`, `lib/report_data.py` `build_payload()`):
  one row per row of `presence_matrix.function.tsv`. A row is flagged novel when its
  `protein_id` appears in any `novelties.<SHORT>.tsv`. It shows presence, E-values,
  targets, TBLASTN hits, annotation, family, outgroup signal, query low-coverage
  counts, and context columns. Written twice: `results/` (offline) and `docs/`
  (`--online`, with the alignment popup).
- `core.html` (`bin/make_core_report.py`, `build_core_payload()`): rows of the
  novelty-direction matrix whose source is an ingroup proteome and whose presence
  fraction over all `IN` + `OUT` columns (self cell included) is
  `>= core_min_frac` (default `0.95`).
- `losses.html` (`bin/make_losses_report.py`, `build_losses_payload()`): recomputes the
  loss predicate on `loss_presence_matrix.function.tsv`: source in the outgroup,
  `out_frac >= outgroup_min_frac` and `in_frac <= loss_ingroup_max_frac`. A TBLASTN
  hit in an ingroup genome sets `tb_hit` and does not remove the row. It adds
  family-level `out_breadth` and `in_retained`.
- `summary.pdf` (`bin/make_pdf_report.py`, not read [UNVERIFIED in detail]); skipped
  when `--pdf_report false`.
- `report.html` landing page and `alignment.html` (`bin/make_index_report.py`, not
  read [UNVERIFIED in detail]).

### Step 10. Report-only context search (CONTEXT_SEARCH), novelty only

`workflows/context_search.nf`, `bin/context_presence.py`.

When the config has `NEAR_INGROUP` or `BROAD_OUTGROUP` rows, the novelty candidates
and their paralogs are searched against those proteomes. The same filters 1-3 are
applied. Here filter 1 uses `--evalue` (passed as `--default-evalue`, line 226). The
result (`context_presence.tsv`, `context_presence.evalues.tsv`) appears only in
`novelties.html`. It never changes `candidates.txt`. With no such rows the step is a
no-op.

---

## 4. Denominators and edge cases

### 4.1 What counts as N

- N_query: number of `IN` rows (novelty) or `OUT` rows (loss) in the config.
- N_other: number of `OUT` rows (novelty) or `IN` rows (loss).
- `NEAR_INGROUP` and `BROAD_OUTGROUP` rows are not in any denominator
  (`build_presence_matrix.py` lines 300-304).
- The source proteome is in the numerator and the denominator of `query_frac`. Its
  cell is always 1.
- `core.html` uses N = N_in + N_out as the denominator.

Required count of other query-group proteomes with a surviving hit, for
`min_frac = 0.75`: the protein needs `ceil(0.75 × N_query) − 1` other query-group
proteomes with a hit. Worked values from the formula: N_query = 2 needs 1 of 1;
N_query = 3 needs 2 of 2 (2/3 = 0.67 < 0.75); N_query = 4 needs 2 of 3;
N_query = 5 needs 3 of 4.

### 4.2 A protein with no surviving hit

It has no matrix row. It cannot be a candidate. It does not appear in any report.
The pipeline does not report true orphans (proteins with no detectable homolog in any
other proteome).

### 4.3 Small N

- N_query = 1 (novelty, `other_max_frac = 0.0`): the only possible targets are
  outgroup proteomes. A protein with a surviving hit has an outgroup cell equal to 1
  and fails. A protein with no hit has no row. So the novelty direction gives zero
  candidates by construction. The comment in `modules/mmseqs_cluster.nf` names this
  case ("a single-ingroup-species run where nothing cleared...").
- The same holds in the loss direction for N_out = 1 when `loss_ingroup_max_frac = 0`.
- N_other = 0: `other_frac` is 0, so every row that meets `min_frac` is a candidate.
  The launch check (Section 2.2) requires DNA in both groups, so an empty group stops
  the run before this point.
- With small N_query, each proteome changes `query_frac` in large steps (4.1). One
  missed hit in one proteome can remove a candidate.

### 4.4 Loss direction asymmetries

- Novelty requires strict absence: `other_max_frac` is hard-coded `0.0` in
  `workflows/search.nf`. No param changes it.
- Loss allows near-absence: `--loss_ingroup_max_frac` (default `0.0`) can be raised.
  `bin/make_losses_report.py` accepts values in [0, 1).
- Presence thresholds are separate params with equal defaults (`ingroup_min_frac` and
  `outgroup_min_frac`, both 0.75).
- The loss self-search and paralog filter use outgroup paralogs.
- Filter 3, when on, judges hits to ingroup proteomes in the loss direction.
- The loss direction has no `MAKE_NOVELTIES` step, no outgroup signal class, no
  query low-coverage column, and no context search in its report. Its other-evidence
  file is produced but not consumed.
- Loss TBLASTN searches ingroup genomes. The result is shown as a flag.

---

## 5. Parameter table

Defaults come from `nextflow.config` unless the row says "hard-coded".

| Param | Default | Step | Effect | Justification |
|---|---|---|---|---|
| `--run_tool` | `phmmer` | 1, 3 | Search tool for pairwise and self search. | No justification found. |
| `--parse_evalue` | `0.01` | 1, 2 | DIAMOND/BLAST `--evalue`; parse cutoff (`<=`). Also bounds which paralog hits enter `best_ev`. | No justification found ("loose noise ceiling", config comment). |
| `DEFAULT_EVALUE` (script) | `1e-5`, hard-coded | 4 (filter 1) | Hit kept when E < 1e-5. Not wired to `--evalue`. | No justification found. Convention. |
| `--evalue` | `1e-5` | 6, 10 | TBLASTN `-evalue` and summary cutoff (`<=`); context-search filter 1. | No justification found. Convention. |
| `--diamond_sensitivity` | `very-sensitive` | 1, 10 | DIAMOND sensitivity mode for pairwise and context search. | Data: config comment cites a 3-clade benchmark in NovInvenio_Investigations `notes/diamond-sensitivity/README.md` (not read here). |
| Self-search `-E/--evalue` | `100`, hard-coded | 3 | Loose ceiling so a paralog is found. | Design rationale in `modules/self_search.nf`; no data found for the value. |
| Self-search hit cap | `2` (diamond, blast), none (phmmer), hard-coded | 3 | Only the best non-self hit is used. | Design rationale only. |
| DIAMOND_SELF sensitivity | `--very-sensitive`, hard-coded | 3 | Paralog detection. | Data: one case (N. crassa HEX-1/eIF-5A) in `modules/self_search.nf` and `.living/decisions.md` 2026-09-10. |
| `--paralog_competition_scope` | `target` | 4 (filter 2) | Paralog must beat the query on the same target protein (`target`) or anywhere in the target proteome (`proteome`). | Data: one curated case (HEX-1/eIF5A) in `bin/build_presence_matrix.py` docstring. |
| `--paralog_rescue_evalue` | `1e-20` | 4 (filter 2) | Keeps a hit with E ≤ 1e-20 even if the paralog wins. | Data: pezizo_set1 comparison of floor vs delta (floor removes 213 candidates at 57.3% TBLASTN-breadth hit-rate). I found no sweep of the floor value itself. |
| `--paralog_rescue_delta` | `null` (off) | 4 (filter 2) | Keeps a hit the paralog beats by < delta orders of magnitude. | Data: measured non-selective on pezizo_set1 (47.1% vs 48.2%). |
| `--other_coverage_floor_qcov` | `null` (off); 15 recommended when on | 4 (filter 3) | Other-group hit with qcov < F counts as absent. | Data: pezizo_set1 BUSCO 1:1 pairs, 0.38% false rejection at 15, 1.18% at 20 (config comment). |
| `--ingroup_min_frac` | `0.75` | 4, 8, 9 | Novelty presence threshold in the ingroup. | No justification found. |
| novelty `other_max_frac` | `0.0`, hard-coded | 4 | Strict outgroup absence. | Definition of novelty; no data. |
| `--outgroup_min_frac` | `0.75` | 4 (loss), 9 | Loss presence threshold in the outgroup. | No justification found. |
| `--loss_ingroup_max_frac` | `0.0` | 4 (loss), 9 | Maximum ingroup fraction for a loss. | Convention (strict absence). |
| `--core_min_frac` | `0.95` | 9 | Core report threshold over all proteomes. | No justification found. |
| Candidate clustering | `--min-seq-id 0.3 -c 0.8 --cov-mode 0`, hard-coded | 5 | Defines TBLASTN representatives and report families. | No data found. ADR-0001 states the trade-off only. |
| `--other_signal_qcov` | `50` | 8, 9 | Threshold T for `broad` vs `domain_only`. | "Observed median, not validated" (config comment). |
| SwissProt search | `--evalue 1e-5 --max-target-seqs 1`, hard-coded | 7 | Best SwissProt hit. | No justification found. |
| Pfam search | no threshold flag, hard-coded | 7 | All domains hmmsearch reports. | No justification found. |
| `--pfam_hmm`, `--swissprot_dmnd`, `--modelorgs_config`, `--uniprot_index` | `null` | 7 | Enable each annotation source. | Not applicable. |
| `--report_sequences` | `novelties` | 9 | Which rows embed a sequence. | File-size rationale (config comment). |

---

## 6. Known limitations and caveats

1. **True orphans are not reported.** A protein needs at least one surviving hit to
   get a matrix row (4.2).
2. **One ingroup proteome gives no novelties** with the default strict absence (4.3).
3. **The search is not exhaustive.** Default `run_tool` is phmmer. DIAMOND and BLAST
   run without a hit-count setting, so the tools' own per-query target limits apply
   [UNVERIFIED: exact defaults]. A truncated hit list can hide a homolog or the
   paralog's hit that filter 2 needs.
4. **Absence is a protein-annotation call.** An outgroup gene that is missing from
   the outgroup protein FASTA reads as absent. TBLASTN reports such cases but does not
   remove the candidate. `.living/decisions.md` (2026-09-22) records that on
   pezizo_set1 36.6% of 1479 candidate × outgroup cells with TBLASTN-only evidence were
   annotation gaps (not re-measured here).
5. **TBLASTN covers only the representative.** Members inherit the representative's
   result. A member's own outgroup hit can be missed (ADR-0001). Coverage is a union
   of HSPs from any locus.
6. **Filter 2 uses one paralog per protein.** It uses only the best within-proteome
   hit. Other paralogs are not tested.
7. **Filter 2 compares E-values.** The query and its paralog are compared against the
   same target database, so database size is the same for both. phmmer reports
   per-sequence E-values only; DIAMOND can report E = 0 for strong hits (comment at
   `build_presence_matrix.py` line 354). A query with E = 0 is always rescued by the
   floor arm.
8. **Filter 1 is fixed at 1e-5.** `--evalue` does not change it (Section 7).
9. **Coverage-based evidence needs DIAMOND or BLAST.** With phmmer, filter 3 is
   rejected and the outgroup signal class relies on TBLASTN coverage only.
10. **Outgroup signal threshold is not validated.** T = 50 is an observed median.
11. **Pfam annotation has no significance threshold set in the pipeline.**
    `Pfam_Names` can list weak domain hits [UNVERIFIED: hmmsearch default reporting
    threshold].
12. **Core and non-candidate rows lack Pfam/SwissProt annotation.** The searches run
    on `candidates.fa` only.
13. **`core.html` is the ingroup view.** Rows exist only for ingroup-sourced proteins,
    and one gene can appear once per ingroup proteome (`build_core_payload()` docstring).
14. **Protein IDs are assumed unique across proteomes.** `paralog_of`
    (`load_paralog_info()`), `candidates.fa` headers, the cluster TSV, the TBLASTN
    summary, `make_novelties.py` TBLASTN lookup and the report's novelty flag are keyed
    by `protein_id` alone. I found no check for duplicate protein IDs across proteomes
    [UNVERIFIED: absence of such a check outside the files read]. The mmseqs ID
    restore step stops on a collapsed-ID collision (docstring).
15. **Clustering parameters are hard-coded** and not exposed as params.
16. **Threshold steps are coarse for small N** (4.1, 4.3).
17. **No statistical model.** Calls are threshold rules. No error rate is estimated for
    presence or absence.

---

## 7. Discrepancies and items to verify

Docs versus code:

- `METHOD_DESCRIPTION.md` §1 says filter 1 uses the flat `--evalue`. Code:
  `modules/build_presence_matrix.nf` does not pass `--default-evalue`; the script
  constant `1e-5` applies. `--evalue` changes TBLASTN and the context search only.
  The values match at the default, but a user who sets `--evalue` changes only part
  of the pipeline.
- `METHOD_DESCRIPTION.md` §1 says the default `--paralog-competition-scope` is
  `proteome`. Code default is `target` (`nextflow.config` line 34; argparse default
  in `build_presence_matrix.py` line 224). The module docstring of
  `bin/build_presence_matrix.py` (line 14) also still says "proteome (default)".
- `METHOD_DESCRIPTION.md` does not describe the paralog rescue floor (`1e-20`, on by
  default), the delta arm, filter 3, the query low-coverage report, the other-group
  evidence record, or the outgroup signal class.
- `METHOD_DESCRIPTION.md` says the self-search uses `-E 100 --max-target-seqs 2`.
  Code: phmmer self-search has `-E 100` and no hit cap.
- `METHOD_DESCRIPTION.md` says Pfam uses `hmmscan`. Code uses `hmmsearch`
  (`workflows/annotate.nf`).
- `METHOD_DESCRIPTION.md` says Pfam and SwissProt "annotate every row in the final
  matrix". Code searches `candidates.fa` only; other rows get empty Pfam/SwissProt
  columns.
- `METHOD_DESCRIPTION.md` says the loss direction is "the identical algorithm".
  Code differs in thresholds and downstream steps (Section 4.4).
- `METHOD_DESCRIPTION.md` §2 says the `mmseqs` pathway applies TBLASTN as a
  disqualifying filter. `workflows/summarize.nf` passes `--skip_tblastn_filter` for
  every `cluster_tool` (out of scope, noted for the removal work).
- `METHOD_DESCRIPTION.md` §2 lists `--family_chunk_size` default 200 and
  `hmm_presence_cov` 0.5. `nextflow.config` has 50 and 0.3 (out of scope).
- `bin/make_novelties.py` docstring rule 4 and the `workflows/summarize.nf` header
  comment say novelties are "absent from all outgroup genomes". The pipeline always
  passes `--skip_tblastn_filter`, so TBLASTN hits do not remove rows.
- `modules/self_search.nf` comment calls the pairwise DIAMOND search "diamond's
  default fast mode". The pairwise default is now `very-sensitive`
  (`nextflow.config` line 16).
- `lib/report_data.py` `LOSSES_ROW_FIELDS` comment for `seq` says "the loss matrix
  only ever contains candidate rows". `build_presence_matrix.py` writes every row with
  a surviving hit, and `build_losses_payload()`'s own docstring says the matrix is the
  full table.
- `docs/method_description.md` only includes `METHOD_DESCRIPTION.md`. It has no
  separate content.
- ADR-0001 says representatives are searched "against outgroup genomes". This is true
  for novelties. Losses search ingroup genomes.

Facts I could not determine from the code read:

- DIAMOND, BLASTP, TBLASTN, phmmer and hmmsearch built-in defaults (target limits,
  reporting E-value, low-complexity filtering) that apply because the pipeline does
  not set them.
- Whether `Short` length or protein-ID uniqueness is validated anywhere.
- The bodies of `bin/build_alignment_shards.py`, `bin/make_pdf_report.py`,
  `bin/make_index_report.py`, `lib/model_organisms.py`, `bin/restore_mmseqs_cluster_ids.py`
  (only its docstring), and `Helpers.projectName()` (output directory naming; the
  config comment says the project defaults to the config CSV basename).
- How the `1e-20` floor value was chosen. The record compares the floor to the delta
  arm; I found no sweep of floor values.
- The source of the `0.75`, `0.95`, `0.01` and `1e-5` defaults.
- The data behind the `very-sensitive` default (cited file in the other repository,
  not read).
