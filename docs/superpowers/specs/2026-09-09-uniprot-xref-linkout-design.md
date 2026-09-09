# General UniProt cross-reference linkout capability (design)

Date: 2026-09-09
Status: approved (rev. 2, post-review)

Revision note: rev. 1 was reviewed against the actual codebase and raw
UniProt data (see review summary below). The review found the central
`VEuPathDB` DR-line claim was wrong (would have shipped a dead FungiDB
link — the motivating use case), plus several under-specified extraction
rules and one design choice (interning) that added complexity without
benefit. This revision incorporates all nine required fixes.

## Problem

`lib/report_common.py`'s `externalLinksNode` renders a fixed, hardcoded set
of external links per candidate row (UniProt, AlphaFold, one model-org gene
URL, source-genome-DB link, NCBI Protein search, taxonomy). Every additional
external database (FungiDB, NCBI Gene, RefSeq, KEGG, ...) currently requires
a new hardcoded JS block plus, for FungiDB specifically, an entire
diamond-search-shaped workaround on the NII side
(`NovInvenio_Investigations/studies/fungi/pezizo_set1/bin/build_ncra_fungidb_crosswalk.py`)
to resolve a UniProt accession down to a FungiDB `NCU#####` gene ID via a
2-hop `id_transform: diamond_fasta` lookup — even though the source UniProt
record already carries that ID as a direct cross-reference.

## What the raw data actually contains

Evaluation of the cached Neurospora crassa proteome
(`NovInvenio_Investigations/data/uniprot/UP000001805/UP000001805_367110.dat.gz`,
9,759 records — corrected from rev. 1's 9,760) confirms `DR` cross-reference
lines carry FungiDB, NCBI Gene, RefSeq, KEGG, and Ensembl Fungi accessions,
at coverage: `VEuPathDB` 9,702/9,759, `GeneID` 9,704, `RefSeq` 10,234 (462
records carry more than one line — see §1), `KEGG` 9,553, `EnsemblFungi`
1,054 (10.8% — genuinely sparse, kept anyway since it's the only linkout for
non-FungiDB-indexed fungi). `PDB` (425 records, 4.4%) is **dropped from v1**
(see §1) after checking `UP000002311` (*S. cerevisiae*): `P23638` alone
carries 354 distinct `DR PDB` lines, so "one link per xref" does not scale to
this DB the way it does for the others.

This is not fungus-specific in structure, but the exact DB names are:
bacterial proteomes already in this repo (`UP000000265`, *K. pneumoniae*)
have zero `VEuPathDB`/`EnsemblFungi` lines and instead carry
`EnsemblBacteria` (and only 19% `GeneID` coverage, vs. 99% in Ncra) — the
allow-list below (§1) covers both, so the capability is repo-wide across the
fungal and bacterial studies both live in this repo.

## Scope

Two repos, one linked capability:

- **NovInvenio_Investigations (NII)** — extracts cross-references from the
  cached `.dat.gz` and carries them through its annotation-merge pipeline.
- **NovInvenio (nf_NovInvenio)** — renders them generically in the static
  HTML report.

This is additive: existing `modelorgs.yaml` configs (e.g. Ncra's
`diamond_fasta` self-identity-map workaround) are left in place and
untouched. No forced migration in this v1.

## 1. Extraction — NII `bin/extract_dat_annotations.py`

Extend the existing single-pass `DR`-line parser (already handles
GO/Pfam/InterPro/AlphaFoldDB) with a small, allow-listed set of additional
cross-reference types. Per-DB field layout is **not uniform** — this is
the part rev. 1 got wrong or left unspecified:

| DR type | Raw DR line example | Field taken | Notes |
|---|---|---|---|
| `VEuPathDB` | `VEuPathDB; FungiDB:NCU10683; -.` | field 2, **as-is, not stripped** | The DR value itself is `<Project>:<GeneID>` (e.g. `FungiDB:NCU10683`, or `PlasmoDB:...`/`ToxoDB:...` for other VEuPathDB projects). Store the whole `Project:GeneID` string as the id — do **not** strip the project prefix at extraction time. The renderer (§3) splits it again to pick both the correct host and the correct label, so a non-fungal VEuPathDB proteome renders correctly without a code change here. |
| `GeneID` | `GeneID; 5847462; -.` | field 2 | Bare NCBI Gene ID. |
| `RefSeq` | `RefSeq; XP_001728253.1; XM_001728201.2.` | field 2 **only** | Field 2 is the protein accession (`XP_`/`NP_`/`WP_`), field 3 the transcript/nucleotide accession (`XM_`/`NM_`/`NC_`). Taking both would point half the "RefSeq" links at nucleotide records under a protein-URL template. Keep only the **first** `RefSeq` line per record (462/9,759 records in Ncra alone carry more than one) — same "first hit wins" rule the parser already applies to `AlphaFoldDB` (`extract_dat_annotations.py:147`, `if alphafold_id is None`). |
| `KEGG` | `KEGG; ncr:NCU10683; -.` | field 2, as-is | Already `organism:locus`-shaped; kept intact (see §1's packing rule below — this is exactly why a fixed-arity split does not work here). |
| `EnsemblFungi` | `EnsemblFungi; EAA34304; EAA34304; NCU06768.` (Ncra: protein;protein;gene) or `EnsemblFungi; YML051W_mRNA; YML051W; YML051W.` (Scer: transcript;protein;gene) | **last** field | Field order/count is not stable across organisms (3 fields in some records, more in others) but the last field is consistently the stable gene ID, which is also what Ensembl's `/id/` resolver handles best. Extract via "split on `;`, take the last non-empty token before the trailing `.`" rather than a fixed field index. |
| `EnsemblBacteria` | same shape as `EnsemblFungi` | last field, same rule | Needed for this repo's bacterial studies (Akkermansia, K. pneumoniae) — same extraction logic, different DR type name. |

`PDB` is **not extracted in v1** — dropped rather than capped, since
structure links are already covered by the existing hardcoded AlphaFold link
and a "first N + N more" UI adds complexity for a DB this spec doesn't need
yet (see Deferred).

Sparse/bookkeeping databases beyond the above (STRING, PaxDb, OMA, HOGENOM,
InParanoid, OrthoDB, etc.) remain unextracted — orthology/expression
cross-references, not something a user clicks from a gene-info page.

**Packing rule (new — rev. 1 left this unstated and wrong-by-implication):**
`xrefs` packs `DB:id` pairs joined by `|`, where **only the first colon**
separates `DB` from `id` — `id` itself may contain further colons (`KEGG`,
`VEuPathDB`) and must not be split further at extraction time. This is a
different rule from the existing `go_ids` column, which is fixed-arity
(`GO:0016020:IEA`, always exactly 2 colons, always split into exactly 3
parts — see `report_common.py:302-305`); `xrefs` is variable-arity by
design, so unpacking code must use "split on first colon only"
(`indexOf(':')`), never a fixed `.split(':')`. Verified no allow-listed DR
value contains `|`, so the outer delimiter is safe.

Example packed value:
`VEuPathDB:FungiDB:NCU10683|GeneID:5847462|RefSeq:XP_001728253.1|KEGG:ncr:NCU10683`.

**Parser implementation note:** rather than adding 5 more sequential
`.match()` calls to the existing per-line regex chain (the module's own
docstring flags parser performance as a reason it isn't using Biopython;
Ncra alone has ~250k `DR` lines), dispatch on the DR line's own DB token
first (split once on `;`), then look up which extractor to run — O(1) per
line instead of up to 10 sequential regex attempts.

This directly obsoletes the *need* for `Ncra_self_id_crosswalk.tsv`'s
diamond-identity workaround for any *new* UniProt-backed model-organism
config (`xrefs` now carries a correctly-formed FungiDB reference per
protein) — but does not remove or migrate the existing Ncra config in this
v1. Regenerating `xrefs` requires re-running `extract_dat_annotations.py`
against **all** existing cached `.dat.gz` files, not just Ncra's — every
`studies/*/annotations/*.tsv` currently on disk predates this column.

## 2. Data flow (corrected — rev. 1 had the arrow backwards)

Real order, verified against `NovInvenio_Investigations/bin/sync_reports.sh`:

1. `extract_dat_annotations.py` (NII) produces the per-accession annotation
   TSV directly from the cached `.dat.gz` — independent of any Nextflow run.
2. nf_NovInvenio's `bin/annotate_presence_matrix.py` runs *inside* the
   pipeline and needs **no change**: it already preserves any columns not
   explicitly named in its own output field list
   (`out_fields = list(reader.fieldnames) + extra_cols`,
   `annotate_presence_matrix.py:120`).
3. `sync_reports.sh` runs *after* the pipeline, post-hoc: it invokes NII's
   `merge_uniprot_annotations.py` against `presence_matrix.function.tsv` and
   all `studies/*/annotations/*.tsv` (`sync_reports.sh:90-106`), then
   re-invokes `make_report.py` (line 124) and `make_core_report.py` (line
   139) to regenerate the static HTML from the merged matrix.

Changes:

- `merge_uniprot_annotations.py` gains one more passthrough column,
  `uniprot_xrefs`, using **`a.get("xrefs", "")`** (not `a["xrefs"]`) for the
  new field specifically — the existing code already mixes both styles
  (`a["pfam_ids"]` for older columns vs. `a.get(...)` for later additions at
  `merge_uniprot_annotations.py:76,78,81-82`); every annotation TSV on disk
  needs re-extraction (§1) before this runs, but `.get` keeps the merge from
  hard-failing on any not yet regenerated.
- `report_data.py`'s three field lists (`ROW_FIELDS`, `CORE_ROW_FIELDS`, and
  the losses-report equivalent — genuinely independent lists with different
  membership, not copies of one another, since they're indexed by name via
  each report's own `DATA.fields.forEach` at load time, not by position) each
  gain one more field, `xrefs`, via `row.get('uniprot_xrefs', '') or ''`.
- **No interning** (rev. 1's `xref_sets` pool is dropped). Interning
  (`go_sets`/`ipr_sets`) pays off only when a value repeats across many rows;
  `xrefs` values are 1:1 with the protein (GeneID/RefSeq/VEuPathDB/KEGG are
  all per-gene identifiers), so a pool would add a table entry *and* an
  index per row with nothing to compress — this is exactly the case
  `report_data.py`'s own comment on `ec`/`af` already calls out ("low
  per-row cardinality doesn't justify it"). `xrefs` follows the `ec`/`af`
  precedent: stored inline as a plain string.
- On the render side (all three call sites, not implied generically as in
  rev. 1): `report_template.py`, `core_report_template.py`, and
  `losses_report_template.py` each build an `o` object passed to
  `externalLinksNode` — each needs `xrefs: row[F.xrefs] || ''` added
  alongside its existing `sprot`/`gene_url` fields.

## 3. Renderer — nf_NovInvenio `lib/report_common.py`

Replace the "one hardcoded link per DB" pattern with a small registry plus a
generic loop. `VEuPathDB` needs special handling because its `id` is itself
a `Project:GeneID` pair (§1) — everything else is a flat template:

```js
var VEUPATHDB_HOSTS = {
  FungiDB: "fungidb.org/fungidb",
  PlasmoDB: "plasmodb.org/plasmo",
  ToxoDB: "toxodb.org/toxo"
  // extend as new VEuPathDB projects show up in cached proteomes
};
var XREF_LINK_TEMPLATES = {
  VEuPathDB: {
    render: function(id) {
      var i = id.indexOf(":");
      var project = i === -1 ? "VEuPathDB" : id.slice(0, i);
      var geneId = i === -1 ? id : id.slice(i + 1);
      var host = VEUPATHDB_HOSTS[project] || ("veupathdb.org/" + project.toLowerCase());
      return { label: project + " " + geneId,
               url: "https://" + host + "/app/record/gene/" + encodeURIComponent(geneId) };
    }
  },
  GeneID:    { label: "NCBI Gene", url: function(id){ return "https://www.ncbi.nlm.nih.gov/gene/" + encodeURIComponent(id); } },
  RefSeq:    { label: "RefSeq", url: function(id){ return "https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(id); } },
  KEGG:      { label: "KEGG", url: function(id){ return "https://www.genome.jp/dbget-bin/www_bget?" + encodeURIComponent(id); } },
  EnsemblFungi:    { label: "Ensembl Fungi", url: function(id){ return "https://fungi.ensembl.org/id/" + encodeURIComponent(id); } },
  EnsemblBacteria: { label: "Ensembl Bacteria", url: function(id){ return "https://bacteria.ensembl.org/id/" + encodeURIComponent(id); } }
};
```

`externalLinksNode` parses `o.xrefs` (the `DB:id|DB:id` string, splitting
each entry on the **first** colon only — §1) once, loops over entries, and
for each `DB` present in the registry appends a link (either via the
template's plain `label`/`url(id)` pair, or via `render(id)` for the
`VEuPathDB` special case above). A `DB` not in the registry (a future `DR`
type nobody has added a template for yet) is silently skipped — never shown
as raw/unlinked text. Every interpolated id goes through
`encodeURIComponent` (rev. 1 omitted this; `KEGG`'s `ncr:NCU10683` and any
future PDB-chain-bearing value need it — established convention already used
elsewhere in this file, e.g. `report_common.py:411,419`).

This sits alongside the existing UniProt/AlphaFold/model-org-gene/taxonomy
links; it does not replace them. **FungiDB duplication caveat** (stronger
than rev. 1's "harmless" framing): `genomeDbLink` (`report_common.py:406-420`)
already renders a FungiDB link off `SourceDB=fungidb` using
`geneIdFromProteinId(o.id)` — for a UniProt-sourced proteome that resolves
to the *UniProt accession*, not a real FungiDB gene ID, so that link is
already broken for these rows today. Once the new `xrefs`-derived FungiDB
link is added, a reader would see **two** FungiDB links on the same row, one
dead. Fix: when `o.xrefs` yields a `VEuPathDB` entry whose host matches the
row's `genomeDbLink` host, suppress the `genomeDbLink` render for that row
rather than showing both.

Adding a new database later is one registry entry; no loop-logic change.

## Testing

- `extract_dat_annotations.py`: unit test against a small synthetic `.dat`
  fixture covering: (a) multiple DR types packed into one `xrefs` value, (b)
  a record with 2+ `RefSeq` lines keeping only the first, (c) `EnsemblFungi`
  field-order differing between two synthetic records both resolving to the
  last field, (d) a non-allow-listed DB type producing no `xrefs` entry.
- `report_common.py`: extend `tests/test_report_js_behaviour.py` with a case
  asserting the registry loop (a) renders known DBs including the
  `VEuPathDB` project-split special case, (b) skips unknown DBs, (c)
  URL-encodes a `KEGG`-style id containing a colon. Note:
  `test_report_js_behaviour.py` self-skips when `jsdom` isn't importable
  (confirmed not installed in this environment) — this test won't actually
  execute here without installing it first; flag rather than assume green.
- Manual/integration check: re-run `extract_dat_annotations.py` for all
  proteomes referenced by `pezizo_set1` (not just Ncra), re-run
  `sync_reports.sh`'s merge step, confirm `xrefs` contains a correctly
  split `VEuPathDB:FungiDB:NCU10683`-shaped entry end to end, and spot-check
  the rendered report shows one working FungiDB link (not two, not dead).

## Deferred / out of scope

- Extracting `PDB` cross-references — dropped for v1 rather than capped;
  revisit with an explicit "first N + N more" UI if a study specifically
  needs structure linkouts beyond the existing AlphaFold link.
- Migrating existing `modelorgs.yaml` `diamond_fasta` configs (e.g. Ncra) to
  rely on `xrefs` instead of the diamond-identity-map workaround.
- Extracting the long tail of sparse DR types (STRING, PaxDb, BioCyc, CAZy,
  etc.) — revisit only if a specific study need arises.
- Any UI beyond a flat link list (grouping/collapsing many xrefs per DB) —
  not needed once `PDB` (the one DB with unbounded per-protein multiplicity)
  is dropped from v1.
