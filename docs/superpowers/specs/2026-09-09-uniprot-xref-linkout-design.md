# General UniProt cross-reference linkout capability (design)

Date: 2026-09-09
Status: approved

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
record already carries that exact ID as a direct cross-reference.

Evaluation of the cached Neurospora crassa proteome
(`NovInvenio_Investigations/data/uniprot/UP000001805/UP000001805_367110.dat.gz`,
9,760 records) confirms this: `DR` cross-reference lines already carry
`VEuPathDB:NCU10683` (FungiDB), `GeneID`, `RefSeq`, `KEGG`, `EnsemblFungi`,
and `PDB` accessions directly, at high coverage (VEuPathDB 9,702/9,760,
GeneID 9,704, RefSeq 10,234, KEGG 9,553). This is not Neurospora-specific —
every UniProt-sourced proteome pulled by NII's
`bin/fetch_uniprot_proteome.py` carries the same `DR` line structure, so
extracting these cross-references is a repo-wide capability, not a one-off
fix for Ncra.

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
cross-reference types, chosen from the observed frequency table:

| DR type | Rendered as |
|---|---|
| `VEuPathDB` | FungiDB |
| `GeneID` | NCBI Gene |
| `RefSeq` | RefSeq |
| `KEGG` | KEGG |
| `EnsemblFungi` | Ensembl Fungi |
| `PDB` | PDB |

Sparse/bookkeeping databases (STRING, PaxDb, OMA, HOGENOM, InParanoid,
OrthoDB, etc.) are deliberately not extracted — they are orthology/expression
cross-references, not something a user would click from a gene-info page.

New output column: `xrefs`, packed as `DB:id` pairs joined by `|` — the same
packing convention the file already uses for `pfam_ids`/`go_ids` (no new
delimiter concept introduced). Example:
`VEuPathDB:NCU10683|GeneID:5847462|RefSeq:XP_001728253.1|KEGG:ncr:NCU10683`.
A protein carrying two accessions from the same DB (e.g. two RefSeq entries)
keeps both pairs — no deduplication or truncation.

This directly obsoletes the *need* for `Ncra_self_id_crosswalk.tsv`'s
diamond-identity workaround for any *new* UniProt-backed model-organism
config (`xrefs` already carries `VEuPathDB:NCU10683` per protein) — but does
not remove or migrate the existing Ncra config in this v1.

## 2. Data flow — NII `bin/merge_uniprot_annotations.py` → nf_NovInvenio `bin/annotate_presence_matrix.py` → `lib/report_data.py`

- `merge_uniprot_annotations.py` gains one more passthrough column:
  `uniprot_xrefs` (same pattern as its existing `uniprot_pfam_ids` etc.).
- `report_data.py`'s three near-duplicate `ROW_FIELDS` blocks each gain one
  more positional field, `xrefs`, populated via
  `row.get('uniprot_xrefs', '') or ''` (the same pattern already used for
  `ec`/`af`), interned through a new `xref_sets` string pool (the same
  approach already used for `go_sets`/`ipr_sets`) since the same xref string
  repeats across many rows of a shared model organism/reference proteome.

## 3. Renderer — nf_NovInvenio `lib/report_common.py`

Replace the "one hardcoded link per DB" pattern with a small registry plus a
generic loop:

```js
var XREF_LINK_TEMPLATES = {
  VEuPathDB: { label: "FungiDB", url: function(id){ return "https://fungidb.org/fungidb/app/record/gene/" + id; } },
  GeneID:    { label: "NCBI Gene", url: function(id){ return "https://www.ncbi.nlm.nih.gov/gene/" + id; } },
  RefSeq:    { label: "RefSeq", url: function(id){ return "https://www.ncbi.nlm.nih.gov/protein/" + id; } },
  KEGG:      { label: "KEGG", url: function(id){ return "https://www.genome.jp/dbget-bin/www_bget?" + id; } },
  EnsemblFungi: { label: "Ensembl Fungi", url: function(id){ return "https://fungi.ensembl.org/id/" + id; } },
  PDB:       { label: "PDB", url: function(id){ return "https://www.rcsb.org/structure/" + id; } }
};
```

`externalLinksNode` parses `o.xrefs` (the `DB:id|DB:id` string) once, loops
over entries, and for each `DB` present in the registry appends
`extLink(template.label + " " + id, template.url(id))`. A `DB` not in the
registry (a future `DR` type nobody has added a template for yet) is
silently skipped — never shown as raw/unlinked text. This sits alongside the
existing UniProt/AlphaFold/model-org-gene/taxonomy links; it does not
replace them. FungiDB may appear via both the old `gene_url_template` path
(for existing `modelorgs.yaml` configs) and this new generic path (for
proteins whose own UniProt record carries a `VEuPathDB` DR line) — harmless
duplication across two independent code paths keyed off different payload
fields, not a conflict.

Adding a 7th database later is one registry line; no loop-logic change.

## Testing

- `extract_dat_annotations.py`: unit test against a small synthetic `.dat`
  fixture asserting `xrefs` packs multiple DR types correctly in one pass
  and skips non-allow-listed DB types.
- `report_common.py`: extend `tests/test_report_js_behaviour.py` with a case
  asserting the registry loop renders known DBs and skips unknown ones.
- Manual/integration check: rebuild `pezizo_set1`'s Ncra annotation TSV,
  confirm `xrefs` contains `VEuPathDB:NCU10683`-style entries end to end, and
  spot-check the rendered report shows a working FungiDB link.

## Deferred / out of scope

- Migrating existing `modelorgs.yaml` `diamond_fasta` configs (e.g. Ncra) to
  rely on `xrefs` instead of the diamond-identity-map workaround.
- Extracting the long tail of sparse DR types (STRING, PaxDb, BioCyc, CAZy,
  etc.) — revisit only if a specific study need arises.
- Any UI beyond a flat link list (e.g. grouping/collapsing many xrefs) — not
  needed at current per-protein xref counts (typically ≤6).
