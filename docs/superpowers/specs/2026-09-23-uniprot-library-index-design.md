# UniProt library index and per-run linking

Date: 2026-09-23. Status: design approved in conversation; spec awaiting review.

## Goal

Every protein in a study that can be traced to a UniProt record gets that record's links and annotation in all three reports (novelties, core, losses). The linking runs from the pipeline alone, with no post-run step in NovInvenio_Investigations (NII). Links cover:

- the UniProt entry and the AlphaFold DB structure;
- the Pfam and InterPro domains that UniProt already computed;
- the other cross-reference databases on the record;
- publications linked to the record.

UniProt data comes from a library folder downloaded in advance. A Nextflow download step can be added later. It is out of scope here.

## Background

What already exists:

- `lib/uniprot_dat.py`: a single-pass `.dat` parser (AC/GN/DE/OX/DR lines). It extracts GO, Pfam, InterPro, EC, AlphaFold and 6 allow-listed cross-reference databases.
- `modules/uniprot_xref.nf` / `bin/build_uniprot_refseq_xref.py` (issue #92): per species, maps an NCBI RefSeq protein ID to its UniProt record through the record's `DR RefSeq` line. It runs only when the config row sets `UniProtDatGz`.
- `bin/annotate_presence_matrix.py --uniprot_xref_files` adds the `uniprot_*` columns, and `lib/report_data.py` reads them into the report payload (`go`, `ipr`, `ec`, `af`, `xrefs`).

The gaps:

1. **UniProt-sourced proteomes are not covered in the pipeline.** Their protein IDs are UniProt accessions (`tr|Q7S6W2|Q7S6W2_NEUCR`). Today NII merges their annotation after the run (`sync_reports.sh` → `merge_uniprot_annotations.py`), and only when the study directory has an `annotations/` folder. `sordariales_shallow_cluster` had none, so Q7S6W2's card showed no AlphaFold link there.
2. **Proteomes with other IDs are not covered at all**, such as BFD/funannotate gene IDs.
3. **The UniProt and AlphaFold buttons under "External resources" key off `o.sprot`**, the best SwissProt diamond hit, not the protein's own UniProt accession (`lib/report_common.py` `externalLinksNode`). A TrEMBL protein with no SwissProt hit never gets them.
4. **Only 6 cross-reference databases and no publications are extracted.** A *Neurospora* record carries more than 40 DR databases. Counts over UP000001805 (9759 records): InterPro 26139 lines, GO 24888, PANTHER 13046, OrthoDB 9191, Gene3D 8067, STRING 6646, SUPFAM 6610, PROSITE 5267, CDD 3836, SMART 2891, PDB 426.
5. **The `.dat` file is set by hand per species** (`UniProtDatGz`).

## Decisions

- **Coverage: all proteomes.** Match by protein ID where possible, and by exact protein-sequence identity for everything else.
- **Match scope: own species first, then any fungal proteome.** Take a record from the species' own taxid when one exists. Otherwise accept an identical sequence from any proteome in the library, and label it as such. Gene-database links, which are strain- and locus-specific, are shown only for own-species matches.
- **Approach: a one-time library index plus a per-run lookup step.** Rejected: parsing each species' `.dat` per run, which cannot support the cross-species fallback without scanning the whole library; and querying the UniProt REST API per protein (network access from compute nodes, rate limits, not reproducible).
- **`UniProtDatGz` is deprecated.** It stays accepted for one release as an alias that restricts the own-species proteome, and is then removed. `UNIPROT_XREF` is replaced by `UNIPROT_LINK`.

## Library

Current download: `/bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/`.

- UniProt release 2026_03 (02-Sep-2026), fungal reference proteomes.
- `data/<UP>_<taxid>.dat.gz` and `.fasta.gz`: 1526 proteomes, 17,344,385 canonical proteins, 19 GB.
- `fungi_proteomes_2026_03.csv`: `proteome_id, tax_id, species_name, oscode, n_canonical, n_isoform, file_prefix, lineage`.

The index design assumes only this layout (a folder of `<file_prefix>.dat.gz` plus a proteome CSV with `proteome_id`, `tax_id`, `species_name`, `file_prefix`). A later bacterial or all-eukaryote library with the same layout works unchanged.

## Components

### 1. Index build: `-entry UNIPROT_INDEX`

New workflow `workflows/uniprot_index.nf`, run on its own:

```bash
nextflow run main.nf -entry UNIPROT_INDEX -profile slurm \
    --uniprot_library /bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03 \
    --uniprot_library_csv fungi_proteomes_2026_03.csv \
    --uniprot_index /bigdata/stajichlab/shared/db/Uniprot/Fungi_2026_03/novinvenio_index/v1
```

- **`UNIPROT_PARSE_CHUNK`** runs `bin/uniprot_parse_dat.py` over a batch of proteomes.
  - Measured parse speed of the current parser is 9.8 MB of `.dat.gz` in 1.5 s (UP000001805, 9759 records), about 6.4 MB/s. The whole library then takes about 50 core-minutes. The added fields (sequence hash, publications, wider DR set) may double that.
  - Batches are sized by compressed bytes toward 1–1.5 h per job (`--uniprot_index_chunk_gb`, default 8), which gives about 2–4 jobs.
  - Each proteome writes `records/<UP>.records.tsv.zst`, one row per accession. Columns are listed under "Record fields".
  - Work happens in `$SCRATCH`, and the outputs are copied to `--uniprot_index`.
- **`UNIPROT_SEQ_INDEX`** merges all the records tables into `seq_index.sqlite`:
  - `seq(md5 TEXT, accession TEXT, proteome_id TEXT, taxid INTEGER, reviewed INTEGER)`, indexed on `md5`;
  - `acc(accession TEXT PRIMARY KEY, proteome_id TEXT, taxid INTEGER)`;
  - `refseq(refseq_id TEXT, accession TEXT)`, indexed on `refseq_id`;
  - `proteome(proteome_id TEXT PRIMARY KEY, taxid INTEGER, species_name TEXT, file_prefix TEXT)`.
- **`manifest.json`** records the library path, the UniProt release (read from the library README), the proteome CSV sha256, per-proteome record counts and `.dat.gz` sha256, the index format version (`1`) and the build date.
- The outputs are storeDir-cached on `--uniprot_index`, so the build runs once per release.

### 2. Record fields

These fields are written by `bin/uniprot_parse_dat.py` and implemented in an extended `lib/uniprot_dat.py`. Existing fields and their packing rules stay as they are.

| Field | Source | Notes |
|---|---|---|
| accession | AC (first) | |
| entry_name | ID | e.g. `Q7S6W2_NEUCR` |
| reviewed | ID line status | 1 = Swiss-Prot, 0 = TrEMBL |
| protein_existence | PE | 1–5 |
| taxon_id | OX | |
| gene_name | GN Name/ORFNames | existing rule |
| description | DE RecName/SubName | existing rule |
| ec_numbers | DE EC= | existing rule |
| go_ids | DR GO | existing rule, `GO:x:EVIDENCE` |
| pfam_ids, pfam_names | DR Pfam | existing rule |
| interpro_ids | DR InterPro | existing rule |
| alphafold_id | DR AlphaFoldDB | existing rule |
| xrefs | DR, allow-list | `DB:id\|…`; existing packing and RefSeq de-dup rules; allow-list widened (below) |
| pubs | RX, RT, RP | `PMID;DOI;scope;title` entries joined by `\|`. scope = `proteome` when RP says "NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]", else `protein` |
| seq_md5 | SQ block | MD5 of the upper-case sequence with no whitespace |
| seq_len | SQ line | |

The cross-reference allow-list becomes:
- existing: VEuPathDB, GeneID, KEGG, RefSeq, EnsemblFungi, EnsemblBacteria;
- added: PDB, PANTHER, OrthoDB, STRING, eggNOG, Gene3D, SUPFAM, PROSITE, SMART, CDD, PRINTS, PIRSF, HAMAP, NCBIfam, FunFam, MEROPS, CAZy, ESTHER, TCDB, BRENDA, UniPathway. All of these occur in UP000001805 or UP000002311 (checked 2026-09-23); PhylomeDB occurs in neither and is not included.

Each added database has an explicit first-field or last-field extraction rule, checked against real records in tests. Databases not on the list are dropped, which keeps the tables small.

### 3. Per-run linking: `UNIPROT_LINK`

New `modules/uniprot_link.nf`, process `UNIPROT_LINK`, running `bin/uniprot_link.py`. It runs once per proteome in the config (ingroup and outgroup, all GROUP roles) when `--uniprot_index` is set. When the parameter is unset, the step is skipped with a log note.

Inputs: the proteome FASTA, the meta (`id`, `species`, optional `NCBI_TaxID`, optional deprecated `UniProtDatGz`) and the index directory.

Match order per protein. The first rule that matches wins:

1. **`id`**: the protein ID is a UniProt accession. The pattern is `sp|ACC|…`, `tr|ACC|…`, or a bare accession matching the UniProt accession pattern. It must be present in `acc`.
2. **`refseq`**: the protein ID (first token, version kept) is in `refseq`.
3. **`seq_own`**: an exact `seq_md5` match whose `taxid` equals the row's own taxid. The own taxid comes from `NCBI_TaxID` if set. Otherwise it is looked up by exact `species_name` match in `proteome`, and failing that it is unknown. The deprecated `UniProtDatGz` restricts this rule to that one proteome.
4. **`seq_other`**: an exact `seq_md5` match in any proteome.

Tie rule, when one sequence matches several accessions: prefer the own taxid, then `reviewed = 1`, then the lowest accession in lexical order. The number of candidates is kept in `uniprot_n_matches`.

The full record fields for the chosen accessions come from `records/<UP>.records.tsv.zst`, reading only the proteomes that matched.

Output: `<Short>.uniprot_link.tsv`, keyed by the pipeline's `protein_id`, with these columns:
- the existing columns: `uniprot_accession, uniprot_gene_name, uniprot_description, uniprot_go_ids, uniprot_pfam_ids, uniprot_pfam_names, uniprot_interpro_ids, uniprot_ec_numbers, uniprot_alphafold_id, uniprot_xrefs`;
- new columns: `uniprot_match` (id | refseq | seq_own | seq_other), `uniprot_match_species`, `uniprot_reviewed`, `uniprot_pubs`, `uniprot_n_matches`.

Each task prints a coverage summary to stderr: the proteins per match type and the number with no match.

`main.nf` collects these outputs and passes them to `ANNOTATE`/`LOSS_ANNOTATE` in place of `UNIPROT_XREF`'s output. `bin/annotate_presence_matrix.py` already merges by `protein_id`. It adds the new columns to `UNIPROT_XREF_COLS`.

### 4. Reports

**`lib/report_data.py`** adds four row fields to all three payload builders:
- `umatch` (match type, or ''),
- `usp` (match species, when `seq_other`),
- `urev` (reviewed flag),
- `pubs` (an index into an interned `pub_sets` list, like `go_sets`).

**`lib/report_common.py`, `externalLinksNode`:**
- The UniProt and AlphaFold links use the protein's own matched accession (`uniprot_accession`, `af`). They fall back to `o.sprot` only when there is no match. When both exist and differ, the SwissProt hit is shown as "Best SwissProt hit: <acc>".
- Cross-reference links are rendered through an allow-listed URL-template table, one entry per database in the list above. IDs are URL-encoded, and URLs never come from the data.
- For `seq_other` matches, gene-database links (VEuPathDB, GeneID, KEGG, RefSeq, Ensembl*) are suppressed. The card shows "Identical sequence in <species> (<accession>)". Domain, family, structure and PDB links are kept.
- The Pfam and InterPro chips use the UniProt DR domains when the pipeline's own Pfam scan did not run or found nothing. The pipeline's own `Pfam_Names` take priority when present.
- A "Publications" field lists up to 5 `protein`-scope references as PubMed and DOI links. `proteome`-scope references are collapsed into one "Genome paper" line.
- All text is inserted with `textContent`, via the existing `el()` helper.

The NII post-run merge stays valid for studies run without `--uniprot_index`. For studies run with it, `sync_reports.sh` skips the merge. A follow-up NII change detects a `uniprot_match` column in `presence_matrix.function.tsv`.

## Error handling

- **No match** gives empty `uniprot_*` columns, never an error.
- **`--uniprot_index` set but unusable** is a hard error with a message naming the missing file or the version found. This covers a missing `manifest.json`, `seq_index.sqlite` or `records/`, and an unsupported `format_version`.
- **A proteome listed in the CSV with no `.dat.gz`** is skipped at build time with a warning. It appears in `manifest.json` under `missing`.
- **A malformed record** (no AC, or no SQ) is skipped, and the count is logged per proteome.

## Performance

Measured: 1.5 s per 9.8 MB `.dat.gz` with the current parser (single core).

Expected, to be confirmed on the first build:
- index build: about 1–2 core-hours in 2–4 jobs;
- `seq_index.sqlite`: about 17M rows;
- per-run lookup: one indexed query per protein.

The size of `seq_index.sqlite` and the per-proteome lookup time are recorded in `manifest.json` and in the first run's trace.

## Testing

- `tests/test_uniprot_dat.py` (extended), on fixture records cut from real files: UP000001805 (*Neurospora*, TrEMBL and one Swiss-Prot) and UP000002311 (*S. cerevisiae*, Ensembl isoform tags). It covers:
  - reviewed and protein-existence levels;
  - publication scope;
  - the sequence MD5;
  - each added DR database's field rule;
  - a record with no AlphaFold line.
- `tests/test_uniprot_link.py` on a small SQLite built by the index code from the fixtures. It covers each match rule, the tie rule, own-taxid lookup by `NCBI_TaxID` and by species name, the `UniProtDatGz` alias, and no match.
- `tests/data/uniprot_index/`: a fixture index used by `-profile test`.
- `tests/test_report_data.py`: the new payload fields.
- `tests/test_report_js_behaviour.py` (jsdom):
  - UniProt and AlphaFold links from the own accession, with no SwissProt hit;
  - gene-database links suppressed for `seq_other`;
  - a hostile publication title and a hostile cross-reference ID rendered inert.

## Acceptance

1. The `Fungi_2026_03` index builds, and `manifest.json` lists 1526 proteomes (or names any missing ones).
2. On `sordariales_shallow` and `sordariales_shallow_cluster`, Q7S6W2's card shows UniProt, AlphaFold and FungiDB links, plus the PANTHER/OrthoDB/STRING links present on its record. This works without an NII `annotations/` folder.
3. Per-species coverage is reported for one UniProt-sourced study (sordariales_shallow) and one non-UniProt study (a BFD/funannotate config). The non-UniProt coverage is a measured number, not assumed.

## Out of scope

- Downloading or refreshing the UniProt library from Nextflow.
- Near-identical matching (for example ≥ 95% identity via diamond). Only exact sequence identity is used.
- Isoform FASTA (`*_additional.fasta.gz`). Only canonical sequences are indexed.
- The NII enrichment scripts. They may later read `records/*.records.tsv.zst` instead of their own extraction.
