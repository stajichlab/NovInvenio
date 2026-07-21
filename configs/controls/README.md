# Validation control sets

Hand-curated positive/negative controls for validating the family-profile pathway
(`docs/adr/0002-family-profile-search-pathway.md`, Q8). One file per clade:
`<clade>.controls.csv`. Consumed by `bin/score_controls.py` (Phase 2 — not yet built;
see `todo/cross-method-support-column.md`).

## Purpose

- **positive** controls (`expected_call: novel`) — genes that *should* be flagged as
  lineage-specific novelty candidates. Measures **recall / sensitivity**.
- **negative** controls (`expected_call: core`) — conserved genes that must **not** be
  flagged. Measures **specificity / false-novelty rate**. BUSCO single-copy orthologs
  double as free negatives.

The scorer resolves each row's anchor to a family in a run, reads that family's call from
the presence matrix, and compares it to `expected_call` — one recall + FP number per
parameter-sweep grid-point, so the shipped default is chosen on biology, not just counts.

## Columns

| Column | Meaning |
|---|---|
| `control_id` | stable label. Rows beginning `EXAMPLE_` are placeholders — **delete or replace them**. |
| `class` | `positive` or `negative` |
| `expected_call` | the verdict the family should get: `novel` or `core` |
| `anchor_type` | how the gene is located: `protein_id`, `fasta`, or `busco` |
| `anchor` | the anchor value — a protein ID, a FASTA path under `configs/controls/seqs/`, or a BUSCO id |
| `proteome_short` | the ingroup `Short` a `protein_id` lives in (blank for `fasta`/`busco`) |
| `gene_name` | optional; usually empty for true novelties (clade-restricted genes are often uncharacterized) |
| `expected_origin` | clade/phylostratum where it should arise (sanity for the Phase-4 phylostrata work) |
| `source` | citation/DOI, "prior pairwise run", or "lab unpublished" — provenance |
| `notes` | free text |

## Curating tips

- **Anchor positives by `protein_id` or `fasta`, not `gene_name`** — most true novelties are
  hypothetical proteins with no name.
- Use `anchor_type: fasta` (sequence in `configs/controls/seqs/`) to promote a
  **silver-standard** hit (recovered by the existing pairwise pathway) into a durable
  **gold** control that survives re-annotation and ID changes.
- Negatives are cheapest as `anchor_type: busco` — `fungi_odb12` single-copy orthologs are
  conserved by definition and must never be called novel.
- Placeholder anchors are wrapped in `<...>`; the scorer skips `EXAMPLE_*` rows so the
  template is inert until you curate real entries.
