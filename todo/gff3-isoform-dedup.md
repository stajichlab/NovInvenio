# Dedup/filter alternative splice isoforms in GFF3-derived chrom/start rows

| Field | Value |
|-------|-------|
| **Date** | 2026-08-07 |
| **Author** | Jason Stajich |
| **Priority** | low |
| **Status** | open |
| **Category** | chore |
| **Related** | `lib/gff3_genes.py`, `lib/report_data.py`, CLAUDE.md "GFF3 chrom/start is per-protein-record, not per-gene" |

## Idea

The Chromosome/Start columns added to `novelties.html`/`core.html`/`losses.html` are
resolved per protein/transcript ID (`lib/gff3_genes.py`'s `lookup_gene_position()`), and
nothing in the pipeline currently deduplicates or filters alternative splice isoforms. A
single gene with multiple annotated transcripts appears as multiple report rows — one per
protein/transcript ID — each pointing to the same or a very similar chrom/start.

## Why this matters

Not incorrect, just redundant: a user sorting/filtering by genomic position sees repeated
near-identical rows for multi-isoform genes, which can skew apparent novelty-candidate
counts per locus and clutter the position-sorted view.

## Possible approach

Collapse rows sharing the same (chrom, start) — or the same GFF3 gene-feature ID once
resolved — within a source proteome before emitting `payload['rows']`, keeping one
representative per locus (e.g. longest sequence, or first by protein ID). Needs a decision
on whether this happens in `report_data.py`'s payload builders (report-only, no upstream
pipeline change) or earlier (candidates.txt / presence matrix).
