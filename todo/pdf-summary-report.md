# PDF summary report — publication-quality figures from the run

| Field | Value |
|-------|-------|
| **Date** | 2026-07-22 |
| **Author** | Jason Stajich |
| **Priority** | medium |
| **Status** | open (plan only — not implemented) |
| **Category** | feature / reporting |
| **Related analyses** | novelty / loss / core |
| **Related data** | `results/<project>/presence_matrix.function.tsv`, `tblastn_summary.tsv`, `novelties.*.tsv`, loss-direction equivalents |

## Description

Add a self-contained **PDF summary report** of publication-quality figures, complementary
to the interactive `novelties.html`/`core.html`/`losses.html`. The HTML reports are for
exploration; the PDF is for a PI/collaborator hand-off, a slide, or a paper figure panel —
something printable and embeddable that the interactive canvas heatmap is not.

## Motivation

The pipeline already assembles everything a summary needs in `lib/report_data.py`
(`build_payload()` / `build_core_payload()` / `build_losses_payload()`) — proteomes,
per-protein presence, novelty/loss calls, gene families, annotations, TBLASTN, and now the
cross-method `support` field. A PDF step needs **no new computation**, only rendering. A
static, paginated artifact is the missing output format for sharing off-cluster.

## Proposed approach (research — no implementation yet)

### Tooling decision

Surveyed options for a Python + Nextflow + pixi repo that already ships self-contained HTML:

| Approach | Fit | Notes |
|---|---|---|
| **matplotlib/seaborn → `PdfPages`** | ★★★ figures | pure-Python, publication-grade, multi-page; manual layout |
| **WeasyPrint (HTML+CSS → PDF)** | ★★★ layout | **reuses the existing report CSS/design tokens** → one visual system; embeds matplotlib PNG/SVG; static only (no canvas/JS) |
| Typst + matplotlib | ★★ polish | gorgeous, fast, lighter than LaTeX; adds a `typst` binary |
| Quarto / RMarkdown | ★★ narrative | reproducible narrative+figures; heavy (pandoc/LaTeX-or-typst) |
| ReportLab | ★ programmatic | total control; verbose |
| LaTeX (tectonic) | ★ | gold-standard typesetting; heaviest dep |

**Recommendation:** a two-part pipeline that **reuses `lib/report_data.py`** —
1. `bin/make_figures.py` — matplotlib/seaborn render a fixed figure set from the payload to
   PNG/SVG (embed the `dataviz` skill's palette + light/dark-safe, colour-blind-safe colours,
   matching the HTML reports' evidence colours: series-1 blue = search presence, series-2
   green = TBLASTN).
2. `bin/make_pdf_report.py` — assemble via **WeasyPrint** from an HTML+CSS template that
   reuses the report colour tokens, so the PDF reads as the same visual system as the HTML.
   Fallback: `matplotlib PdfPages` if a figures-only PDF (no HTML dependency) is preferred.

Add `weasyprint` (+ `matplotlib`/`seaborn`) to `pixi.toml [workspace]`. Wire a new
`workflows/pdf_report.nf` (`MAKE_PDF_REPORT`) gated like the HTML reports, publishing
`view/<project>/summary.pdf`.

### Figures to generate (domain-driven)

- Novelty candidates **per ingroup species** (bar).
- **Presence/absence heatmap** of the top-N novelty candidates (static twin of the canvas heatmap).
- **Gene-family size** distribution (histogram) and **Pfam-domain frequency** among novelties (top-N bar).
- **Cross-method concordance** — pairwise vs mmseqs novelty sets (UpSet or Venn), using the
  `support` field added in the report payload.
- **Loss** view: `out_breadth` vs `in_retained` scatter (clean-loss corner highlighted).
- **Annotation-source breakdown** (model-org / Pfam / SwissProt / none).
- A **run-summary header block** (ingroup/outgroup counts, search tool, thresholds) mirroring
  the HTML landing page.

## Acceptance criteria

- [ ] `summary.pdf` generated from an existing results dir with **no re-run** of search/annotation.
- [ ] Figures reuse the reports' colour language and are colour-blind- and print-safe.
- [ ] Self-contained (fonts/images embedded); opens/prints identically off-cluster.
- [ ] Regenerable standalone (like `make_report.py`) from `presence_matrix.function.tsv` + payload inputs.

## Notes

Depends on nothing new analytically — it is a rendering layer over `report_data.py`. Keep
the figure code in `bin/` + `lib/` (a `lib/figures.py` for the shared matplotlib styling) so
it is unit-testable on a fixture payload the way `report_data.py` is. Consult the `dataviz`
skill before writing any chart code. Cross-method concordance figure depends on both
pathways having run (or a two-matrix input like `make_report.py --support_matrix`).
