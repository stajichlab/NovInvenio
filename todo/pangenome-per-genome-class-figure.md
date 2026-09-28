# Per-genome gene-class composition figure in the pangenome report

- **Priority**: medium
- **Status**: in-progress (2026-09-28)
- **Category**: pangenome report / visualization
- **Date**: 2026-09-28
- **Author**: Jason Stajich
- **See**: `bin/pangenome_report_render.py`, `per_strain_summary.tsv`,
  `todo/pangenome-phylogeny-aware-gain-loss.md`; model figure: PPanGGOLiN
  (Gautreau et al. 2020, PLoS Comput Biol 16:e1007732, Fig 3)

## What

One horizontal stacked bar per genome in the pangenome summary report:
core / soft_core / shell / cloud / singleton family counts, read directly
from `per_strain_summary.tsv`. No tree is drawn; the model figure's value
is the per-genome class bars, not the tree.

## Genome order (first version)

1. Samplesheet `GROUP`: ingroup, then outgroup, then any other value.
2. Samplesheet `Species` within each `GROUP`, so a multi-species ingroup
   (e.g. Coccidioides: C. immitis + C. posadasii, all `IN` in the genus
   run) is split into one block per species.
3. Total family count within each block.

Without a samplesheet, genomes are ordered by total family count only.

## Later: order by evolutionary relatedness

Jason is exploring a method to build a phylogeny directly from the genomes.
When a strain tree is available (from that method, or the existing
`--pangenome_species_tree` input used for Dollo polarization, see
`todo/pangenome-phylogeny-aware-gain-loss.md`), order the bars by tree tip
order instead of step 3 above (or instead of steps 2-3). Drawing the tree
itself next to the bars, as in the model figure, is optional at that point.
Open questions for then:
- Does the genome-based tree method run fast enough to be a default
  pipeline step, or stay an optional input?
- Tip names must match samplesheet `Short`; reuse the tip-mapping and
  pruning rules from #185.
