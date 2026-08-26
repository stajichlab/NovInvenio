Pipeline overview
=================

NovInvenio is a Nextflow DSL2 pipeline built around one pivot data
contract — a protein x proteome presence/absence matrix
(``presence_matrix.tsv``) and a candidate list (``candidates.txt``) — that
every downstream stage consumes regardless of which search strategy
produced it.

Stages
------

1. **Search** — pairwise proteome searches, family-profile HMM searches, or
   a two-phase target/screen search, depending on ``--cluster_tool``. See
   :doc:`method_description` for the full comparison.
2. **Cluster** — group candidate proteins (or read them directly off gene
   families, depending on pathway) ahead of genomic validation.
3. **Validate** — TBLASTN cluster/family representatives against outgroup
   genomes.
4. **Annotate** — Pfam (``hmmscan``), SwissProt (``diamond blastp``), and
   model-organism gene-name lookup merged onto the matrix.
5. **Summarize** — one novelty table per ingroup species.
6. **Report** — self-contained, offline-capable HTML reports
   (``novelties.html``, ``core.html``, ``losses.html``) plus a landing page
   in ``view/<project>/``.

Every stage above also runs in the mirrored **loss** direction (outgroup as
query, looking for genes present in the outgroup but absent from the
ingroup), except the ``novelty_discovery`` pathway, where the loss
direction is a deferred future extension.

See the repository's ``CLAUDE.md`` for the full data-flow reference
(channel contracts, module/workflow layout, parameter table) and
:doc:`method_description` for the search/cutoff logic itself.
