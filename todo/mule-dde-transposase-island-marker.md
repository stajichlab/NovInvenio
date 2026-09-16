# Screen accessory islands for MULE/DDE transposase proximity, not just known captains

| Field | Value |
|-------|-------|
| **Date** | 2026-09-16 |
| **Author** | Jason Stajich |
| **Priority** | idea |
| **Status** | open |
| **Category** | analysis |
| **Related analyses** | studies/fungi/coccidioides_pangenome/ (generalized island+Pfam step), studies/fungi/Afumigatus_pangenome/ (captain/SM-backbone marker precedent) |
| **Related data** | — |

## Description

The pangenome island/Pfam pipeline step being generalized (see the
2026-09-16 design work on `BUILD_ISLANDS`/marker-gene handling) currently
only has concrete marker searches for the Afumigatus study's own case:
captain genes (DUF3435/Starship) and SM-backbone (PKS/NRPS,
PF00109+PF00668). For Coccidioides — and future studies generally — there
isn't yet a clear equivalent "known mobile-element marker" to search for.

Worth exploring later: Transposase/DDE-domain families, particularly MULE
(Mutator-like transposase) DNA transposases, as a marker to check for
proximity to accessory islands. DDE transposases are a broad,
well-characterized functional class (multiple Pfam accessions, not one
HMM the way DUF3435 is for Starships specifically), so this isn't a
drop-in reuse of the existing single-HMM marker mechanism — it may need
either a curated small set of DDE/MULE-specific Pfam accessions, or a
broader domain-family screen.

## Motivation

Mobile genetic elements (transposons, not just fungal Starships) are a
plausible mechanistic driver of some accessory-island gene clustering
patterns generally, not just in Aspergillus. If MULE/DDE transposase
domains turn up enriched near or within flagged islands, that's a
biologically meaningful signal about how those islands may have formed
or moved, independent of whether a Starship-specific captain is present.

## Proposed Approach

Not scoped yet — flagged as exploration, not a committed task. Rough
shape once picked up:
- Identify the right Pfam accession(s) for MULE/DDE transposase domains
  (there are several DDE-superfamily Pfam families; needs a literature/
  Pfam-clan check, not just one accession the way DUF3435 was).
- If the generalized N-named-marker-search design (from the 2026-09-16
  island/Pfam architecture review) lands first, this becomes one more
  named marker set (`--pangenome_island_marker_hmm_names`-style), reusing
  the existing generic `MARKER_HMMSEARCH`-style module rather than new
  pipeline code.
- Otherwise, this could start as a study-specific ad hoc script (matching
  how Afumigatus's own captain/SM-backbone searches started before any
  generalization) rather than blocking on the shared-pipeline work.

## Acceptance Criteria

- [ ] Decide which Pfam accession(s) actually represent MULE/DDE
      transposases (not yet identified)
- [ ] Run the search against at least one real dataset's accessory
      islands and report whether proximity/enrichment shows a pattern
- [ ] Decide whether this belongs as a permanent named marker set in the
      generalized pipeline step, or stays study-specific

## Notes

Raised during the 2026-09-16 design discussion for generalizing the
Afumigatus study's manual island/Pfam-enrichment scripts into a real
`nf_NovInvenio` pipeline step for the Coccidioides pangenome study (and
future studies generally). Explicitly deferred by the user ("I think is
an exploration for later") — not part of that design's initial scope.
