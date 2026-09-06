# Species trait data

Supporting input for the targeted config-builder (see
`docs/superpowers/specs/2026-09-05-config-builder-design.md`) — not a
NovInvenio run config, never passed directly as `--config` to `main.nf`.

## Files

- `trait_definitions.yaml` — the controlled vocabulary. Each trait declares
  its own fixed set of legal values. A value that isn't declared here is
  a hard error when `traits.csv` is loaded, not a near-miss best guess.
- `traits.csv` — the actual species-to-trait assignments, one row per
  (species, trait) pair: `Species,trait,value,source,notes`.
- `funguild_reference.csv` — a raw reference dump, **not** curated
  input to `traits.csv` yet. Pulled 2026-09-05 from the `funguild`
  table in `/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb`
  (a real FUNGuild run over the Fungi_5k assembly set), deduplicated to
  one row per distinct `(taxon, growthForm, guild, trophicMode,
  confidenceRanking)` combination with an `n_assemblies` count of how
  many assemblies support that call (1,871 rows; no genus had
  conflicting calls across its assemblies as of this pull). `taxon` is
  almost always genus-level, occasionally a full binomial — FUNGuild's
  own reference database resolves at genus level in most cases, so
  treat any imported trait as a genus-level inference unless verified
  otherwise, not a species-specific fact. Whether/how this gets mapped
  into `trait_definitions.yaml`'s controlled vocabulary (FUNGuild's
  `trophicMode` is only a 3-way Pathotroph/Saprotroph/Symbiotroph
  split, coarser than a plant-pathology-style
  biotroph/necrotroph/hemibiotroph distinction) is still an open
  question — see "Not yet decided" below.

## Adding data

1. **Using an existing trait** on a new species: add a row to
   `traits.csv` using a `value` already declared for that `trait` in
   `trait_definitions.yaml`. `source` is optional (a DOI, a paper
   citation, or blank for "known to the lab"); `notes` is free text.
2. **A new value for an existing trait**: add it to that trait's
   `values` map in `trait_definitions.yaml` first, then use it in
   `traits.csv`.
3. **A brand new trait**: add a new top-level key under `traits:` in
   `trait_definitions.yaml` (with a `description` and its `values` map)
   before referencing it in any `traits.csv` row.

Keep trait names and values as short lowercase snake_case tokens
(matches the rest of this repo's `--Short`/config conventions) —
`spore_motility: motile`, not `"Spore Motility": "Motile"`.

## Scope for now

This is a small, hand-curated seed intentionally scoped to the species
relevant to the first targeted-study batch (Mucor circinelloides,
Phycomyces blakesleeanus, Basidiobolus meristosporus, Spizellomyces
punctatus, and their near neighbors) — not an attempt at a general
fungal trait database. Grow it study-by-study as new comparisons need
new traits or new species, rather than front-loading traits that aren't
in use yet.

`ontology_term` in `trait_definitions.yaml` is optional provenance
citing a matching term from an established phenotype/anatomy ontology
(PATO, the Fungal Anatomy Ontology, GO, FYPO, OMP) — populate it only
when you've actually looked up the real term ID; leave it blank rather
than guess. Nothing in v1 validates against a live ontology or depends
on this field being filled in.

## Reliability caveat: trophically diverse genera

A genus-level FUNGuild call is weakest exactly where a genus spans very
different ecological strategies across its species — flagged explicitly
in `traits.csv`'s `notes` (search for "CAUTION") for `Aspergillus`,
`Fusarium`, `Alternaria`, `Penicillium`, `Trichoderma`,
`Colletotrichum`, `Cladosporium`, `Curvularia`, `Bipolaris`, `Candida`,
`Talaromyces`, `Chaetomium`, and `Pyrenophora`. E.g. `Aspergillus
nidulans`'s genus-level composite guild (drawn from the whole genus,
which also contains the very different *A. fumigatus*) should not be
trusted for a real study without species-level verification. Treat
these rows as a starting point to check, not settled fact.

## Known vocabulary gaps (not yet fixed)

Found while populating real data — flagged rather than force-fit into
an existing value:

- **Carotenoid pigmentation** (Mucor, Phycomyces) isn't representable —
  `pigmentation`'s `melanized`/`non_melanized` values only cover melanin.
  Both species are recorded `non_melanized` (their hyphae are hyaline)
  with a note that their carotenoid-based sporangial colour isn't
  captured by any current value.
- **Glomerospores** (Rhizophagus irregularis's AM-fungal spore type) fit
  none of `asexual_reproduction`'s values (`conidia`/`sporangiospores`/
  `budding`/`fragmentation`) — deliberately left unassessed for that
  species rather than stretched into a wrong bucket.
- **Chytrid zoospore release** was stretched onto `sporangiospores` for
  Spizellomyces punctatus and Batrachochytrium dendrobatidis (both
  produce motile zoospores from a zoosporangium, structurally closer to
  "enclosed in a sporangium" than to conidia/budding/fragmentation) —
  flagged in both rows' notes as a judgment call, not a clean match.

## How a trait drives species selection

Resolved 2026-09-05 — see the "Trait mode" section of the design spec
for the full reasoning. Short version: `mode: trait` in a study spec is
`mode: nearest`, filtered — it takes the same lineage-scoped candidate
pool `nearest` already uses, keeps only candidates with the exact
`(trait, value)` given, then ranks the survivors by the same
phylogenetic-proximity function `nearest` uses. Trait matching filters;
proximity orders. It never overrides phylogeny outright, and a trait
shared only with something phylogenetically distant (the original
motivating example — an early-diverging fungus and an animal both
having motile zoospores) is out of scope for the automatic mode and
stays a manual `mode: explicit` pick with a `reason:` string citing the
relevant row(s) here.
