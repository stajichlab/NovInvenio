# Choosing in/out groups: the framework, and this directory's studies

Each `.yaml` file here is a batch spec for `bin/build_targeted_configs.py`
(see `docs/superpowers/specs/2026-09-05-config-builder-design.md` for the
full design). This file documents two things: the decision process for
picking a focal species' ingroup companions and outgroup pool, and the
specific hypothesis each study in this directory tests.

## The framework: how to choose an ingroup and outgroup

1. **Pick the focal species and the question.** What biology are you
   actually asking about — a shared trait, a lineage's overall novelty,
   or a phenotype that cuts across phylogeny?

2. **Pick the companion-selection mode based on where the shared biology
   actually sits in the tree:**

   - **`mode: nearest`** — the trait (or "the whole lineage") tracks
     phylogeny. Picks the closest real relatives at a given `scope_rank`
     (default `ORDER`; narrow to `FAMILY` for a tighter test, widen if the
     family is too sparse — see the "real data is sparse" note below).
     This is a **lineage-novelty** test: "what's specific to this species'
     immediate neighborhood," not "what's specific to a trait."
   - **`mode: trait`** — the trait recurs among *close* relatives (same
     `scope_rank`-scoped pool `mode: nearest` would use), and you have
     real `config_support/traits/` data for more than one candidate.
     Filters to the trait, then ranks by proximity — see
     `dimorph_onygenales.yaml` for a worked example (dimorphism recurs
     within Onygenales). This does **not** search the whole pool; a trait
     shared only with something phylogenetically distant needs
     `mode: explicit` instead (see below), or several independent
     `mode: nearest` studies compared post-hoc (see the salt-tolerance
     studies).
   - **`mode: explicit`** — the companions you want aren't phylogenetic
     neighbors at all (a manual override, or a genuinely cross-kingdom
     comparison via `--extra-pool`). Always give a `reason:` — it's the
     only place "why" gets recorded for this mode.

3. **Pick the outgroup pool to match the actual contrast you want**, not
   just "the rest of Dikarya" by default:
   - A broad Dikarya reference panel (this directory's shared
     `dikarya_reference_v1`: *Saccharomyces cerevisiae*, *Neurospora
     crassa*, *Aspergillus nidulans*, *Schizosaccharomyces pombe*,
     *Coprinopsis cinerea*) answers "is this novel across most of Fungi."
   - A **trait-matched or morphology-matched outgroup** (see
     `zygo_filamentous.yaml`'s filamentous-only panel, or
     `dimorph_onygenales.yaml`'s non-dimorphic same-order dermatophytes)
     controls for a confound and answers a sharper question: "novel
     relative to fungi that share this one property but not the trait
     under study."
   - A **different early-diverging lineage as the outgroup** (see
     `zoosporic_zygomycete.yaml`) asks "is this specific to my focal
     lineage, or shared broadly across all early-diverging fungi."

4. **When a trait is convergent across distant lineages** (this
   directory's three `salt_*.yaml` studies), `mode: trait` can't unite
   them in one study by design. Run one independent `mode: nearest` study
   per lineage against the *same* outgroup pool, then compare the
   resulting novelty-candidate lists for overlapping Pfam domains or
   functional annotations post-hoc — the convergence signal is in the
   comparison across studies, not in any single study.

5. **Real data is sparse for early-diverging fungi.** A focal species'
   family or even order may have only one or zero other representatives
   with `repr_assignments.tsv` coverage (`Batrachochytrium
   dendrobatidis`'s own family has zero; its order has exactly one). An
   `n:` target is a maximum, not a guarantee — `mode: nearest`/`trait`
   silently return however many real candidates exist. Check first (a
   quick `grep`/`csv` scan of `config_support/master_pool.csv`'s `Lineage`
   column for the rank you're scoping to) rather than assuming `n` will
   be met.

## The studies in this directory

| Batch | Focal(s) | Mode | Outgroup | Hypothesis |
|---|---|---|---|---|
| `dimorph_onygenales.yaml` | Coccidioides immitis | `trait` (`cell_morphology=dimorphic`, scope `ORDER`) | non-dimorphic Onygenales dermatophytes | Genes correlated with the evolution of thermal dimorphism, using same-order non-dimorphic relatives as the sharpest possible negative control. |
| `pathogen_cryptococcus.yaml` | Cryptococcus neoformans | `nearest` (scope `FAMILY`) | broad Dikarya | Candidate virulence-factor genes specific to the pathogenic Cryptococcus/Cryptococcaceae clade. |
| `pathogen_candida_auris.yaml` | Candidozyma auris | `nearest` (scope `FAMILY`) | broad Dikarya | Genes specific to the emerging multidrug-resistant Candidozyma/Metschnikowiaceae pathogen clade. |
| `salt_hortaea.yaml` | Hortaea werneckii | `nearest` (scope `FAMILY`) | broad Dikarya | Replicate 1/3 of the convergent-halotolerance test (Dothideomycetes lineage). |
| `salt_debaryomyces.yaml` | Debaryomyces hansenii | `nearest` (scope `FAMILY`) | broad Dikarya | Replicate 2/3 (Saccharomycotina lineage) — note its own family-mates include non-halotolerant Candida species; that's expected for `mode: nearest`. |
| `salt_rhodotorula.yaml` | Rhodotorula mucilaginosa | `nearest` (scope `FAMILY`) | broad Dikarya | Replicate 3/3 (Pucciniomycotina/Basidiomycota lineage) — phylogenetically farthest from the other two, the strongest convergence test if all three overlap. |
| `zygo_filamentous.yaml` | Mucor circinelloides, Phycomyces blakesleeanus | `nearest` (scope `FAMILY`/`ORDER`) | filamentous-only Dikarya (excludes yeasts and dimorphic Cryptococcus) | Genes specific to Mucoromycotina, controlling for gross cell morphology so the signal isn't just "filamentous vs. yeast." |
| `zoosporic_dikarya.yaml` | Batrachochytrium dendrobatidis | `nearest` | broad Dikarya | Baseline early-diverging-fungi-vs-crown-fungi contrast. |
| `zoosporic_zygomycete.yaml` | Batrachochytrium dendrobatidis | `nearest` | Mucor, Phycomyces, Basidiobolus | Separates chytrid-*specific* genes from genes shared broadly across all early-diverging (non-Dikarya) fungi — compare against `zoosporic_dikarya.yaml`. |
| `zoosporic_opisthokont_dikarya.yaml` | Batrachochytrium dendrobatidis | `explicit` (+ `--extra-pool config_support/animal_pool.csv`) | broad Dikarya | The original motivating question: genes shared between zoosporic fungi and animals/their closest non-animal relative (Monosiga brevicollis), absent from Dikarya — candidate ancestral-opisthokont genes lost along the Dikarya stem. |
| `mucoromycota_focal_v1.yaml` | Mucor, Phycomyces, Basidiobolus, Batrachochytrium | `nearest`/`explicit` | broad Dikarya | The original worked example from the design spec — one study per early-diverging lineage vs. a shared reference panel. |

Every batch above was run against the real `config_support/master_pool.csv`
(built from `Fungi_BFD`) and, where relevant, `config_support/animal_pool.csv`
(via `--extra-pool`) to confirm it renders a real config before being
committed here — see each file's own header comment for render-time details
(species substitutions forced by real data availability, scope-rank choices,
etc.).
