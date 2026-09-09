<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-09-08

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-20 | Implementation workflow: tickets → branch → PR, mycelium Stop hook (`mycelium-stop-check.sh`) is disabled |
| decisions.md | 2 entries (large — read selectively) | 2026-09-08 | ADR-0002 grilling resolutions (running — folds into the ADR when complete), Report presentation: one skin registry, one linkout builder, one landing-page design |
| learnings.md | 12 entries (large — read selectively) | 2026-09-08 | Claude Code owns settings.local.json — put hooks in settings.json, mycelium Stop hook false-blocks interactive (grilling) sessions, famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes, mycelium generate_index.py needs Python 3.10+ — silently fails under system python3 (3.9), Session-log stubs vs substance: hooks create logs, agents must write the findings |
| log/ | 30 sessions | 2026-09-08 | novinvenio (30) |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-09-08 (heuristic)

## Tag clusters

- **mycelium** (4 entries) — L-1, L-2, L-4, L-5
- **nextflow** (4 entries) — L-7, L-8, L-10, L-11
- **slurm** (4 entries) — L-3, L-8, L-10, L-11
- **mmseqs** (3 entries) — L-6, L-12, D-1
- **tooling** (3 entries) — L-1, L-2, L-4
- **accessibility** (2 entries) — L-9, D-2

## Most recent (10)

- [2026-09-08] L-12: diamond cluster's UniProt ID handling and CLI quirks (vs mmseqs2)
- [2026-09-06] L-11: `task.index` in an output filename silently defeats `-resume` caching across separate invocations
- [2026-09-05] L-8: Two distinct BUILD_CHUNK "failure" signatures on the broader-grid sweeps — one benign, one from manual queue intervention
- [2026-09-05] L-9: The report pages' colour duplication was hiding two real defects
- [2026-09-05] L-10: A live-checkout `-profile slurm` job reads whatever `bin/`/`lib/` is on disk *right now*, not what it started with
- [2026-09-05] D-2: Report presentation: one skin registry, one linkout builder, one landing-page design
- [2026-09-04] L-7: Fixing a bin/ script doesn't bust Nextflow's `-resume` cache by itself
- [2026-08-06] L-4: mycelium generate_index.py needs Python 3.10+ — silently fails under system python3 (3.9)
- [2026-08-06] L-5: Session-log stubs vs substance: hooks create logs, agents must write the findings
- [2026-08-06] L-6: Investigation-run observations (benchmarks & cross-method divergences)

## By tag

- `mycelium`: L-1, L-2, L-4, L-5
- `nextflow`: L-7, L-8, L-10, L-11
- `slurm`: L-3, L-8, L-10, L-11
- `mmseqs`: L-6, L-12, D-1
- `tooling`: L-1, L-2, L-4
- `accessibility`: L-9, D-2
- `caching`: L-7, L-11
- `colour-blind`: L-9, D-2
- `findings`: L-5, L-6
- `gotcha`: L-7, L-11
- `grilling`: L-2, D-1
- `hooks`: L-1, L-2
- `linkouts`: L-9, D-2
- `pfam`: L-6, L-9
- `report`: L-9, D-2
- `resume`: L-7, L-11
- `skins`: L-9, D-2
- `wcag`: L-9, D-2
- `adr-0002`: D-1
- `avx2`: L-3
- `benchmarking`: L-12
- `bin-scripts`: L-7
- `bioconda`: L-3
- `claude-code`: L-1
- `cluster`: L-3
- `cluster-tool`: L-12
- `concordance`: L-6
- `configuration`: L-1
- `container`: L-10
- `correction`: L-10
- `debounce`: L-2
- `diamond`: L-12
- `docker-build`: L-10
- `duplication`: L-9
- `environment`: L-4
- `exitReadTimeout`: L-8
- `false-failure`: L-8
- `famsa`: L-3
- `gene-family`: D-1
- `git`: L-10
- `github-pages`: D-2
- `gpfs`: L-8
- `hardware-compatibility`: L-3
- `hmmer`: D-1
- `hmmscan`: L-6
- `id-normalization`: L-12
- `in-progress`: D-1
- `index`: L-4
- `interactive-session`: L-2
- `knowledge-capture`: L-5
- `logging`: L-5
- `loss`: L-6
- `novelty`: L-6
- `pairwise`: L-6
- `param-sweep`: L-11
- `pixi`: L-3
- `presence-calling`: L-7
- `process`: L-5
- `python`: L-4
- `queue-management`: L-8
- `race-condition`: L-10
- `scatter-gather`: L-11
- `sessions`: L-5
- `settings`: L-1
- `silent-failure`: L-4
- `simd`: L-3
- `singleton`: L-6
- `stop-hook`: L-2
- `task-hash`: L-7
- `uniprot`: L-12
- `version`: L-4

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 "$(cat .mycelium/plugin-root)/skills/core/scripts/recall_lessons.py" --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
