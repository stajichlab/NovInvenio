<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-07

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-20 | Implementation workflow: tickets → branch → PR, mycelium Stop hook (`mycelium-stop-check.sh`) is disabled |
| decisions.md | 1 entry | 2026-08-07 | ADR-0002 grilling resolutions (running — folds into the ADR when complete) |
| learnings.md | 6 entries | 2026-08-07 | Claude Code owns settings.local.json — put hooks in settings.json, mycelium Stop hook false-blocks interactive (grilling) sessions, famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes, mycelium generate_index.py needs Python 3.10+ — silently fails under system python3 (3.9), Session-log stubs vs substance: hooks create logs, agents must write the findings |
| log/ | 10 sessions | 2026-08-07 | novinvenio (10) |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-07 (heuristic)

## Tag clusters

- **mycelium** (4 entries) — L-1, L-2, L-4, L-5
- **tooling** (3 entries) — L-1, L-2, L-4
- **findings** (2 entries) — L-5, L-6
- **grilling** (2 entries) — L-2, D-1
- **hooks** (2 entries) — L-1, L-2
- **mmseqs** (2 entries) — L-6, D-1

## Most recent (10)

- [2026-08-06] L-4: mycelium generate_index.py needs Python 3.10+ — silently fails under system python3 (3.9)
- [2026-08-06] L-5: Session-log stubs vs substance: hooks create logs, agents must write the findings
- [2026-08-06] L-6: Investigation-run observations (benchmarks & cross-method divergences)
- [2026-07-21] L-3: famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes
- [2026-07-20] L-1: Claude Code owns settings.local.json — put hooks in settings.json
- [2026-07-20] L-2: mycelium Stop hook false-blocks interactive (grilling) sessions
- [2026-07-20] D-1: ADR-0002 grilling resolutions (running — folds into the ADR when complete)

## By tag

- `mycelium`: L-1, L-2, L-4, L-5
- `tooling`: L-1, L-2, L-4
- `findings`: L-5, L-6
- `grilling`: L-2, D-1
- `hooks`: L-1, L-2
- `mmseqs`: L-6, D-1
- `adr-0002`: D-1
- `avx2`: L-3
- `bioconda`: L-3
- `claude-code`: L-1
- `cluster`: L-3
- `concordance`: L-6
- `configuration`: L-1
- `debounce`: L-2
- `environment`: L-4
- `famsa`: L-3
- `gene-family`: D-1
- `hardware-compatibility`: L-3
- `hmmer`: D-1
- `hmmscan`: L-6
- `in-progress`: D-1
- `index`: L-4
- `interactive-session`: L-2
- `knowledge-capture`: L-5
- `logging`: L-5
- `loss`: L-6
- `novelty`: L-6
- `pairwise`: L-6
- `pfam`: L-6
- `pixi`: L-3
- `process`: L-5
- `python`: L-4
- `sessions`: L-5
- `settings`: L-1
- `silent-failure`: L-4
- `simd`: L-3
- `singleton`: L-6
- `slurm`: L-3
- `stop-hook`: L-2
- `version`: L-4

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
