<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-21

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-20 | Implementation workflow: tickets → branch → PR, mycelium Stop hook (`mycelium-stop-check.sh`) is disabled |
| decisions.md | 1 entry | 2026-07-21 | ADR-0002 grilling resolutions (running — folds into the ADR when complete) |
| learnings.md | 3 entries | 2026-07-21 | Claude Code owns settings.local.json — put hooks in settings.json, mycelium Stop hook false-blocks interactive (grilling) sessions, famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-21 (heuristic)

## Tag clusters

- **grilling** (2 entries) — L-2, D-1
- **hooks** (2 entries) — L-1, L-2
- **mycelium** (2 entries) — L-1, L-2
- **tooling** (2 entries) — L-1, L-2

## Most recent (10)

- [2026-07-21] L-3: famsa (bioconda) requires AVX2 — SIGILLs on the cluster's Abu Dhabi nodes
- [2026-07-20] L-1: Claude Code owns settings.local.json — put hooks in settings.json
- [2026-07-20] L-2: mycelium Stop hook false-blocks interactive (grilling) sessions
- [2026-07-20] D-1: ADR-0002 grilling resolutions (running — folds into the ADR when complete)

## By tag

- `grilling`: L-2, D-1
- `hooks`: L-1, L-2
- `mycelium`: L-1, L-2
- `tooling`: L-1, L-2
- `adr-0002`: D-1
- `avx2`: L-3
- `bioconda`: L-3
- `claude-code`: L-1
- `cluster`: L-3
- `configuration`: L-1
- `debounce`: L-2
- `famsa`: L-3
- `gene-family`: D-1
- `hardware-compatibility`: L-3
- `hmmer`: D-1
- `in-progress`: D-1
- `interactive-session`: L-2
- `mmseqs`: D-1
- `pixi`: L-3
- `settings`: L-1
- `simd`: L-3
- `slurm`: L-3
- `stop-hook`: L-2

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
