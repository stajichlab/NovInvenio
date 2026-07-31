# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

<!-- Document any project-specific convention overrides here. -->

## Implementation workflow: tickets → branch → PR

When moving from planning to implementation, do not edit code on `main`. Turn the work
into GitHub issues first (one issue = one PR-sized change, sourced from `todo/` items +
the relevant `docs/adr/`), branch per ticket (`NNN-short-slug` off `main`), implement with
small commits referencing `#NNN`, open a PR (`Closes #NNN`, links the ADR/todo), and on
merge flip the `todo/TODO_REGISTRY.md` row to `complete`. Only create issues/branches/PRs
once the user has asked to start implementing; confirm before the first outward-facing
action. Canonical text lives in `CLAUDE.md` → Development Guide → "Implementation workflow
(tickets → branch → PR)"; the `gh` mechanics are in `docs/agents/issue-tracker.md`.

**Source**: requested 2026-07-20; the user referred to this as the "to-tickets" flow (not
a defined slash command — it is a documented workflow).

**Tags**: workflow, git, github-issues, pull-request, process, implementation

## mycelium Stop hook (`mycelium-stop-check.sh`) is disabled

`.claude/settings.json`'s `Stop` array no longer registers `mycelium-stop-check.sh`; only
the non-blocking `mycelium-data-lineage-stop.sh` remains. Disabled because the hook
repeatedly false-blocked interactive/multi-turn sessions (grilling, back-and-forth Q&A)
even when `.living/decisions.md`/`learnings.md` were demonstrably updated after its own
debounce timestamp — see `[[stop-hook-false-blocks-interactive-sessions]]` in
`.living/learnings.md` for the diagnosis. There is no reliable signal inside a `Stop` hook
to distinguish "interactive session" from "batch session" (it only sees file mtimes and a
debounce timestamp), so this is a **blanket** disable, not a conditional exemption, despite
being requested as one — the honest tradeoff, not a workaround.

**Consequence**: automatic enforcement of the "update `.living/` after significant work"
protocol is off. The protocol itself is unchanged — still follow it every session — this
just removes the hard stop that used to catch a skipped one.

**Reinstate when**: the upstream hook gets a real interactive-session signal (e.g. an
`active_skill` field in the Stop payload, or keying "reflected" off any post-session-start
`.living/*` mtime instead of the Bash-re-armed reminder), or this repo decides the
false-block nuisance is worth tolerating for the enforcement.

**Source**: requested 2026-07-20, after the hook false-blocked several turns of the
ADR-0002 grilling session (see conversation + `.living/decisions.md` grilling entries).

**Tags**: mycelium, hooks, stop-hook, tooling, workflow
