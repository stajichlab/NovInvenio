# Rename novelty_discovery/novelty_screen GROUP labels for clarity

| Field | Value |
|-------|-------|
| **Date** | 2026-07-28 |
| **Author** | Jason Stajich |
| **Priority** | low |
| **Status** | complete |
| **Category** | chore |
| **Related** | `todo/novelty-discovery-screen.md`, issues #24-#29 |

## Idea

The four `GROUP` values introduced for the two-phase targeted novelty pipeline
(`--cluster_tool novelty_discovery`) are terse and could be more self-explanatory:

| Current | Proposed | Role |
|---|---|---|
| `TARGET` | `DISCOVERY_TARGET` | The genome(s) to find novelties in |
| `DISC_OUT` | `DISCOVERY_OUT` | Small reference panel for the phase-1 absence call |
| `NEAR_IN` | `NEAR_INGROUP` | Close relatives of the target, same clade |
| `BROAD_OUT` | `BROAD_OUTGROUP` | Distant lineages, outside the target clade |

Requested as aliases or a rename — either accept both old and new spellings during a
transition, or rename outright with existing configs updated.

## Why this matters

The current names are ambiguous out of context — `TARGET` and `NEAR_IN` in particular
don't read as obviously novelty_discovery-specific roles next to the classic `IN`/`OUT`
labels, and `DISC_OUT` isn't an obvious abbreviation on first read.

## Where this touches

- `lib/config_parser.py` — `GROUPS`, `INGROUP_ROLES`, `OUTGROUP_ROLES`,
  `get_target()`/`get_disc_out()`/`get_near_in()`/`get_broad_out()`
- `main.nf` — the `samples_ch.filter { meta.group == 'TARGET' }`-style channel filters
- `bin/novelty_presence_matrix.py`, `bin/novelty_screen.py` — group-string comparisons
- `configs/novelty_discovery_example.csv`, `configs/sordario.csv` — existing configs using
  the current labels
- `CLAUDE.md`'s "Two-Phase Targeted Novelty Pipeline" section and
  `todo/novelty-discovery-screen.md`'s "Config Format" table — documentation

## Status

Complete — implemented as a rename with backward-compatible aliases. `lib/config_parser.py`
now normalizes `TARGET`/`DISC_OUT`/`NEAR_IN`/`BROAD_OUT` (via `GROUP_ALIASES`) to the
canonical `DISCOVERY_TARGET`/`DISCOVERY_OUT`/`NEAR_INGROUP`/`BROAD_OUTGROUP` spellings in
`parse_config()`; `main.nf` mirrors the same alias table for its own Groovy-side CSV
parsing. Existing config CSVs (e.g. `configs/sordario.csv`) using the old labels keep
working unchanged. `get_target()`/`get_disc_out()`/`get_near_in()`/`get_broad_out()` were
renamed to `get_discovery_target()`/`get_discovery_out()`/`get_near_ingroup()`/
`get_broad_outgroup()` to match.
