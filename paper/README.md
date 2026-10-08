# paper/: working folder for the NovInvenio manuscript

Nothing here is a finished manuscript. These are drafts to be checked and rewritten.

| File | What it is |
|---|---|
| `OUTLINE.md` | Section plan, claims, planned figures, and for each result where its evidence is and whether it must be rerun. |
| `methods_pairwise.md` | Methods draft for the pairwise novelty pipeline (novelty and loss directions), written from the code. |
| `methods_pangenome.md` | Methods draft for the pangenome pipeline, written from the code. |
| `flow_pairwise.md`, `flow_pangenome.md` | Information flow: Mermaid charts and one table row per step (inputs, outputs, what is filtered or tested). |
| `OPEN_ISSUES_FOR_PAPER.md` | Problems found in review that change what the paper can claim, each with a status. Read this before writing results. |

## How the methods drafts were made

Two language-model agents read the code at commit `3b045c6` (`main` on 2026-10-08) and wrote the
methods and flow files. They were told to trace each statement to code and to tag anything they
did not read directly with `[UNVERIFIED]`. A maintainer must still check every number and rule
against the code before it goes into a manuscript. Treat the drafts as a map of the code, not as
evidence.

Known state of the code at that commit that the drafts may describe:
- Fixes from PR #224 (`--evalue` reaching filter 1, wording) and PR #221/#223 (labels, losses signal) may not be merged yet.
  Check the discrepancy sections of each methods file against `main`.
- The `mmseqs` and `novelty_discovery` pathways are scheduled for removal (issue #219) and are not described.

## Conventions

- Plain, short sentences. One idea per sentence. No idioms.
- Numbers come from a named file. Counts of proteins are not counts of gene families.
- "Absent" means "no qualifying hit", never proof of absence.
- Re-run the pipeline that produced a published number before citing it.

## Updating

When a PR changes behaviour described in a methods file, update the file in the same PR and
move the row in `OPEN_ISSUES_FOR_PAPER.md`.
