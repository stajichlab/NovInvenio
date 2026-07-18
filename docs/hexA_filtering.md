# Why hexA and hex-1 don't group into one family (pezizo5_fungi)

**Short answer:** hexA (`Afu5g08830`, *A. fumigatus*) and hex-1 (`NCU08332`,
*N. crassa*) are true homologs and match each other superbly in the raw search
(e-value **6.6e-69**). They still fail to become a family because the
**paralog-competition filter** (filter 2 in `bin/build_presence_matrix.py`)
discards their cross-proteome presence calls. HEX-1 is an evolutionary
derivative of eIF5A, and in most genomes the *real* eIF5A ortholog hits the
other genome slightly harder than the HEX-1 ortholog does — so the filter
attributes the hit to the shared eIF5A domain and zeroes it out. That drops each
protein below the ingroup presence threshold, so neither is ever written to
`candidates.txt`, and mmseqs only clusters candidates. The grouping step never
sees them.

## Cast

| Role | Afum | Ncra |
|---|---|---|
| HEX-1 (Woronin body protein) | `Afu5g08830-T-p1` (hexA) | `NCU08332-t26_1-p1` (hex-1) |
| eIF5A paralog (ancestral fold) | `Afu1g04070-T-p1` | `NCU05274-t26_1-p1` |

Config `configs/pezizo5_fungi.csv`: ingroup = {Amega, Ncra, Afum, Ztri, Cimm}
(5 species), `ingroup_min_frac = 0.75` (need ≥ 4 of 5).

## The two filters (`bin/build_presence_matrix.py`)

A raw hit becomes a "present" call only if it passes **both**:

1. **Paralog-cutoff filter** — hit e-value must beat the query protein's own
   within-genome paralog e-value (from the self-vs-self search).
2. **Paralog-competition filter** — if the query's paralog hits the *same target
   proteome* with a better (lower) e-value than the query itself does, the hit is
   disqualified. Rationale: such a hit is better explained by the conserved
   domain shared with the paralog than by the query protein.

### Filter 1 passes easily

Self-search paralog cutoffs (`self_hits/*.paralog_cutoffs.tsv`):

- `Afu5g08830` → paralog `Afu1g04070`, cutoff **4.2e-11**
- `NCU08332` → paralog `NCU05274`, cutoff **8.0e-09**

The hexA↔hex-1 hit (6.6e-69 / 7.8e-70) is orders of magnitude better than either
cutoff, so filter 1 keeps it.

### Filter 2 kills the cross calls

Filter 2 compares, for each `(query, target proteome)`, the query's best hit
against its eIF5A paralog's best hit into the same proteome:

| Direction | HEX-1 ortholog hit | eIF5A paralog hit into same proteome | Verdict |
|---|---|---|---|
| hexA → **Ncra** | 6.6e-69 (→ hex-1) | `Afu1g04070` → Ncra = **8.7e-70** | eIF5A wins → **dropped** |
| hex-1 → **Afum** | 7.8e-70 (→ hexA) | `NCU05274` → Afum = **7.5e-70** | eIF5A wins (barely) → **dropped** |
| hexA → **Amega** | 8.6e-49 | `Afu1g04070` → Amega = **3.2e-74** | eIF5A wins → **dropped** |
| hexA → **Cimm** | **6.7e-152** | `Afu1g04070` → Cimm = 2.8e-82 | hexA wins → **kept** |
| hexA → **Ztri** | — | — | hexA wins → kept |

The Afum↔Ncra reciprocal miss is decided by margins as thin as 7.5e-70 vs
7.8e-70 — effectively a tie between the HEX-1 and eIF5A signals. Cimm survives
only because its HEX-1 copy is an unusually strong match (6.7e-152) that
out-scores even the highly conserved eIF5A.

## Downstream consequence

After filter 2, the presence-matrix rows are:

```
protein_id           source  Afum Amega Ccin Cimm Cneo Mcir Ncra Nirr Scer Spom Ztri
Afu5g08830-T-p1  (hexA)  Afum   1    0    0    1    0    0    0    0    0    0    1
NCU08332-t26_1-p1 (hex-1) Ncra   0    0    0    1    0    0    1    0    0    0    1
```

Each protein registers presence in only **3 of 5** ingroup species
(3/5 = 0.60 < 0.75). Both fail `ingroup_min_frac`, so:

1. Neither is written to `candidates.txt`.
2. `EXTRACT_CANDIDATES` never pulls them into `candidates.fa`.
3. `MMSEQS_CLUSTER` only clusters `candidates.fa` → the two proteins are never
   presented to clustering and **cannot** be placed in a gene family.

The "why don't they group into a family" question is therefore downstream of the
real cause: they are filtered out *before* clustering. Even the family machinery
aside, both are correctly non-novelties — HEX-1 is a Pezizomycotina-wide Woronin
body protein, not lineage-specific — but the run under-counts their true ingroup
breadth (should be ~5/5) because the eIF5A-derived signal trips the competition
filter.

## Levers if you want these captured

- **Lower `ingroup_min_frac`** (e.g. 0.6) — blunt; admits many weak candidates.
- **Relax filter 2** so a paralog only disqualifies a hit when it beats it by a
  meaningful margin (e.g. paralog e-value must be lower by ≥ N orders of
  magnitude), instead of the current strict `<`. This is the targeted fix: the
  Afum↔Ncra misses are near-ties, and a margin requirement would keep them.
- **Filter 2 by target identity, not just target proteome** — disqualify only
  when the paralog out-scores the query *on the same target protein*. Here the
  eIF5A paralog wins on a *different* target (the genome's eIF5A gene), while
  hexA is still the best hit to the HEX-1 gene. Requiring the competition to be
  on the same target would preserve the HEX-1 call.

Both filter changes touch `bin/build_presence_matrix.py` lines 132–141 and would
need regression checks against the existing candidate counts.
