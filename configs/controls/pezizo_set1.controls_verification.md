# pezizo_set1 positive-control verification (2026-09-10)

Companion to `pezizo_set1.controls.csv`. Records the empirical check behind each row —
real `diamond blastp --very-sensitive` search of each Neurospora crassa (Ncra) anchor
protein against all 11 `pezizo_set1` proteome databases (`results/pezizo_set1/search_cache/*.dmnd`),
not a literature claim taken on faith. Gene list originally supplied by the user
(Jason Stajich) from Neurospora genetics literature; NCU IDs resolved via
`config_support/modelorgs/Neurospora_crassa_gene_names_FungiDB.csv` (NII repo) where the
UniProt `GN=` tag didn't already carry the common name.

Ingroup = Amega, Ncra, Afum, Ztri, Cimm (Pezizomycotina). Outgroup = Nirr, CneoH99, Ccin,
Mcir, Scer, Spom. `ingroup_min_frac` default is 0.75 (>=4/5).

## Result table (best e-value per species, `-` = no hit even at e<=1)

| gene | Amega | Ncra | Afum | Ztri | Cimm | Nirr | CneoH99 | Ccin | Mcir | Scer | Spom |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hex-1 | 4.6e-52 | 4.7e-121 | 1.6e-70 | 8.9e-83 | 6.8e-74 | 1.2e-07 | 3.2e-10 | 1.5e-11 | 3.4e-12 | 1.7e-13 | 7.6e-12 |
| wsc | 3.6e-79 | 5.0e-213 | 3.7e-88 | 1.0e-88 | 4.2e-92 | - | 3.5e-03 | 4.1e-48 | 4.4e-46 | 2.3e-06 | 3.2e-03 |
| lah (NCU02793) | 1.1e-59 | 0.0 | 4.1e-111 | 1.6e-114 | 4.3e-117 | - | - | - | - | - | - |
| so/soft | 1.1e-300 | 0.0 | 0.0 | 0.0 | 0.0 | 6.5e-04 | 4.4e-41 | 3.5e-54 | 2.4e-58 | 1.9e-41 | - |
| ada-1 | 3.7e-83 | 0.0 | 1.9e-73 | 1.6e-62 | 4.9e-77 | - | 1.8e-01 | - | 5.7e-02 | - | - |
| ham-5 | 4.3e-134 | 0.0 | 5.7e-241 | 3.0e-238 | 1.6e-208 | - | - | 4.8e-01 | - | 2.9e-01 | - |
| ham-8 | 6.3e-71 | 0.0 | 8.0e-148 | 1.0e-134 | 2.3e-133 | - | - | - | - | - | - |
| spa-1 | 7.4e-12 | 7.8e-303 | 7.4e-17 | 8.9e-49 | 2.5e-37 | - | - | - | - | - | - |
| spa-9 | - | 0.0 | 2.0e-42 | 5.4e-01 | 1.6e-42 | - | - | - | - | - | - |
| ham-11 | - | 0.0 | 1.3e-06 | - | 3.9e-06 | - | - | - | - | - | - |
| spa-18 | 6.8e-14 | 0.0 | 8.1e-38 | - | 8.5e-29 | 6.8e-08 | - | - | - | - | - |

## Verdicts

| gene | ingroup frac | outgroup contamination | verdict |
|---|---|---|---|
| hex-1 | 5/5 | weak hit in all 6 (eIF-5A paralog cross-reactivity) | in `controls.csv`, flagged: only trust once the paralog-competition filter (post `--very-sensitive` self-search fix, 2026-09-10) is confirmed disqualifying these |
| lah | 5/5 | 0/6 | in `controls.csv`, clean |
| ada-1 | 5/5 | noise-level in 2/6 (e-01, e-02) | in `controls.csv`, clean |
| ham-5 | 5/5 | noise-level in 2/6 (e-01) | in `controls.csv`, clean |
| ham-8 | 5/5 | 0/6 | in `controls.csv`, clean |
| spa-1 | 5/5 | 0/6 | in `controls.csv`, clean |
| spa-18 | 4/5 (missing Ztri) | weak hit in Nirr (6.8e-08) | in `controls.csv`, flagged (same paralog-risk class as hex-1) |
| spa-9 | 3/5 (no hit Amega, marginal Ztri) | 0/6 | **excluded** — fails `ingroup_min_frac=0.75` for this exact species set regardless of broader Pezizomycotina literature support |
| ham-11 | 3/5, and weak (e-06) where present | 0/6 | **excluded**, same reason |
| wsc | 5/5 | **strong** hits in Ccin (4.1e-48) and Mcir (4.4e-46) | **excluded from controls.csv** — contradicts clean lineage restriction; likely a broader WSC-domain family issue, not a paralog the competition filter would catch |
| so/soft | 5/5 | **strong** hits in CneoH99/Ccin/Mcir/Scer (e-41 to e-58) | **excluded from controls.csv** — contradicts literature claim of Basidiomycota/yeast absence; needs its own investigation before reuse anywhere |

## Caveats

- `-` (no hit at e<=1 under `--very-sensitive`) is strong but not absolute evidence of
  true absence — same caveat as the hex-1/eIF-5A investigation this session started from.
- This check used the anchor's *own* best hit per proteome only; it does not run the
  pipeline's actual paralog-competition filter, so "outgroup contamination" here means
  "diamond finds *some* homology," not "the pipeline would call it present" — for hex-1
  and spa-18 specifically, whether they're usable controls depends on that filter
  actually firing correctly on a real run.
- `wsc`/`so` results were a surprise relative to general Pezizomycotina-vs-outgroup
  literature statements; not deeply investigated further (e.g. whether the outgroup hit
  is to a real ortholog, a shared domain, or a same-genome paralog) — flagged, not solved.
