# Methods: pangenome cluster-profiling pipeline (`pangenome.nf`)

Source of truth: the code at this checkout (origin/main). Every statement below was read
in code unless it carries the tag [UNVERIFIED]. File paths are quoted. Line references are
approximate. Defaults are from `"nextflow.config"` (params block, lines ~222-477) unless
stated otherwise.

---

## 1. Overview

The pipeline takes a set of annotated genomes in two groups: an ingroup and an outgroup.
It answers these questions:

1. Which gene families exist across the strains, and how often does each family occur
   (core, soft core, shell, cloud, singleton)?
2. Is a family absent because the gene is missing, or only because the annotation missed
   it? (optional genome-level TBLASTN rescue)
3. Which pairs of accessory families occur together across strains more often than
   chance, after a control for clade structure?
4. Is each co-occurring pair physically linked on the genome (within a gene window), or
   not linked ("trans")? Is a linked pair near a mobile-element marker gene ("captain")?
5. Which runs of accessory genes ("accessory islands") hold a linked co-occurring pair?
6. Which Pfam domains are enriched in island families against the accessory background?
7. Which unlinked co-occurring families form modules (Leiden communities), and do the
   module genes cluster in the genome?
8. At each island locus, which strains carry the locus, which lack it, and is a lack a
   real DNA deletion or a gene-model difference?
9. Is the pangenome open or closed (accumulation curve and Heaps' law fit)?
10. Does assembly quality confound accessory-gene counts?

Entry point: `"pangenome.nf"`. Subworkflow: `"workflows/pangenome_profile.nf"`
(`PANGENOME_PROFILE`). Processes: `"modules/pangenome/*.nf"`. Scripts:
`"bin/pangenome_*.py"`. Shared logic: `"lib/"`.

---

## 2. Inputs

### 2.1 Samplesheet (`--pangenome_samplesheet`, required)

CSV with header `GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup`
(`"pangenome.nf"` header; `"lib/config_parser.py"` `parse_config`).

| column | use in this pipeline |
|---|---|
| `GROUP` | Role. Rows with `GROUP == --pangenome_ingroup_label` (default `IN`) are the ingroup. Rows with `GROUP == --pangenome_outgroup_label` (default `OUT`) are the outgroup. `"pangenome.nf"` drops every other row before it builds the strain channel. |
| `Species` | Species label. Used by the island locus view (per-species counts, ranks "species" and "within") and by figure ordering. Not used by any statistical test upstream of the locus view. |
| `Strain` | Read by the parser. Not used by any pangenome step that was read. |
| `Protein` | Protein FASTA file name. `"pangenome.nf"` resolves it under `data_dir`, then `pep/`, then `proteins/`. |
| `DNA` | Genome FASTA file name. `"pangenome.nf"` resolves it under `data_dir`, then `dna/`, `genome/`, `scaffolds/`. |
| `GFF3` | GFF3 file name. Required. Must exist at `<gff3_dir>/<GFF3>`. The run stops at parse time if it is empty or missing. |
| `Short` | Strain ID. It prefixes every protein and contig ID (`Short|id`). It is the column key of the presence matrix. |
| `TaxonGroup` | Clade label. It is the stratification variable of the co-occurrence test and the clade count of pair classification. Empty cells are filled from Mash clades (section 3.5). |

### 2.2 Group roles

- Ingroup: the strain set for all frequency bins (`bin`), the co-occurrence test, and
  the clade stratification.
- Outgroup: the reference for gain/loss polarity (`direction_a`, `direction_a_freq`,
  Dollo tree columns) and for the outgroup bins (`bin_out`).
- Both groups: clustering, presence matrix, rescue, gene positions, Mash dereplication,
  pair classification, islands, locus view, accumulation curve.

`"lib/config_parser.py"` maps group aliases (`NEAR_IN`, `BROAD_OUT`) to canonical names.
`"pangenome.nf"` compares the raw `GROUP` string, without that alias map.

### 2.3 `data_dir` layout (`--pangenome_data_dir`, required)

- `pep/` (or flat, or `proteins/`): protein FASTA files.
- `dna/` (or flat, or `genome/`, `scaffolds/`): genome FASTA files.
- `gff3/`: GFF3 files. Override with `--pangenome_gff3_dir`.

Several scripts do not use the flat/subdir search. They use a fixed path:
`GENE_POSITIONS` reads `<data_dir>/pep/<Protein>`; `DEREPLICATE`, `ASSIGN_CLADES`,
`ASSEMBLY_QUALITY_QC` and `ISLAND_DNA_CHECK` read `<data_dir>/dna/<DNA>`
(`"bin/pangenome_build_gene_positions.py"`, `"bin/pangenome_dereplicate_strains.py"`,
`"bin/pangenome_assign_clades.py"`, `"bin/pangenome_assembly_quality_qc.py"`,
`"bin/pangenome_island_dna_check.py"`). See section 7.

### 2.4 Optional inputs

- `--pangenome_captain_hmm` (an HMM file), or `--pangenome_captain_hmm_name` plus
  `--pangenome_pfam_hmm` (fetch that model by name from a Pfam database).
- `--pangenome_island_pfam_hmm` (Pfam-A.hmm). This gates Pfam scan, domain enrichment,
  module domains, the island locus view, clinker panels and the island synteny page.
- `--pangenome_marker_names` / `--pangenome_marker_hmm_paths` (named marker HMMs).
- `--pangenome_pfam2go` (local pfam2go file).
- `--pangenome_species_tree` (rooted Newick, tips = `Short` IDs of all matrix strains).

---

## 3. Steps

### 3.1 Samplesheet parse and input checks

- Purpose: build the per-strain channel and fail early on bad inputs.
- Algorithm:
  1. Read the CSV. Keep rows whose trimmed `GROUP` is the ingroup or outgroup label.
  2. Resolve the protein and DNA files (section 2.1).
  3. Require a GFF3 value and an existing GFF3 file per row.
  4. Validate locus-view, clinker and backend parameters (integer and range checks).
- Tool: Nextflow (`"pangenome.nf"`).
- Output: channel `[meta(id=Short, group, taxon), protein_fa, dna_fa]`; a second channel
  with ingroup DNA only.

### 3.2 Proteome and genome prefixing, concatenation

- Purpose: make every sequence ID carry its strain.
- Algorithm: rewrite the first token of each FASTA header to `<Short><id_sep><id>`. Keep
  the rest of the header. Concatenate all prefixed proteomes into one file.
- Processes: `PREFIX_PROTEOME`, `PREFIX_GENOME`, `CONCAT_PROTEOMES`
  (`"modules/pangenome/prefix_and_cluster.nf"`); script `"bin/pangenome_prefix_fasta.py"`.
- Parameter: `pangenome_id_sep`, default `|`.
- Outputs: `<Short>.prefixed.pep.fa`, `<Short>.prefixed.dna.fa`, `all_strains.pep.fa`.
- Note: the pipeline does not collapse isoforms. The docstring of
  `"bin/pangenome_cluster_backend.py"` says the input should be isoform-collapsed.
  No isoform step is wired into `"workflows/pangenome_profile.nf"`.

### 3.3 Tier-1 clustering

- Purpose: define the gene family. Family ID = the cluster representative sequence ID,
  verbatim.
- Algorithm (mmseqs backend, default): `mmseqs easy-cluster all_strains.pep.fa tier1 tmp
  --min-seq-id 0.9 -c 0.8 --cov-mode 0 --cluster-reassign --threads N`
  (`"bin/pangenome_cluster_backend.py"` `run_mmseqs_cluster`).
- Algorithm (diamond backend): `diamond cluster -d all_strains.pep.fa --approx-id 90
  --member-cover 80`. Then `"bin/verify_diamond_cluster_ids.py"` checks every ID. The
  representatives (column 1, unique) are extracted with awk to build
  `tier1_rep_seq.fasta`.
- Process: `CLUSTER_TIER1`. Label `high_cpu`.
- Parameters: `pangenome_cluster_backend` = `mmseqs`; `pangenome_tier1_min_id` = `0.9`;
  `pangenome_tier1_cov` = `0.8`. For diamond the fractions are multiplied by 100.
- Outputs: `tier1_cluster.tsv` (rep, member), `tier1_rep_seq.fasta`.
- The mmseqs cluster mode and sensitivity are mmseqs defaults; the code does not set
  them [UNVERIFIED which defaults apply].
- A tier-2 "superfamily" mode exists in the script (`mmseqs-tier2`, `diamond-tier2`). It is
  not wired into the workflow.

### 3.4 Presence matrix

- Purpose: family x strain table of calls.
- Algorithm:
  1. Read `tier1_cluster.tsv` into {member: rep}. Build {rep: [members]}. Singleton
     clusters are kept as families.
  2. Split each member ID on the first `id_sep` to get the strain.
  3. Columns = ingroup + outgroup strains from the samplesheet.
  4. Cell = `present` when the strain has >= 1 member; otherwise `absent`.
  5. Copy number = member count per (family, strain), in a sidecar file.
  6. Hard error if every member has an unknown strain prefix.
- Process: `PRESENCE_MATRIX`; script `"bin/pangenome_build_presence_matrix.py"`; data
  class `"lib/pangenome_matrix.py"` `PresenceMatrix`.
- Outputs: `presence_matrix.tsv` (`family`, then one column per strain, values
  `present|genome_only|absent`); `presence_matrix.copy_number.tsv`
  (`family, strain, copies`).
- Rule used by every later step: `is_present` is true for `present` and `genome_only`
  (`"lib/pangenome_matrix.py"` line ~113). No step before the locus view separates the two.

### 3.5 Gene positions

- Purpose: protein -> genomic position per strain.
- Algorithm (`"bin/pangenome_build_gene_positions.py"`):
  1. For each ingroup/outgroup strain, parse GFF3 `CDS` rows.
  2. Key each CDS on `protein_id=`. If absent, key on `Parent=` (split on comma).
  3. Collapse all CDS rows of one ID to (contig, min start, max end).
  4. Keep only IDs found in the strain's protein FASTA headers (`<data_dir>/pep/<Protein>`).
  5. If < 50% of the strain's proteins match, try a third join: FASTA `GN=` to GFF3
     `locus_tag=` (1:1 names only). Use it if it matches more.
  6. If the match fraction is still < 0.50, stop the run. If unmatched > 2%, warn.
- Process: `GENE_POSITIONS` (runs right after `PRESENCE_MATRIX`).
- Output: `gene_positions.tsv.zst` (`Short, protein_id, contig, start, end`).

### 3.6 Optional rescue pass (TBLASTN)

Gate: `pangenome_rescue_enable` = `true`. With `false`, the unrescued matrix is used
everywhere and empty stub files replace the rescue outputs.

- Purpose: upgrade an `absent` cell to `genome_only` when the strain's genome holds the
  gene but its annotation does not.
- Algorithm:
  1. `EXTRACT_ABSENT_QUERIES` (`"bin/pangenome_extract_absent_family_queries.py"`): for
     each strain, write the representative sequences of every family whose cell is
     `absent` in that strain (`<Short>.absent.fa`). Both ingroup and outgroup strains get
     queries.
  2. `MAKE_STRAIN_GENOME_DB`: `makeblastdb -dbtype nucl` on that strain's prefixed genome.
  3. `TBLASTN_PER_STRAIN`: `tblastn -outfmt "6 std qcovs" -evalue 1e-10
     -max_target_seqs 200`, query = the strain's absent reps, db = the strain's own genome.
  4. `RESCUE_PASS` (`"bin/pangenome_rescue_pass.py"`):
     a. Threshold: keep an HSP row when `pident >= 90.0` and `qcovs >= 80.0`. Note:
        `pident` is per HSP; `qcovs` is the per-subject query coverage reported by BLAST.
     b. Structural filter (`pangenome_rescue_structural_filter` = `true`), in this order:
        - overlap: reject the HSP if its subject span overlaps a GFF3 gene (from
          `gene_positions`) in that strain that belongs to a different family;
        - short rep: reject if the family representative is < `150` aa;
        - hotspot: bin each surviving HSP by `floor(sstart / 2000)` per (strain, contig);
          reject the HSP if its bin holds >= `3` distinct families. The bins are fixed,
          not sliding windows.
     c. A (family, strain) pair with >= 1 surviving HSP is a rescue hit.
     d. Apply: set the cell to `genome_only` only if it is `absent`. Never downgrade
        `present`.
     e. Stop the run if input files were given but zero rows were parsed, or if every
        hit had an unknown strain/family.
  5. `EXTRACT_RESCUE_POSITIONS` (`"bin/pangenome_extract_rescue_positions.py"`): for each
     `genome_only` cell, keep the single highest-bitscore HSP that passes the
     identity/coverage thresholds. Position = (contig, min(sstart, send)). This step does
     not apply the structural filter (see section 7).
- Parameters: `pangenome_rescue_evalue` 1e-10; `pangenome_rescue_max_target_seqs` 200;
  `pangenome_rescue_min_pident` 90.0; `pangenome_rescue_min_qcov` 80.0;
  `pangenome_rescue_min_rep_length` 150; `pangenome_rescue_hotspot_window` 2000;
  `pangenome_rescue_hotspot_min_families` 3.
- Outputs: `presence_matrix.rescued.tsv`; `rescue_funnel.tsv` (`rows_parsed`,
  `rows_passed_threshold`, `rows_rejected_overlap`, `rows_rejected_short_rep`,
  `rows_rejected_hotspot`, `applied`, `skipped`); `rescue_positions.tsv`
  (`Short, family, contig, start`); per-strain `<Short>.tblastn.tsv.zst`.

### 3.7 Mash dereplication

Gate: `pangenome_dereplicate` = `true`.

- Purpose: pick one representative per group of near-identical genomes. Representatives
  form the frequency and co-occurrence denominators.
- Algorithm:
  1. `MASH_SKETCH_ALL`: `mash sketch` and `mash dist -t` over all ingroup + outgroup
     genomes. Sketch size and k-mer size are mash defaults; the code passes neither
     [UNVERIFIED: mash defaults k=21, s=1000; `"docs/pangenome-assumptions.md"` sec. 3
     states s=1000].
  2. `DEREPLICATE` (`"bin/pangenome_dereplicate_strains.py"`): build a graph with an edge
     where distance `< 0.001` (strict). Groups = connected components (union-find, i.e.
     single linkage).
  3. Assembly stats per strain (n_contigs, N50, total length) from `<data_dir>/dna/<DNA>`.
  4. Representative = the member with the highest N50; tie -> lexically last `Short`.
- Parameter: `pangenome_mash_threshold` = `0.001`.
- Output: `strain_inventory.tsv` (`Short, n_contigs, n50, total_length, dedup_group,
  is_representative`).
- A group can mix ingroup and outgroup strains. The code does not prevent this.

### 3.8 Clade assignment and TaxonGroup fill

Gate: `pangenome_assign_clades` = `true`.

- Purpose: provide a stratification label when `TaxonGroup` is empty.
- Algorithm:
  1. `MASH_SKETCH_INGROUP`: mash sketch/dist on ingroup genomes only.
  2. `ASSIGN_CLADES` (`"bin/pangenome_assign_clades.py"`): symmetrize the distance
     matrix; classical PCoA (double-centering, eigendecomposition, keep up to 4 positive
     axes); k-means (`scipy.cluster.vq.kmeans2`, `minit="++"`, seed 0) for k in 2..10;
     choose the k with the highest mean silhouette. With `pangenome_clades_k_fixed` set,
     skip the sweep.
  3. `FILL_TAXON_GROUP` (`"bin/pangenome_fill_taxon_group.py"`): write `clade_<k>` into
     every empty `TaxonGroup` cell. Never overwrite a non-empty cell.
- Parameters: `pangenome_clades_k_range` `'2,10'`; `pangenome_clades_k_fixed` `null`;
  `pangenome_clades_n_pcoa` `4`.
- Outputs: `clade_assignments.tsv` (`Short, mash_clade, pcoa1..`),
  `samplesheet.with_clades.csv` (the "effective samplesheet").
- Outgroup strains get no Mash clade. Their `TaxonGroup` stays as given.
- All non-rep ingroup strains also get a label; only reps enter the tests.

### 3.9 Frequency bins

- Purpose: class each family by how many strains of each group carry it.
- Process `FREQUENCY_BINS`; script `"bin/pangenome_frequency_bins.py"`.
- Strain sets (per group G):
  - `reps(G)`: G strains in the matrix, and `is_representative == 1` when dereplication is
    on. All G strains when it is off.
  - `nonreps(G)`: the other G strains in the matrix.
- Rules:
  - `count_G(f)` = number of `reps(G)` with `is_present(f)`.
  - `freq_G(f) = count_G(f) / |reps(G)|`.
  - If `count_G >= 1`, apply `assign_bin` in this order:
    `count_G <= 1` -> `singleton`; `freq >= 0.95` -> `core`; `freq >= 0.90` ->
    `soft_core`; `freq >= 0.15` -> `shell`; else `cloud`.
  - If `count_G == 0` and present in any `nonreps(G)` strain -> `nonrep_only`.
  - Else if present in any strain (rep or non-rep) of the other group -> `outgroup_only`
    (ingroup side) or `ingroup_only` (outgroup side).
  - Else -> `absent`.
- Outgroup binning: computed only when `|reps(OUT)| >= pangenome_outgroup_min_bin_strains`
  (default 3). Otherwise the `*_out` columns are `-`.
- Parameters: `pangenome_core_cutoff` 0.95; `pangenome_softcore_cutoff` 0.90;
  `pangenome_shell_cutoff` 0.15; `pangenome_outgroup_min_bin_strains` 3.
- Output: `frequency_table.tsv` (`family, frequency, strain_count, bin, frequency_out,
  strain_count_out, bin_out`).
- Note: `singleton` is decided by count (`<= 1`) before any frequency cutoff. With a
  very small ingroup (e.g. 1 rep) every present family is `singleton`.

### 3.10 Co-occurrence test

- Process `COOCCURRENCE`; script `"bin/pangenome_cooccurrence.py"`. Label `high_cpu`.
- Strain set S = ingroup strains, restricted to reps when dereplication is on.
- Eligible families: `bin` in {`shell`, `cloud`} AND carried by >= `pangenome_min_strain_count`
  (5) strains of S. `soft_core`, `core`, `singleton` and the zero-count classes are not tested.
- Hypotheses: all unordered pairs of eligible families. m = E(E-1)/2.
- Per pair: 2x2 table over S (both, a only, b only, neither) from bitsets.
- Screen (speed only): a continuity-corrected one-sided normal approximation
  `z = (both - E[both] - 0.5) / sqrt(var)`, `E[both] = r1 c1 / n`,
  `var = r1 c1 (n-r1)(n-c1) / (n^2 (n-1))`, `p_screen = 0.5 erfc(z / sqrt 2)`; 1.0 for
  `z <= 0` or a zero margin. Pairs with `p_screen > screen_alpha` are dropped without an
  exact test. `screen_alpha` is clamped to `>= fdr_alpha`. Default `screen_alpha` = 0.2.
- Exact test: one-sided Fisher exact test, `alternative="greater"`
  (`scipy.stats.fisher_exact`) on screen survivors. Only positive association
  (co-occurrence) is tested. Mutual exclusion is not tested.
- Multiple testing: Benjamini-Hochberg (`scipy.stats.false_discovery_control`,
  method `bh`) over m hypotheses. Screened-out pairs enter as p = 1.0. A pair is kept
  when `fdr_q < pangenome_fdr_alpha` (0.05, strict `<`).
- Within-clade stratified exact test (FDR survivors only):
  - Clade of a strain = its `TaxonGroup` in the effective samplesheet; empty -> one
    shared `"unknown"` stratum.
  - Statistic: `both`, the overlap count.
  - Null: within each clade c, `X_c ~ Hypergeom(n_c, |a_c|, |b_c|)`, independent.
    Clades with `|a_c| = 0` or `|b_c| = 0` are skipped.
  - `permutation_p = P(sum_c X_c >= both_observed)`, by convolution of the per-clade PMFs.
    It is exact. No simulation runs. The name "permutation" is historical.
  - No multiple-testing correction is applied to `permutation_p`.
  - A warning is printed when > 50% of S has an empty `TaxonGroup`.
- Polarity (family A of the pair only):
  - Outgroup set = all outgroup strains that are matrix columns (reps and non-reps).
  - `direction_a` (strict): outgroup count = all -> `loss`; = 0 -> `gain`; otherwise or
    no outgroup -> `ambiguous`.
  - `asymmetry_a` = outgroup count / outgroup total.
  - `direction_a_freq`: `asymmetry_a >= 0.9` -> `loss`; `<= 0.1` -> `gain`; else
    `ambiguous`. Requires `0 <= gain_max < loss_min <= 1`.
  - Optional Dollo columns with `--pangenome_species_tree` (`"lib/ancestral_states.py"`):
    the tree must hold every matrix strain; extra tips are pruned. Single gain at the
    MRCA of all carriers (all matrix strains, reps and non-reps); each maximal all-absent
    clade below it = one loss; a multifurcation counts its absent children as one loss.
    `direction_a_tree`: `gain` if the gain clade lies inside the ingroup set, `loss` if it
    spans outside and some ingroup strain lacks the family, else `conserved`,
    `outgroup_only`, `absent`. The ingroup set here is S (reps only).
- Other columns: `jaccard = both / (both + a_only + b_only)` over S;
  `clade_composition` = counts of family A carriers in S per clade label.
- Parameters: `pangenome_min_strain_count` 5; `pangenome_fdr_alpha` 0.05;
  `pangenome_screen_alpha` 0.2; `pangenome_polarity_loss_min_frac` 0.9;
  `pangenome_polarity_gain_max_frac` 0.1; `pangenome_species_tree` null.
- Output: `cooccurring_pairs.tsv.zst` (`family_a, family_b, jaccard, fisher_p, fdr_q,
  permutation_p, direction_a, clade_composition, asymmetry_a, direction_a_freq`, plus
  `direction_a_tree, gain_node, n_loss_events, loss_clades, asr_method` with a tree).
  Only FDR-significant pairs are written. `fisher_p`/`fdr_q` are written with 3
  significant digits (`%.2e`); `permutation_p` with 4 decimals (`%.4f`).

### 3.11 Family positions

- Purpose: per-strain gene-order ranks for each family.
- Algorithm (`"bin/pangenome_build_family_positions.py"`):
  1. Per strain, merge annotated genes (`gene_positions`) and rescue positions
     (`rescue_positions`, when rescue is on).
  2. Sort by (contig, start). Rank = the running index over the whole sorted list.
     The index does not restart per contig. Ranks are compared only within one contig.
  3. Map annotated genes to families through `tier1_cluster.tsv`. Unclustered genes keep
     their rank slot but are not written.
  4. A multi-copy family gets one row per copy.
- Process: `FAMILY_POSITIONS`. Output: `family_positions.tsv.zst`
  (`Short, family, contig, rank`).

### 3.12 Captain marker search (optional)

- Run only when a captain HMM is given.
- `HMMFETCH_CAPTAIN` (`hmmfetch <pfam> <name>`) when only a name is given.
- `CAPTAIN_HMMSEARCH`: `hmmsearch --tblout -E 1e-3` against `all_strains.pep.fa`.
- Each hit protein is mapped to its tier-1 family per strain.
- Parameter: `pangenome_captain_hmm_evalue` 1e-3. Output: `captain_vs_study.tblout`.

### 3.13 Pair classification

- Process `PAIR_CLASSIFICATION`; script `"bin/pangenome_pair_classification.py"`;
  `"lib/pangenome_synteny.py"` `linkage_fraction`.
- Input: every FDR-significant pair from section 3.10.
- Strain set: every strain in `family_positions` (ingroup and outgroup, reps and
  non-reps, annotated and rescued positions).
- Rules, in order:
  1. `n_co` = strains with >= 1 position for both families. If `n_co < 5` ->
     `insufficient_data` (linkage reported as 0.0).
  2. `linkage_fraction` = fraction of the `n_co` strains where some copy of A and some
     copy of B are on the same contig with `|rank_A - rank_B| <= k` (k = 10).
  3. `linkage >= 0.5`: if a captain family lies within k ranks of A or of B on the same
     contig, in ANY strain that has positions for A or B (not only co-carrying strains),
     -> `starship_explained`; else -> `unexplained_physical`.
  4. `linkage <= 0.05`: if `permutation_p < 0.05` AND the number of distinct labels in
     `clade_composition` >= 2 -> `trans`; else `trans_unconfirmed`.
  5. Otherwise -> `ambiguous_linkage`.
- `clade_composition` comes from family A carriers only. An `"unknown"` label counts as a
  clade.
- Parameters: `pangenome_pair_class_k` 10; `pangenome_pair_class_physical_threshold` 0.5;
  `pangenome_pair_class_trans_threshold` 0.05; `pangenome_pair_class_min_co_carrying` 5;
  `pangenome_pair_class_perm_alpha` 0.05; `pangenome_pair_class_min_clades` 2.
- Output: `pair_classification.tsv.zst` (`family_a, family_b, classification,
  linkage_fraction, jaccard, fisher_p, fdr_q, permutation_p, direction_a,
  clade_composition, asymmetry_a, direction_a_freq`). The Dollo tree columns are not
  passed through.

### 3.14 Accessory islands

- Process `BUILD_ISLANDS`; script `"bin/pangenome_build_islands.py"`. Runs on every run.
- Definitions:
  - Core family: ingroup `bin` in {`core`, `soft_core`}. Every other `bin` value is
    non-core. A family missing from the table defaults to core.
  - Accessory island (per strain): a maximal run of consecutive non-core families in
    rank order, on one contig. All strains in `family_positions` are scanned, including
    outgroup strains, using the ingroup `bin`.
  - "Significant" island: size >= 2 genes AND it holds both members of >= 1 pair
    classified `starship_explained`, `unexplained_physical` or `ambiguous_linkage`.
  - Size filter: `pangenome_island_min_size` (2).
- Deduplication across strains: key = (set of member families, set of supporting
  classifications). Gene order and copy number are not in the key. `n_strains` = number
  of per-strain islands with that key. `example_strain` = the first strain that produced
  the key, in `family_positions` file order. `island_size` and `member_families` come from
  that first island.
- Markers: for each named marker, `has_<name>` = `Y` if any member family has a hit
  (`MARKER_HMMSEARCH`, `hmmsearch --tblout -E 1e-5`).
- Output: `significant_islands.tsv` (`n_strains, example_strain, island_size,
  member_families, n_supporting_pairs, classifications, has_<marker>...`), sorted by
  `n_strains` descending.
- No test statistic is computed for an island. "Significant" means only the
  pair-containment rule above.

### 3.15 Leiden trans modules

- Process `LEIDEN_MODULES`; script `"bin/pangenome_detect_trans_modules.py"`.
- Graph: nodes = families in `trans` pairs; edges = `trans` pairs; weight =
  `min(300, -log10(fdr_q))` (300 when `fdr_q <= 0`).
- Algorithm: `leidenalg.find_partition`, `RBConfigurationVertexPartition`, resolution 1.0,
  seed 0. Modules are renumbered by size (0 = largest). Zero `trans` edges -> empty
  header-only outputs.
- Outputs: `family_modules.tsv` (`family, module_id, module_size`), `module_summary.tsv`.
- `MODULE_DOMAINS` (Pfam branch only): per module with >= 2 families, count families
  with >= 1 domain row in `pfam.domtblout` and list the top 15 domain names. It applies no
  extra i-Evalue filter beyond the scan's `--domE`.
- `MODULE_NEIGHBORHOOD` (`"lib/pangenome_neighborhood.py"`): see section 3.16.

### 3.16 Module genomic neighbourhood

- Purpose: test whether module genes sit close together within each strain's assembly.
- Genes: annotated genes only (`gene_positions`), all strains (ingroup and outgroup).
  Rank = 0-based order by start within (strain, contig).
- For each module with >= 2 families and each strain with >= 2 member genes, count
  pairs of member genes from different families:
  `total`; `cross_contig` (excluded); `adjacent` (same contig, rank gap < 11, excluded);
  `informative` (same contig, rank gap >= 11); `colocalized` (informative and starts within
  100 kb).
- Statistics: `obs_frac = sum(colocalized) / sum(informative)`;
  `obs_frac_total = sum(colocalized) / sum(total)`. Sums are over strains.
- Null: in each scored strain, draw the same number of genes, without replacement, from
  the strain's pool (accessory genes = families whose ingroup `bin` is not core/soft_core,
  plus the module genes). Recompute. 200 draws. Seed 0.
- `p_empirical = (1 + #{null >= obs}) / (1 + n_valid_perm)`. `effect_ratio = obs / null mean`.
- No multiple-testing correction across modules.
- Parameters: `pangenome_neighborhood_max_kb` 100; `pangenome_neighborhood_min_gene_gap`
  11; `pangenome_neighborhood_n_perm` 200; `pangenome_neighborhood_seed` 0;
  `pangenome_neighborhood_null_pool` `accessory`; `pangenome_module_min_size` 2.
- Output: `module_neighborhood.tsv`.

### 3.17 Pfam scan of background representatives and domain enrichment

Gate: `pangenome_island_pfam_hmm` set.

- `SELECT_BACKGROUND_REPS` (`"bin/pangenome_select_background_reps.py"`): background =
  families with ingroup `bin` in {`shell`, `cloud`}. Write their tier-1 representative
  sequences.
- `HMMPRESS_PFAM` once; `FAMILY_PFAM_SCAN`: `hmmscan --domtblout --domE 1e-3` per chunk of
  3000 sequences; `MERGE_PFAM_DOMTBLOUT` concatenates to `pfam.domtblout`.
- `DOMAIN_ENRICHMENT` (`"bin/pangenome_domain_enrichment.py"`):
  1. Domain hits: keep domtblout rows with i-Evalue (column 13) <= 1e-3. Family -> set of
     domain names.
  2. Island set I = union of `member_families` over ALL rows of `significant_islands.tsv`,
     intersected with the background B. One pooled test, not one test per island.
  3. Tested domains = domains seen in >= 1 family of I.
  4. Per domain d, 2x2 over families of B: [[in I with d, in I without d], [not in I with d,
     not in I without d]]. One-sided Fisher exact test, `alternative="greater"`.
  5. BH-FDR over tested domains. The log counts `fdr_q < 0.05` as significant.
- Optional `PFAM2GO` (`"bin/pangenome_pfam2go.py"`): add GO IDs per Pfam accession.
- Outputs: `pfam.domtblout`, `island_pfam_enrichment.tsv` (`domain, pfam_accession,
  pfam_url, n_with_domain_in_islands, n_with_domain_in_background, n_island_families,
  n_background_families, fisher_p, fdr_q`).
- Parameters: `pangenome_pfam_domain_evalue` 1e-3; `pangenome_pfam_chunk_size` 3000;
  `pangenome_pfam2go` null.

### 3.18 Report tables

- Process `REPORT_TABLES`; script `"bin/pangenome_report_tables.py"`. Runs on every run.
- `islands_with_domains.tsv`: island rows plus `pfam_domains` (union of member domains)
  and a locus span in the example strain (`locus_id` = `strain:contig:start-end`; `-` when
  no member resolves or members fall on > 1 contig).
- `island_size_distribution.tsv`, `classification_counts.tsv`, `marker_summary.tsv`.
- `per_strain_summary.tsv`: per matrix strain, families present (present or genome_only),
  tallied by the strain's own group's class (`bin_out` for outgroup strains when the
  outgroup is binned). Outlier flag: modified z-score `0.6745 (x - median) / MAD` of the
  `singleton` count; median and MAD from the group's reps; `|z| > 3.5` -> outlier;
  `-` when n < 3 or MAD = 0.
- `group_class_overlap.tsv`: 7 x 7 counts of (ingroup class, outgroup class); empty when
  the outgroup is not binned.

### 3.19 Diagnostics and assembly-quality QC

- `ASSEMBLY_QUALITY_QC` (`"bin/pangenome_assembly_quality_qc.py"`), every run:
  - Strains: all ingroup strains with a DNA file (reps and non-reps).
  - Per strain: n_contigs, N50, total length; `accessory_present` = families present whose
    ingroup `bin` is not `core` (so `soft_core` counts as accessory here);
    `core_missing`; private-protein terminus fraction: proteins of `singleton` families
    with exactly one carrier among all matrix strains, located within 1000 bp of a contig
    end, vs the same fraction for all proteins.
  - Spearman rho (average ranks) and partial Spearman controlling total length, for
    {accessory_present, core_missing} x {n_contigs, N50}. No p-values are computed.
  - WARNING when |rho| (raw or partial) for accessory_present > 0.3.
  - Outputs: `assembly_quality_vs_content.tsv`, `assembly_quality_correlations.tsv`,
    `assembly_quality_report.md`. It never removes a strain.
- `DIAGNOSTICS` (`"bin/pangenome_diagnostics.py"`):
  - `rescue_redundancy`: triggered when `rows_rejected_overlap / rows_passed_threshold > 0.5`.
  - `assembly_quality_confound`: triggered when max |rho| for accessory_present > 0.3.
  - Clade validity and gain/loss skew: always `not_computed`.
  - `pangenome_strict` = false: advisory only. True: a triggered diagnostic fails the run.
  - Outputs: `diagnostics.tsv`, banner `.md` and `.html`.

### 3.20 Openness and accumulation curve

- In `REPORT_RENDER` (`"bin/pangenome_report_render.py"`).
- Strains: all matrix columns (ingroup and outgroup, reps and non-reps).
- Accumulation: 20 random strain orders (`random.Random(0)`); for each order, pan size =
  families present in >= 1 strain so far; core size = families present in all strains so
  far (100% rule, not the 95% cutoff). Mean and SD per N.
- Heaps fit: least squares of `log(pan_mean)` on `log(N)`, N = 1..n. `gamma` = slope,
  `kappa = exp(intercept)`, R^2 in log space. `is_open = gamma < 1`.
- Core decay: `core(N) = core_inf + A exp(-N / tau)` (`scipy.optimize.curve_fit`).
- Output: `report/pangenome_openness.tsv` (`kappa, gamma, r_squared, is_open, core_inf, tau`).
- See section 6 on the openness rule.

### 3.21 Report render

- `REPORT_RENDER` writes `report/report.md` plus PNG and PDF figures: frequency histogram,
  per-group class composition, group overlap heatmap / shared split / UpSet, presence-absence
  raster (all matrix strains), accumulation curve, classification counts, island sizes,
  per-genome class composition, top enriched domains, module neighbourhood table.
- "Top islands (by size)": islands with `n_strains >= 2` (`pangenome_top_islands_min_strains`).
- Enriched-domain highlight uses `fdr_q < 0.05` (`--fdr_threshold`, not a Nextflow param).
- The diagnostics banner is prepended.

### 3.22 Island locus view

Gate: `pangenome_island_pfam_hmm` set (runs inside that branch).

- Script `"bin/pangenome_island_loci.py"`; logic `"lib/island_locus.py"`.
- Strain set: all matrix columns. Species = samplesheet `Species` (original samplesheet).
  N50 = `assembly_quality_vs_content.tsv` (ingroup only).
- Loci (`group_loci`): islands with a resolved `locus_id`, sorted by size descending. An
  island joins the locus of a strictly larger island when >= 0.5 of its families are in
  that island (most shared families wins; ties to the earlier one).
- Candidates: loci whose summed variant `n_strains` (capped at the strain count) is >= 2;
  top 200 by that proxy (or by size when rank = `size`).
- Carrier of a locus (Ruling R3): one contig holds every member family within a rank span
  `<= (#families - 1) + k`.
- Exemplar: among carriers, tier `full` (>= 5 genes on both sides of the block on its
  contig), else `short_flanks` (>= 3), else `contig_end` (most flank genes). Within a tier:
  highest N50; strains without N50 (outgroup) rank after, by gene count of the contig;
  then name.
- Columns: up to F flank genes left, every gene between the block ends, up to F flank genes
  right, in the exemplar's order (F = 5, or 3 for `short_flanks`).
- Cell states per strain and column (k = 10): `in place` (a copy with a copy of a different
  column family within k ranks on the same contig); `elsewhere` (a copy without such a
  neighbour, or present in the matrix with no position); `absent`; `rescue (in place/
  elsewhere)` when the matrix call is `genome_only`; `contig break` (absent/elsewhere,
  outside the in-place column range, and the nearest in-place copy is within k ranks of a
  contig end).
- Flanks intact: some in-place left-flank copy and some in-place right-flank copy on one
  contig within `n_locus + 2k` ranks.
- Row class: not intact -> `uninformative`; >= 80% of locus columns absent -> `empty`;
  all locus columns in place -> `full`; else `partial` (also when columns are only
  `elsewhere`).
- DNA presence check (`pangenome_locus_dna_check` = true): see section 3.23.
- Ranks (per candidate, after DNA calls): carriers = `full`, `model_difference`, and
  `partial` with no missing column; losses = `empty` and `partial` with a missing column
  (missing = `DNA absent` with the check, `absent` without). `uninformative` excluded.
  - `whole_annot` = min(empty_n, full) if empty_n >= 10 and full >= 2, else -1.
  - `whole_dna` = same with full + model_difference.
  - `presence` = min(losses, carriers) if losses >= 10 and carriers >= 2, else -1.
  - `species` = max_s f_s - min_s f_s over species with n_s >= 10 (hard-coded), if >= 2
    such species and the difference >= 0.95; else -1. f_s = losses_s / n_s.
  - `within` = max over species with n_s >= 20 and 0.05 <= f_s <= 0.95 of
    min(car_s, loss_s); else -1.
  - empty_n = DNA-confirmed empty count with the check on.
- Drawn set: top 20 (score >= 0) under whole_annot, whole_dna, species, within; then fill in
  presence order up to 100 loci. Output order = presence order; keys `L001`, ...
- Outputs: `island_loci.json`, `island_regions.tsv`.
- These ranks are descriptive scores. No p-value is computed.

### 3.23 DNA presence check

- `ISLAND_DNA_TARGETS` (pass 1 of the same script) writes work lists for every candidate
  locus, batched 50 loci per file.
- Checked strains: flank-intact strains with >= 1 locus column not in place.
- Query: exemplar DNA from the first to the last locus gene; rescue-only exemplar columns
  use the TBLASTN HSP span.
- Target: the strain's DNA from its innermost in-place left-flank gene start to its
  innermost in-place right-flank gene end (flank genes included).
- `ISLAND_DNA_CHECK` (`"bin/pangenome_island_dna_check.py"`): `blastn -task megablast
  -query <exemplar locus> -subject <target> -outfmt "6 qstart qend pident"`. No e-value
  option is passed [UNVERIFIED: blastn default e-value 10]. Per exemplar gene: merge HSPs
  with pident >= 90; coverage = covered bp / gene length; `present` if coverage >= 80%,
  else `absent`. Target < 50 bp -> `absent`. Missing genome -> `unchecked`.
- Pass 2 (`ISLAND_LOCI`): an `absent`/`elsewhere`/`rescue elsewhere` locus cell becomes
  `absent, DNA present` or `absent, DNA absent`. Row class: >= 80% DNA absent -> `empty`;
  all in place -> `full`; all in place or DNA present (>= 1 DNA present) ->
  `model_difference`; else `partial`.
- Parameters: `pangenome_locus_dna_min_id` 90; `pangenome_locus_dna_min_cov` 80;
  `pangenome_locus_dna_batch` 50; `pangenome_locus_candidates` 200.

### 3.24 Clinker synteny panels and GenBank slices

Gate: `pangenome_clinker` = true (inside the Pfam branch).

- Strain choice per drawn locus (`select_clinker_strains`): the exemplar; then per row
  class (full, partial, empty, model_difference) and per species, the flank-intact strain
  with the best quality key; then more full/partial strains by quality; cap 12.
  `uninformative` strains are never chosen (the exemplar is always kept).
- Region: on the flank-pair contig, from the lowest to the highest in-place column copy
  within `n_locus + 2k` ranks, plus 5 genes each side, clipped at contig ends.
- `ISLAND_GBK_SLICE` (`"bin/pangenome_island_gbk_slice.py"`, `"lib/genbank_slice.py"`):
  rebuild ranks; check anchor families match (exit non-zero if not); split a region into
  blocks at gene-free gaps > 20000 bp (`pangenome_clinker_max_gap`; 0 = no split); write
  `<locus>/<strain>.gbk` (gene + CDS, translation, `family=` note), `groups.csv`
  (locus_tag, family), `clinker_order.txt`, and `island_slices.tsv`.
- `ISLAND_CLINKER` (`"bin/pangenome_island_clinker.py"`): `clinker <gbks> -gf groups.csv
  -j <cpus> -p <out>.html -ufo`; 50 loci per task; strip embedded sequences when
  `pangenome_clinker_slim` = true (`"lib/clinker_html.py"`). A failed locus is skipped.
  How clinker scores links between genes is clinker's own method [UNVERIFIED].
- `ISLAND_SYNTENY` (`"bin/pangenome_island_synteny.py"`, `"lib/island_synteny.py"`):
  `island_synteny.html`; islands with a resolved locus and `n_strains >= 2`, top 50 by
  size; strains collapsed to haplotype rows; embeds the loci JSON and links the clinker pages.

---

## 4. Denominators

| statistic | strain set | source |
|---|---|---|
| `frequency`, `strain_count`, `bin` | ingroup reps (all ingroup if dereplication off) | `"bin/pangenome_frequency_bins.py"` |
| `nonrep_only` (ingroup side) | ingroup non-reps | same |
| `outgroup_only` (ingroup side) | all outgroup strains in matrix (reps and non-reps) | same |
| `frequency_out`, `strain_count_out`, `bin_out` | outgroup reps; only if >= 3 reps | same |
| co-occurrence eligibility (`min_strain_count`), 2x2 table, Fisher p, BH m, `jaccard`, `clade_composition`, stratified p | ingroup reps | `"bin/pangenome_cooccurrence.py"` |
| `direction_a`, `asymmetry_a`, `direction_a_freq` | ALL outgroup strains in matrix (reps and non-reps) | same, lines ~528-543 |
| Dollo carriers | all matrix strains (ingroup + outgroup, reps + non-reps) | same, `dollo_for` |
| Dollo `direction_a_tree` ingroup set | ingroup reps | same |
| `min_co_carrying`, `linkage_fraction` | all strains with positions (ingroup + outgroup, reps + non-reps) | `"bin/pangenome_pair_classification.py"` |
| captain evidence | any strain with positions for A or B | same |
| island construction, island `n_strains` | all strains in `family_positions`; core status from ingroup `bin` | `"bin/pangenome_build_islands.py"` |
| Pfam enrichment | families (unit), background = ingroup `bin` shell+cloud | `"bin/pangenome_domain_enrichment.py"` |
| module neighbourhood | all strains in `gene_positions`; accessory pool from ingroup `bin` | `"lib/pangenome_neighborhood.py"` |
| assembly QC | all ingroup strains with DNA (reps + non-reps) | `"bin/pangenome_assembly_quality_qc.py"` |
| per-strain outlier median/MAD | reps of the strain's own group | `"bin/pangenome_report_tables.py"` |
| accumulation curve, Heaps fit, presence raster | all matrix strains (ingroup + outgroup, reps + non-reps) | `"bin/pangenome_report_render.py"` |
| locus view counts and ranks | all matrix strains | `"bin/pangenome_island_loci.py"` |

Inconsistencies found:

1. Outgroup polarity (`direction_a*`) uses all outgroup strains, but `bin_out` uses
   outgroup reps only. Near-duplicate outgroup strains weight the polarity fraction.
2. The co-occurrence test uses ingroup reps. Pair classification, which gates `trans`,
   uses all strains, including outgroup strains and non-reps.
3. The accumulation curve and Heaps fit use all strains, including the outgroup and
   non-reps. The frequency bins use ingroup reps. The curve therefore mixes two groups
   and duplicated genomes.
4. "Accessory" has three definitions: islands and neighbourhood use not in {core,
   soft_core}; assembly QC uses not `core`; Pfam background uses {shell, cloud}.
5. Islands are built in outgroup strains with the ingroup core definition.
6. The Dollo carrier set includes non-reps, but the Dollo ingroup set excludes them.
   A gain clade that contains a non-rep ingroup strain is not labelled `gain`.

---

## 5. Parameters

"Doc" = `"docs/pangenome-assumptions.md"` (section number). "Config" = the inline comment
in `"nextflow.config"`.

| param | default | step | effect | justification |
|---|---|---|---|---|
| `pangenome_ingroup_label` / `_outgroup_label` | `IN` / `OUT` | 3.1 | group roles | label (Doc sec. 10) |
| `pangenome_id_sep` | `|` | 3.2 | ID prefix separator | label (Doc sec. 10) |
| `pangenome_cluster_backend` | `mmseqs` | 3.3 | clustering tool | data: ADR-0003 ARI 0.94 on a 5-strain subset (Doc sec. 1) |
| `pangenome_tier1_min_id` | 0.9 | 3.3 | family identity cutoff | negative result: sweep shows it does not move the N50 confound (Doc sec. 1) |
| `pangenome_tier1_cov` | 0.8 | 3.3 | family coverage cutoff | negative result, plus a mechanism against it (Doc sec. 1) |
| `pangenome_rescue_enable` | true | 3.6 | run TBLASTN rescue | theory+practical; needs structural filter (Doc sec. 2) |
| `pangenome_rescue_evalue` | 1e-10 | 3.6 | tblastn e-value | no justification found (Doc sec. 2) |
| `pangenome_rescue_max_target_seqs` | 200 | 3.6 | tblastn hit cap | convention with stated reason (Doc sec. 2) |
| `pangenome_rescue_min_pident` | 90.0 | 3.6 | HSP identity | negative result (Doc sec. 2) |
| `pangenome_rescue_min_qcov` | 80.0 | 3.6 | query coverage | negative result (Doc sec. 2) |
| `pangenome_rescue_structural_filter` | true | 3.6 | overlap/short/hotspot filter | data: 77.6% of rescuable cells overlap another family (Doc sec. 2) |
| `pangenome_rescue_min_rep_length` | 150 | 3.6 | min rep length (aa) | data (Doc sec. 2) |
| `pangenome_rescue_hotspot_window` | 2000 | 3.6 | hotspot bin size (bp) | data (Doc sec. 2) |
| `pangenome_rescue_hotspot_min_families` | 3 | 3.6 | hotspot family count | no justification found for the value (Doc sec. 2) |
| `pangenome_strict` | false | 3.19 | fail on diagnostics | convention; not in Doc |
| `pangenome_rescue_redundancy_threshold` | 0.5 | 3.19 | diagnostic trip | no justification found; not in Doc |
| `pangenome_dereplicate` | true | 3.7 | use reps as denominator | structural (Doc sec. 3) |
| `pangenome_mash_threshold` | 0.001 | 3.7 | dedup distance | data: plateau below 0.002 (Doc sec. 3) |
| `pangenome_assign_clades` | true | 3.8 | fill TaxonGroup | structural; known failure on single species (Doc sec. 3) |
| `pangenome_clades_k_range` | `2,10` | 3.8 | k sweep | data for range only (Doc sec. 3) |
| `pangenome_clades_k_fixed` | null | 3.8 | fixed k | override (Doc sec. 3) |
| `pangenome_clades_n_pcoa` | 4 | 3.8 | PCoA axes | data (Doc sec. 3) |
| `pangenome_core_cutoff` | 0.95 | 3.9 | core bin | convention; not measured (Doc sec. 4) |
| `pangenome_softcore_cutoff` | 0.90 | 3.9 | soft-core bin | convention (Doc sec. 4) |
| `pangenome_shell_cutoff` | 0.15 | 3.9 | shell/cloud split | convention (Doc sec. 4) |
| `pangenome_outgroup_min_bin_strains` | 3 | 3.9 | bin outgroup | judgement (Doc sec. 4) |
| `pangenome_min_strain_count` | 5 | 3.10 | co-occurrence eligibility | no justification found (Doc sec. 5) |
| `pangenome_fdr_alpha` | 0.05 | 3.10 | BH threshold | convention (Doc sec. 5) |
| `pangenome_screen_alpha` | 0.2 | 3.10 | prefilter | structural (Doc sec. 5) |
| `pangenome_polarity_loss_min_frac` | 0.9 | 3.10 | loss call | data: Coccidioides calibration (Doc sec. 5) |
| `pangenome_polarity_gain_max_frac` | 0.1 | 3.10 | gain call | no justification found; symmetric choice (Doc sec. 5) |
| `pangenome_species_tree` | null | 3.10 | Dollo columns | optional input; not in Doc |
| `pangenome_captain_hmm` / `_hmm_name` / `pangenome_pfam_hmm` | null | 3.12 | captain search | structural (Doc sec. 7) |
| `pangenome_captain_hmm_evalue` | 1e-3 | 3.12 | hmmsearch -E | convention (Doc sec. 7) |
| `pangenome_pair_class_k` | 10 | 3.13 | linkage window (genes) | no justification found (Doc sec. 6) |
| `pangenome_pair_class_physical_threshold` | 0.5 | 3.13 | physical call | no justification found (Doc sec. 6) |
| `pangenome_pair_class_trans_threshold` | 0.05 | 3.13 | trans call | no justification found (Doc sec. 6) |
| `pangenome_pair_class_min_co_carrying` | 5 | 3.13 | min co-carriers | no justification found (Doc sec. 6) |
| `pangenome_pair_class_perm_alpha` | 0.05 | 3.13 | stratified p gate | convention (Doc sec. 6) |
| `pangenome_pair_class_min_clades` | 2 | 3.13 | min clade labels | Doc sec. 6 cites an FPR sweep, but that sweep is about stratification depth, not this gate (see sec. 7) |
| `pangenome_island_pfam_hmm` | null | 3.17 | enable Pfam branch | structural (Doc sec. 7) |
| `pangenome_island_min_size` | 2 | 3.14 | min island size | no justification found (Doc sec. 7) |
| `pangenome_pfam_domain_evalue` | 1e-3 | 3.17 | domain i-Evalue | convention (Doc sec. 7) |
| `pangenome_pfam_chunk_size` | 3000 | 3.17 | scan chunking | compute (Doc sec. 7) |
| `pangenome_marker_names` / `_hmm_paths` | `''` | 3.14 | marker searches | structural (Doc sec. 7) |
| `pangenome_marker_evalue` | 1e-5 | 3.14 | hmmsearch -E | convention (Doc sec. 7) |
| `pangenome_pfam2go` | null | 3.17 | GO mapping | structural (Doc sec. 7) |
| `pangenome_accumulation_permutations` | 20 | 3.20 | curve orders | no justification found (Doc sec. 7) |
| `pangenome_accumulation_seed` | 0 | 3.20 | seed | structural (Doc sec. 7) |
| `pangenome_leiden_resolution` | 1.0 | 3.15 | module granularity | data says too coarse at genus scale (Doc sec. 8) |
| `pangenome_leiden_seed` | 0 | 3.15 | seed | structural (Doc sec. 8) |
| `pangenome_module_min_size` | 2 | 3.15/3.16 | min module size | no justification found (Doc sec. 8) |
| `pangenome_neighborhood_max_kb` | 100 | 3.16 | colocalization distance | data: same calls at 25-200 kb (Config; spec 2026-09-19 B1 pass) |
| `pangenome_neighborhood_min_gene_gap` | 11 | 3.16 | exclude adjacent pairs | theory: outside k=10; not swept (Config) |
| `pangenome_neighborhood_n_perm` | 200 | 3.16 | null draws | no justification found |
| `pangenome_neighborhood_seed` | 0 | 3.16 | seed | structural |
| `pangenome_neighborhood_null_pool` | `accessory` | 3.16 | null pool | data: all-genes null gave spurious 1.2-1.3x (Config; spec 2026-09-19) |
| `pangenome_qc_terminus_window_bp` | 1000 | 3.19 | terminus window | no justification found ("this script's own choice", script docstring) |
| `pangenome_qc_rho_warn_threshold` | 0.3 | 3.19 | rho warning | no justification found |
| `pangenome_top_islands_min_strains` | 2 | 3.21/3.22 | min carriers for display and locus candidates | data: 62% single-strain islands (Doc sec. 9) |
| `pangenome_viz_top_islands` | 50 | 3.24 | islands drawn | payload size (Doc sec. 9) |
| `pangenome_locus_flank` / `_flank_min` | 5 / 3 | 3.22 | flank genes | no justification found (Doc sec. 9) |
| `pangenome_locus_k` | 10 | 3.22 | in-place window | no justification found (Doc sec. 9) |
| `pangenome_locus_empty_frac` | 0.8 | 3.22 | empty-site rule | no justification found (Doc sec. 9) |
| `pangenome_locus_containment` | 0.5 | 3.22 | locus grouping | no justification found (Doc sec. 9) |
| `pangenome_locus_rank` / `pangenome_top_loci` | `presence` / 100 | 3.22 | page sort / drawn count | design choice (Doc sec. 9) |
| `pangenome_locus_per_rank` | 20 | 3.22 | per-rank winners | no justification found (Doc sec. 9) |
| `pangenome_locus_poly_min_strains` / `_min_frac` / `_max_frac` | 20 / 0.05 / 0.95 | 3.22 | rank "within" | no justification found (Doc sec. 9) |
| `pangenome_locus_fixed_diff` | 0.95 | 3.22 | rank "species" | no justification found (Doc sec. 9) |
| `pangenome_locus_dna_check` | true | 3.23 | DNA check | data: empty sites were model differences in spot checks (Doc sec. 9) |
| `pangenome_locus_dna_min_id` / `_min_cov` | 90 / 80 | 3.23 | DNA present rule | no justification found (Doc sec. 9) |
| `pangenome_locus_dna_batch` | 50 | 3.23 | batching | compute (Doc sec. 9) |
| `pangenome_locus_candidates` | 200 | 3.22 | loci scored | compute (Doc sec. 9) |
| `pangenome_clinker` | true | 3.24 | clinker panel | structural |
| `pangenome_clinker_max_strains` | 12 | 3.24 | strains per figure | compute (Doc sec. 9) |
| `pangenome_clinker_slim` | true | 3.24 | strip sequences | data: render test (Doc sec. 9) |
| `pangenome_clinker_batch` | 50 | 3.24 | batching | compute (Doc sec. 9) |
| `pangenome_clinker_max_gap` | 20000 | 3.24 | split regions | no justification found; round number (Doc sec. 9) |

Hard-coded values (not params): species-rank minimum n_s = 10 (`"lib/island_locus.py"`
`score_species`); informative minimums 10 empty / 2 full; modified-z cutoff 3.5; Leiden
weight cap 300; GFF3 hard-error match fraction 0.50 and warn 0.02; unstratified warning at
> 50% empty TaxonGroup; module top 15 domains.

---

## 6. Statistical caveats

1. **`permutation_p` is not corrected for multiple testing.** It is computed only on BH
   survivors. Pair classification then gates `trans` on raw `permutation_p < 0.05`.
   The pairs tested can number in the thousands or more. No second correction runs.
2. **The stratified p-value is exact, not a permutation estimate.** The code computes a
   sum of independent per-clade hypergeometric variables. The column name says
   "permutation". The methods text should say "exact within-clade conditional test".
3. **Strata can mix label schemes.** The stratum is the `TaxonGroup` string. Curated
   labels and Mash labels (`clade_<k>`) can co-exist after `FILL_TAXON_GROUP`. Empty
   labels form one `"unknown"` stratum. The code only warns when > 50% of the tested
   strains are unlabelled. The `min_clades` gate counts `"unknown"` as a clade.
4. **The `min_clades` gate counts clades of family A carriers only.** It ignores family B.
5. **Mash clades on one species can be noise.** Doc sec. 3 reports split-half ARI
   0.02-0.26 within one species. The code has no automated check; the diagnostic
   `clade_validity` is `not_computed`.
6. **Pseudo-replication.** Dereplication uses mash distance < 0.001 with single linkage.
   Strains above that distance but still near-identical remain separate reps. The
   Fisher test and BH treat every rep as independent. Clade stratification is the only
   control. Outgroup polarity, pair classification, islands, the accumulation curve and
   the locus view use all strains, including non-reps.
7. **Island "significance" is a label, not a test.** An island is "significant" when it
   contains one pair classified `starship_explained`, `unexplained_physical` or
   `ambiguous_linkage`. `ambiguous_linkage` (0.05 < linkage < 0.5) qualifies. No island
   p-value exists. The pair-level BH applies to co-occurrence, not to linkage.
8. **Physical linkage has no null.** `linkage_fraction` is compared with fixed thresholds
   (0.5, 0.05). It is computed over all strains, not the tested strains.
9. **Pfam enrichment is pooled.** One test per domain compares the union of all island
   families with the shell+cloud background. It tests only domains seen in island
   families. Families, not genes or strains, are the unit. Island families are part of
   the background.
10. **Only positive association is tested.** The Fisher test is one-sided (`greater`).
    Mutually exclusive families are not detected.
11. **The screen is assumed conservative.** The code states that the normal-approximation
    screen never discards a pair that the exact test would keep at `screen_alpha`. No
    proof or test of that claim was found [UNVERIFIED]. BH itself stays valid because
    screened pairs enter as p = 1.
12. **Neighbourhood p-values have a floor and no correction.** The floor is 1/201. The
    design spec (2026-09-19, B1 pass) notes that effects of 1.02-1.04x reach that floor.
    No correction across modules.
13. **`genome_only` equals `present` in all tests.** Rescue calls count as presence in
    bins, Fisher tests, polarity and linkage.
14. **Openness rule.** The code labels the pangenome open when the fitted exponent of pan
    size vs N is < 1. Pan size cannot fall as N grows, so gamma >= 0. A value >= 1 means
    linear or faster growth. The rule therefore labels almost every curve "open". The
    usual Heaps' law reading (Tettelin et al. 2008) is on a different quantity or a
    different threshold [UNVERIFIED literature convention]. The fit also includes the
    outgroup and non-rep strains.
15. **Assembly-quality rho has no p-value.** The warning uses |rho| > 0.3 only.
16. **Locus ranks are scores, not tests.** Thresholds (10 losses, 2 carriers, 0.95
    difference) are fixed. DNA presence uses identity and coverage cutoffs only, with the
    blastn default e-value.
17. **Leiden output depends on resolution and seed.** Doc sec. 8 reports that 1.0 gives
    4-5 giant modules at genus scale. Edge weights use `fdr_q`, which is rounded to three
    significant digits in the input file.

---

## 7. Discrepancies and items to verify

Docs vs code:

- Doc sec. 5 lists `pangenome_pair_class_min_clades` as "the co-occurrence null's
  stratification depth". In code it is not. The stratification uses every distinct
  `TaxonGroup` label (`exact_stratified_pvalue`). `min_clades` is a separate gate on
  the count of labels among family A carriers.
- Doc sec. 6 and sec. 11 item 3 record a decision to stop using `min_clades` as a `trans`
  filter when clades are noise. Code always applies it (`classify_pair`).
- Doc sec. 3 says `pangenome_dereplicate` "does not gate any compute today". Code skips
  `MASH_SKETCH_ALL` and `DEREPLICATE` when it is false.
- Doc sec. 4 says bins are a "fraction of ingroup strains". Code uses ingroup reps when
  dereplication is on (as `"nextflow.config"` says).
- Doc does not list `pangenome_species_tree`, `pangenome_neighborhood_*`,
  `pangenome_strict`, `pangenome_rescue_redundancy_threshold`,
  `pangenome_qc_terminus_window_bp`, `pangenome_qc_rho_warn_threshold`.
- Doc sec. 11 item 4 says no step distinguishes `genome_only`. The locus view does
  (rescue cell codes). All statistics do not.
- `"docs/superpowers/specs/2026-09-28-pangenome-group-bins-design.md"` says "Nothing here
  is implemented". The code implements it (`bin_out`, `nonrep_only`, `outgroup_only`).
- `"docs/superpowers/specs/2026-09-28-outgroup-evidence-coverage-design.md"` concerns the
  novelty pipeline (`main.nf`), not `pangenome.nf`. Not used here.
- Pair classification docstring says captain evidence needs "at least one co-carrying
  strain". Code accepts any strain that has positions for A or B and a captain within k.
- `"modules/pangenome/prefix_and_cluster.nf"` comment says diamond is "Still not validated
  against a real multi-strain run". `"pangenome.nf"` comment and Doc sec. 1 say it was
  validated (ARI 0.94). The two comments disagree.
- Report render docstring says "gamma >= 1 indicates a closed one". See caveat 14.

Code issues to verify (read in code; not run):

- `EXTRACT_RESCUE_POSITIONS` picks the best-bitscore HSP without the structural filter.
  The position of a `genome_only` cell can come from an HSP that `RESCUE_PASS` rejected.
- `DEREPLICATE` builds `dna_paths` from every samplesheet row with a DNA value, but the
  mash matrix holds only ingroup and outgroup strains. A row with another `GROUP` and a
  DNA file would raise a `KeyError` at `dedup_group_id[short]` (inferred, not run).
- Fixed paths `<data_dir>/pep` and `<data_dir>/dna` (section 2.3) disagree with the flat
  layout that `"pangenome.nf"` accepts. A flat `data_dir` would fail in these steps
  (inferred, not run).
- `REPORT_RENDER`, `ISLAND_LOCI` and `ISLAND_SYNTENY` receive the original samplesheet,
  not the clade-filled one. Only `Species` and `GROUP` are read there, so the effect is
  likely nil [UNVERIFIED].
- `"pangenome.nf"` filters `GROUP` without the alias map that `parse_config` applies.
- The mash sketch size and k-mer size are not set (mash defaults) [UNVERIFIED values].
- The isoform-collapse requirement in `"bin/pangenome_cluster_backend.py"` is not met by
  any wired step. Multi-isoform proteomes would inflate copy numbers and paralog counts
  [UNVERIFIED effect].
- `choose_drawn` can return more than `top_loci` loci when 4 x `per_rank` > `top_loci`.
  The validation in `"pangenome.nf"` checks only `per_rank <= top_loci`.
- Islands dedup key ignores gene order; two strains with the same family set in a
  different order merge. `island_size` comes from the first strain only.
- `n_strains` for an island can count one strain twice if it holds two islands with the
  same key.
- I could not determine the mmseqs cluster mode and sensitivity defaults, or clinker's
  link-scoring method, from this repository [UNVERIFIED].
