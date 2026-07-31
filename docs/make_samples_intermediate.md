# Building a `config_support/` sample pool from a BFD `samples.csv`

`bin/make_config.py --samples` expects a master pool CSV in NovInvenio's own
format: `Species,Strain,Protein,DNA,Short,Lineage`. External taxonomy tables —
e.g. `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv`
(`ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,
SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,TRANSL_TABLE,LOCUSTAG`) —
use a different schema and will fail with a `KeyError` (e.g. on `'Species'`)
if passed to `make_config.py` directly.

## Converting with `bin/convert_bfd_samples.py`

`bin/convert_bfd_samples.py` bridges this gap:

- Reads the BFD `samples.csv`, optionally filtered by `--order` or
  `--class-name`.
- For each row, builds the funannotate output dirname as `<SPECIES>_<STRAIN>`
  (spaces → underscores) and looks for
  `--annotation-dir/<dirname>/predict_results/<dirname>.proteins.fa` /
  `.scaffolds.fa`. This naming convention was confirmed against real BFD
  `genome_annotation/` directories, e.g.
  `Anthracina_ramosa_CGMCC_3.16372/predict_results/...`.
- Rows with no matching annotation directory are skipped and reported on
  stderr — many genomes in a given clade may still be mid-annotation.
- Builds `Lineage` as `Fungi;PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;FAMILY;GENUS`
  (empty segments dropped).
- Generates a `Short` ID from genus+species prefixes, deduplicated on
  collision.
- `--link-dir data` symlinks the matched `.proteins.fa`/`.scaffolds.fa` files
  into `data/pep` and `data/dna` so `make_config.py --data-dir data` can
  validate them.

Output is an **intermediate sample pool**, not a run config — write
`--output` into `config_support/`, not `configs/`.

## Worked example: Chaetothyriales

1. First attempt filtered to just the target order:

   ```bash
   python3 bin/convert_bfd_samples.py \
       --input /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \
       --annotation-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation \
       --order Chaetothyriales \
       --link-dir data \
       --output config_support/Chaetothyriales_samples.csv
   ```

   Only 161 of 218 Chaetothyriales genomes matched (57 skipped — no
   `predict_results` directory yet). `make_config.py` then failed: with the
   pool filtered down to one order, there was no sibling clade left to
   auto-derive an outgroup from.

2. Re-ran with a broader pool so `make_config.py` could auto-derive outgroups
   (siblings within the same class):

   ```bash
   python3 bin/convert_bfd_samples.py \
       --input /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \
       --annotation-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation \
       --class-name Eurotiomycetes \
       --link-dir data \
       --output config_support/Chaetothyriales_samples.csv
   ```

   1100 samples matched across Eurotiomycetes.

3. Fed the pool into `make_config.py` to produce the actual run config:

   ```bash
   python3 bin/make_config.py \
       --ingroup-taxon Chaetothyriales \
       --max-per-outgroup-taxon 3 \
       --seed 42 \
       --output configs/Chaetothyriales.csv \
       --data-dir data \
       --samples config_support/Chaetothyriales_samples.csv
   ```

## Directory convention

- `configs/` — run configs (`--config` for `main.nf`) and the repo's own
  master pool (`configs/samples.csv`).
- `config_support/` — intermediate sample pools (e.g.
  `convert_bfd_samples.py` output) passed back into `make_config.py`'s
  `--samples`. Never passed directly as `--config` to `main.nf`.
