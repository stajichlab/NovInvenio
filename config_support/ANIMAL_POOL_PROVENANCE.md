# `animal_pool.csv` provenance and regeneration

`config_support/animal_pool.csv` is a `--extra-pool` for
`bin/build_targeted_configs.py` (same schema as the main master pool:
`Species,Strain,ProteinPath,DNAPath,Lineage,NCBI_TaxID`), covering species
outside `Fungi_BFD`'s scope for the chytrid-vs-animal shared-gene comparison
(see `configs/batches/chytrid_animal_v1.yaml`).

The actual protein FASTAs it points at live under `db/animal_outgroup/`,
which is git-ignored (large reference data, not source — same convention as
`db/modelorgs/*.fasta`). Regenerate them with:

```bash
mkdir -p db/animal_outgroup && cd db/animal_outgroup

# Real RefSeq assemblies used (2026-09-05):
#   Drosophila melanogaster: GCF_000001215.4_Release_6_plus_ISO1_MT
#   Caenorhabditis elegans:  GCF_000002985.6_WBcel235
#   Mus musculus:            GCF_000001635.27_GRCm39
#   Monosiga brevicollis:    GCF_000002865.3_V1.0  (closest living non-animal
#                            relative of animals -- included so a signal
#                            shared with animals specifically can be told
#                            apart from one shared with the broader
#                            Opisthokonta/holozoan clade animals sit inside)

for spec in \
    "Drosophila_melanogaster GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT" \
    "Caenorhabditis_elegans GCF/000/002/985/GCF_000002985.6_WBcel235" \
    "Mus_musculus GCF/000/001/635/GCF_000001635.27_GRCm39" \
    "Monosiga_brevicollis GCF/000/002/865/GCF_000002865.3_V1.0" \
; do
    set -- $spec
    name=$1; path=$2; acc=$(basename "$path")
    base="https://ftp.ncbi.nlm.nih.gov/genomes/all/$path/$acc"
    curl -sL "${base}_protein.faa.gz" | gunzip > "${name}.proteins.fa"
    curl -sL "${base}_feature_table.txt.gz" -o "${name}_feature_table.txt.gz"
done

# Collapse each to one representative (longest) protein per gene, using the
# feature table's real gene<->protein mapping -- NOT a guessed FASTA-header
# heuristic (RefSeq header text/format is not a reliable isoform-grouping
# key across species/annotation pipelines). See bin/collapse_isoforms.py.
for name in Drosophila_melanogaster Caenorhabditis_elegans Mus_musculus Monosiga_brevicollis; do
    ../../bin/collapse_isoforms.py \
        --feature-table "${name}_feature_table.txt.gz" \
        --protein-fasta "${name}.proteins.fa" \
        --output "${name}.proteins.collapsed.fa"
done
```

`animal_pool.csv`'s `ProteinPath` columns point at the `.proteins.collapsed.fa`
outputs. `DNAPath` is deliberately blank for all four species — this pathway
is protein-level presence/absence, not TBLASTN genome validation (which the
schema already permits leaving unset, same as several `Fungi_BFD` species with
no `DNA` in `configs/samples.csv`), and downloading full genome assemblies
(the mouse genome alone is gigabytes) is out of scope for what this
comparison actually needs.

`Lineage` values are real animal/holozoan taxonomy filled into the same
7-field `PHYLUM;SUBPHYLUM;CLASS;SUBCLASS;ORDER;FAMILY;GENUS` slots
`lib/lineage.py`'s `RANK_NAMES` expects (blank where a rank isn't confidently
assigned, same convention as the fungal master pool) — this is what makes
`load_master_pool()` accept the file at all, not an attempt to place these
species within Fungi_BFD's own taxonomy. Because `PHYLUM` always differs from
every fungal focal's, `lib/targeted_selection.py::candidate_pool()` will
never surface them as a `mode: nearest`/`trait` candidate (by construction,
not a bug) — reference them with `mode: explicit` instead.
