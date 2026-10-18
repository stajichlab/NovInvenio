// SPECIES_TREE -- phyling (hmmsearch + hmmalign + treeness/RCV filter) + FastTree species
// tree over every proteome in the config, rooted on the OUT rows by bin/species_tree.py.
// Docs: docs/superpowers/specs/2026-10-17-species-tree-design.md
//
// Never fails the run. The script always exits 0 and always writes all four outputs;
// species_tree.status.txt line 1 is `ok`, `skipped: <reason>` or `failed: <reason>`.
// phyling is not in the container image. The process loads the `phyling/2.3` environment
// module when `module` exists (HPCC), else expects `phyling` on PATH. Not phyling/2.3.1: its
// `filter` and `tree` subcommands fail ("No module specified"), phyling/2.3 (2.3.1.post23) works.
// Proteomes are staged as <Short>.faa so the Newick tip labels equal config Short.
process SPECIES_TREE {
    label 'med_cpu'
    tag "species_tree"
    container "ghcr.io/stajichlab/novinvenio:${params.container_version}"
    publishDir "${params.outdir}/${Helpers.projectName(params)}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(shorts), path(proteomes, stageAs: 'prot_?/*')
    path(config_csv)

    output:
    path("species_tree.svg"),         emit: svg
    path("species_tree.rooted.nwk"),  emit: nwk
    path("species_tree.markers.tsv"), emit: markers
    path("species_tree.status.txt"),  emit: status

    script:
    def enabled   = (params.species_tree && params.species_tree.toString() != 'false') ? 'yes' : 'no'
    def markerdir = params.phyling_markerset_dir ? "${params.phyling_markerset_dir}/${params.species_tree_markerset}" : ''
    """
    set +e
    : > species_tree.svg
    : > species_tree.rooted.nwk
    : > species_tree.markers.tsv
    status() { echo "\$1" > species_tree.status.txt; exit 0; }

    SHORTS=(${shorts.join(' ')})
    FILES=(${proteomes instanceof List ? proteomes.join(' ') : proteomes})
    [ "${enabled}" = yes ] || status "skipped: disabled"
    [ "\${#SHORTS[@]}" -ge 4 ] || status "skipped: fewer than 4 species"
    MARKERDIR="${markerdir}"
    [ -n "\$MARKERDIR" ] || status "failed: phyling_markerset_dir is not set"
    [ -d "\$MARKERDIR/hmms" ] || status "failed: markerset not found: \$MARKERDIR/hmms"

    command -v module >/dev/null 2>&1 && module load phyling/2.3 >/dev/null 2>&1
    command -v phyling >/dev/null 2>&1 || status "failed: phyling not found"

    mkdir -p phyling_in
    for i in "\${!SHORTS[@]}"; do
        f="\${FILES[\$i]}"
        # phyling reads only plain or BGZF FASTA, so decompress plain gzip
        case "\$f" in
            *.gz) gzip -dc "\$f" > "phyling_in/\${SHORTS[\$i]}.faa" ;;
            *)    ln -s "\$(readlink -f "\$f")" "phyling_in/\${SHORTS[\$i]}.faa" ;;
        esac
    done

    {
        phyling align -I phyling_in -m "\$MARKERDIR/hmms" -o phyling_align -t ${task.cpus} &&
        phyling filter -I phyling_align -n ${params.species_tree_markers} -o phyling_filter -t ${task.cpus} &&
        phyling tree -I phyling_filter -M ft -o phyling_tree
    } > phyling.log 2>&1
    [ \$? -eq 0 ] || status "failed: phyling error: \$(grep -h ERROR phyling.log | tail -1 | cut -c1-200)"

    TREE=phyling_tree/final_tree.nw
    [ -s "\$TREE" ] || status "failed: phyling tree wrote no final_tree.nw"
    TSV=phyling_filter/treeness.tsv
    [ ! -s "\$TSV" ] || cp "\$TSV" species_tree.markers.tsv

    species_tree.py --newick "\$TREE" --config ${config_csv} \\
        --out_prefix species_tree 2> species_tree.err
    rc=\$?
    if [ \$rc -ne 0 ]; then
        : > species_tree.svg; : > species_tree.rooted.nwk
        status "failed: \$(tail -1 species_tree.err | cut -c1-200)"
    fi
    status "ok"
    """
}
