// SPECIES_TREE -- phyling (hmmsearch + hmmalign + treeness/RCV filter) + FastTree species
// tree over every proteome in the config, rooted on the OUT rows by bin/species_tree.py.
// Docs: docs/superpowers/specs/2026-10-17-species-tree-design.md
//
// Called from main.nf only when `species_tree` is true and `phyling_markerset_dir` is set;
// otherwise main.nf feeds tracked stub status files (assets/species_tree_stub) to the report.
// Never fails the run. The script always exits 0 and always writes all four outputs;
// species_tree.status.txt line 1 is `ok`, `skipped: <reason>` or `failed: <reason>`;
// an `ok` file also has `markers: <N>`, `rooting: method=... n_out=... n_out_outside=...
// n_non_out_inside=...` and `markerset: <name>` lines.
// Runs in phyling's own image (ghcr.io/stajichlab/phyling:<species_tree_container_version>).
// Image contents (phyling, FastTree, aster/ASTRAL, muscle, python3) are taken from the phyling
// Dockerfile. The image was NOT pulled or run when this was written; only phyling 2.5.0 installed
// with pip was run on real data. bin/species_tree.py needs only the python3 stdlib.
// Outside a container (no -profile docker/singularity) the script falls back to the `phyling/2.3`
// environment module when `module` exists. Not phyling/2.3.1: its `filter` subcommand fails
// ("No module specified").
// phyling_markerset_dir must be on a path the container can read (UCR: /bigdata is bound).
// Proteomes are staged as <Short>.faa so the Newick tip labels equal config Short.
process SPECIES_TREE {
    label 'med_cpu'
    tag "species_tree"
    container "ghcr.io/stajichlab/phyling:${params.species_tree_container_version}"
    publishDir "${params.outdir}/${Helpers.projectName(params)}", mode: 'copy'
    errorStrategy 'ignore'
    time '2.h'

    input:
    tuple val(shorts), path(proteomes, stageAs: 'prot_?/*')
    path(config_csv)

    output:
    path("species_tree.svg"),         emit: svg
    path("species_tree.rooted.nwk"),  emit: nwk
    path("species_tree.markers.tsv"), emit: markers
    path("species_tree.status.txt"),  emit: status

    script:
    def markerdir = params.phyling_markerset_dir ? "${params.phyling_markerset_dir}/${params.species_tree_markerset}" : ''
    """
    set +e
    : > species_tree.svg
    : > species_tree.rooted.nwk
    : > species_tree.markers.tsv
    status() { echo "\$1" > species_tree.status.txt; exit 0; }

    SHORTS=(${shorts.join(' ')})
    FILES=(${proteomes instanceof List ? proteomes.join(' ') : proteomes})
    [ "\${#SHORTS[@]}" -ge 4 ] || status "skipped: fewer than 4 species"
    MARKERDIR="${markerdir}"
    [ -d "\$MARKERDIR/hmms" ] || status "failed: markerset not found: \$MARKERDIR/hmms"

    # only when phyling is not already on PATH (it is, inside the container)
    command -v phyling >/dev/null 2>&1 || { command -v module >/dev/null 2>&1 && module load phyling/2.3 fasttree/2.1.11 >/dev/null 2>&1; }
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
    # Extra status lines for the report (line 1 stays the status word):
    #   markers: <N>   alignments (files or symlinks) kept by `phyling filter` (falls back to the requested N)
    #   rooting: method=<m> n_out=<n> n_out_outside=<k> n_non_out_inside=<m>   from bin/species_tree.py
    NM=\$(find phyling_filter -maxdepth 1 \\( -type f -o -type l \\) \\( -name '*.fa' -o -name '*.faa' -o -name '*.fas' -o -name '*.fasta' -o -name '*.aln' -o -name '*.mfa' \\) 2>/dev/null | wc -l)
    [ "\$NM" -gt 0 ] || NM=${params.species_tree_markers}
    ROOTING=\$(python3 -c "import json;d=json.load(open('species_tree.rooting.json'));print('method=%s n_out=%s n_out_outside=%s n_non_out_inside=%s'%(d['method'],d['n_out'],d['n_out_outside'],d.get('n_non_out_inside',0)))" 2>/dev/null)
    {
        echo ok
        echo "markers: \$NM"
        [ -z "\$ROOTING" ] || echo "rooting: \$ROOTING"
        echo "markerset: ${params.species_tree_markerset}"
    } > species_tree.status.txt
    exit 0
    """
}
