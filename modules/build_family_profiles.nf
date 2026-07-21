// BUILD_FAMILY_PROFILES — turn ingroup gene families into a concatenated HMM database
// (ADR-0002 Q4/Q6). For each family with >= params.family_min_members members: align
// with famsa, build a profile with hmmbuild (HMM name = the family representative id, so
// the domtblout query name ties back to the cluster tsv). Singletons and sub-threshold
// families are dropped by extract_family_seqs.py. The per-family loop runs inside one
// process to avoid a Nextflow task per family (thousands at order scale).
process BUILD_FAMILY_PROFILES {
    label 'high_cpu'
    tag "build_profiles"
    publishDir { "${params.outdir}/${Helpers.projectName(params)}/families" }, mode: 'copy'

    input:
    path(cluster_tsv)
    path(ingroup_fa)

    output:
    path("family_profiles.hmm"), emit: profiles
    path("families.tsv"),        emit: families

    shell:
    '''
    if [ ! -s !{cluster_tsv} ]; then
        touch family_profiles.hmm families.tsv
        exit 0
    fi

    mkdir -p fams
    extract_family_seqs.py \
        --cluster-tsv !{cluster_tsv} \
        --fasta !{ingroup_fa} \
        --min-members !{params.family_min_members} \
        --outdir fams

    : > family_profiles.hmm
    for fa in fams/fam_*.faa; do
        [ -e "$fa" ] || continue
        base=$(basename "$fa" .faa)
        rep=$(awk -F'\\t' -v b="$base" '$1==b{print $2}' fams/families.tsv)
        famsa -t !{task.cpus} "$fa" "${fa}.aln" 2>/dev/null
        hmmbuild --cpu !{task.cpus} -n "$rep" "${fa}.hmm" "${fa}.aln" > /dev/null
        cat "${fa}.hmm" >> family_profiles.hmm
    done
    cp fams/families.tsv families.tsv
    '''
}
