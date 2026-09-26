/**
 * Pipeline helper methods available to all Nextflow process closures.
 * Groovy files in lib/ are automatically loaded by Nextflow.
 */
class Helpers {
    /**
     * Return the project name for output directory construction.
     *
     * Fallback order: params.project, then params.pangenome_project, then the
     * basename of params.config (main.nf's --config), then the basename of
     * params.pangenome_samplesheet (pangenome.nf's --pangenome_samplesheet,
     * extension stripped the same way as --config), then the literal 'output'
     * if none of those are set.
     *
     * pangenome.nf used to try to plug this gap itself by assigning
     * `params.project = ...` at the top of its workflow{} block when neither
     * --project nor --pangenome_project was given. That assignment silently
     * no-ops on Nextflow 26.04.6: once a param is declared in
     * nextflow.config's `params {}` block (both `project` and
     * `pangenome_project` are, with default null), a later `params.x = y`
     * assignment in the script/workflow body does not persist -- a read of
     * params.x immediately afterward still sees the original config-declared
     * value. Undeclared params are unaffected (assigning a param
     * nextflow.config never mentions works normally). Verified with a
     * minimal local reproduction: a two-line nextflow.config (`params {
     * project = null }`) plus a workflow{} that assigns `params.project =
     * "world"` and immediately re-reads it comes back null. This is why a
     * real run without --pangenome_project printed "defaulting
     * params.project to 'null'" and then published under
     * <outdir>/output/pangenome/ instead of
     * <outdir>/<samplesheet-basename>/pangenome/ -- the assignment never
     * took, so every downstream publishDir/storeDir closure calling
     * projectName(params) fell through to the 'output' literal regardless of
     * what the warning printed. Fix: derive the samplesheet-basename
     * fallback here, inside projectName() itself, instead of relying on a
     * script-side mutation of params.project (issue #194).
     */
    static String projectName(params) {
        if (params.project)               return params.project
        if (params.pangenome_project)      return params.pangenome_project
        if (params.config)                 return new File(params.config.toString()).name.replaceFirst(/\.[^.]+$/, '')
        if (params.pangenome_samplesheet)  return new File(params.pangenome_samplesheet.toString()).name.replaceFirst(/\.[^.]+$/, '')
        return 'output'
    }

    /**
     * Return the absolute directory that docs/<project> reports should publish under:
     * the "docs" sibling of params.outdir, resolved against the launch directory if
     * params.outdir is relative. A plain relative "docs/..." publishDir instead resolves
     * against whatever directory `nextflow run` was launched from, which silently diverges
     * from params.outdir whenever a caller launches from an isolated subdirectory (see
     * the isolated-launch-dir pattern in the top-level run_*.sh scripts / README).
     */
    static String docsDir(params, launchDir) {
        def outdirFile = new File(params.outdir.toString())
        def outdirAbs = outdirFile.isAbsolute() ? outdirFile : new File(launchDir.toString(), outdirFile.path)
        return new File(outdirAbs.canonicalFile.parentFile, 'docs').path
    }

    /**
     * Boolean value of a params flag (issue #191). Nextflow 26 passes a command-line
     * `--flag false` to the script as the String "false", which is truthy in Groovy,
     * so `if (params.flag)` stays on. The same value from -params-file arrives as a
     * Boolean. The CLI value overrides the config after it is read, so this cannot be
     * fixed in nextflow.config; read every boolean param through this function.
     * A bare `--flag` arrives as the String "true".
     */
    static boolean asBool(v) {
        if (v == null) return false
        if (v instanceof Boolean) return v
        def s = v.toString().trim().toLowerCase()
        if (s in ['true', 't', 'yes', 'y', '1']) return true
        if (s in ['false', 'f', 'no', 'n', '0', '', 'null']) return false
        throw new IllegalArgumentException("not a boolean value: '${v}'")
    }
}
