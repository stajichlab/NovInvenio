/**
 * Pipeline helper methods available to all Nextflow process closures.
 * Groovy files in lib/ are automatically loaded by Nextflow.
 */
class Helpers {
    /**
     * Return the project name for output directory construction.
     * Uses params.project if set; otherwise derives it from the config CSV basename.
     * Falls back to 'output' only if neither is available (should never happen in normal runs).
     */
    static String projectName(params) {
        if (params.project) return params.project
        if (params.config)  return new File(params.config.toString()).name.replaceFirst(/\.[^.]+$/, '')
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
}
