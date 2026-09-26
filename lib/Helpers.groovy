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
