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
}
