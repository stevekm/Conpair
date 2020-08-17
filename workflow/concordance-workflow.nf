// workflow for running concordance in parallel for all tumors vs. a single list of normals
nextflow.enable.dsl=2

include { run_concordance } from './run-concordance.nf'

log.info("----------------")
log.info("workflow params:")
log.info("${params}")
log.info("----------------")

workflow {
    // need to load all the tumor file paths
    tumors = Channel.fromPath("${params.tumors_list}").splitCsv().map { [file(it[0]), file("${params.normals_list}"), file("${params.markers_txt}")] }

    run_concordance(tumors)

    run_concordance.out.benchmarks | collectFile(name: 'benchmarks.tsv', storeDir: "${params.output_dir}")
    run_concordance.out.concordance_vals | collectFile(name: 'concordance.tsv', keepHeader: true, storeDir: "${params.output_dir}") 
}
