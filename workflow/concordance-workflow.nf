// workflow for running concordance in parallel for all tumors vs. a single list of normals
nextflow.enable.dsl=2

include { run_concordance } from './run-concordance.nf'
include { plot_concordance_distribution } from './plot_concordance_distribution.nf'
include { plot_benchmarks } from './plot_benchmarks.nf'
include { filter_concordance } from './filter_concordance.nf'

log.info("----------------")
log.info("workflow params:")
log.info("${params}")
log.info("----------------")

workflow {
    // need to load all the tumor file paths
    tumors = Channel.fromPath("${params.tumors_list}").splitCsv().map { [file(it[0]), file("${params.normals_list}"), file("${params.markers_txt}")] }

    run_concordance(tumors)

    run_concordance.out.benchmarks | collectFile(name: 'benchmarks.tsv', storeDir: "${params.output_dir}") | set { benchmarks_tsv }
    run_concordance.out.concordance_vals | set { tumor_concordance_vals }

    tumor_concordance_vals | collectFile(name: 'concordance.all.tsv', keepHeader: true, storeDir: "${params.output_dir}") | set { concordance_tsv }

    plot_concordance_distribution(concordance_tsv)
    plot_benchmarks(benchmarks_tsv)
    filter_concordance(tumor_concordance_vals)

    filter_concordance.out.filtered_concordance | collectFile(name: 'concordance.filtered.tsv', keepHeader: true, storeDir: "${params.output_dir}") 
}
