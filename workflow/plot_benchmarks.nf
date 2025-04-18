process plot_benchmarks {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(benchmarks_tsv)

    output:
    path "${aggr_table}"
    path "${time_pairs_threads_pdf}"
    path "${time_per_pair_pdf}"

    script:
    aggr_table = "aggregate_time_per_pair.tsv"
    time_pairs_threads_pdf = "time_pairs_threads.pdf"
    time_per_pair_pdf = "time_per_pair.pdf"
    """
    plot_benchmarks.R "${benchmarks_tsv}"
    """
}
