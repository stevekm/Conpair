process plot_concordance_distribution {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(concordance_tsv)

    output:
    path "${dist_pdf}"

    script:
    dist_pdf = "concordance_dist.pdf"
    """
    plot_concordance.R "${concordance_tsv}"
    """
}
