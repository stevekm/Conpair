process filter_concordance {
    // keep only high concordance values
    publishDir "${params.output_dir}/concordance", mode: 'copy'

    input:
    path(concordance_tsv)

    output:
    path "${output_file}", emit: filtered_concordance

    script:
    output_file = "${concordance_tsv}".replaceFirst(/.tsv$/, ".filtered.tsv")
    """
    # save headers
    head -1 "${concordance_tsv}" > "${output_file}"

    # get all rows with concordance >0.50; skip non-numeric values, skip rows where num_markers_used < 10
    tail -n +2 "${concordance_tsv}" | \
    awk '{ if(\$1+0 == \$1) print \$0 }' | \
    awk '{if(\$2 > 10) print \$0 }' | \
    awk '{if(\$1 > 0.50) print \$0 }' >> "${output_file}"
    """
}
