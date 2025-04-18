process run_concordance {
    publishDir "${params.output_dir}/concordance", mode: 'copy'

    input:
    tuple path(tumor_pileup), path(normals_list), path(markers)

    output:
    path "${output_file}", emit: concordance_vals
    path "${benchmarks_file}", emit: benchmarks

    script:
    output_file = "${tumor_pileup}.concordance.tsv"
    benchmarks_file = "${tumor_pileup}.benchmarks.tsv"
    """
    echo "${tumor_pileup}" > tumors.txt
    run.py concordance \
    --tumors-list tumors.txt \
    --normals-list "${normals_list}" \
    --markers "${markers}" \
    --output-file "${output_file}" \
    --threads "${task.cpus}" \
    --save-benchmarks \
    --benchmarks-file "${benchmarks_file}"
    """
}
