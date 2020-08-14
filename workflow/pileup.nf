process pileup {
    // make a GATK pileup from a bam file
    publishDir "${params.output_dir}/pileup", mode: 'copy'

    input:
    tuple path(bam), path(bai), path(markers)

    output:
    path "${output_file}", emit: pileup

    script:
    output_file = "${bam}".replaceFirst(/.bam$/, ".pileup")
    """
    java \
    -Xmx12g \
    -jar "${params.gatk_jar}" \
    -T Pileup \
    -R "${params.ref_fasta}" \
    -I "${bam}" \
    -L "${markers}" \
    -o "${output_file}" \
    -verbose \
    -rf DuplicateRead \
    --filter_reads_with_N_cigar \
    --filter_mismatching_base_and_quals
    """
}
