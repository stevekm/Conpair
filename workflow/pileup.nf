process pileup {
    // make a GATK pileup from a bam file
    publishDir "${params.output_dir}/pileup", mode: 'copy'
    errorStrategy 'ignore' // see known errors listed below

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
    -rf BadCigar \
    --filter_reads_with_N_cigar \
    --filter_mismatching_base_and_quals
    """
}
//  ##### ERROR MESSAGE: SAM/BAM file SAMFileReader{....bam} is malformed: read starts with deletion. Cigar: 1D74M. Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar
// ##### ERROR MESSAGE: SAM/BAM file SAMFileReader{...duplex.bam} appears to be using the wrong encoding for quality scores: we encountered an extremely high quality score of 87; please see the GATK --help documentation for options related to this error
