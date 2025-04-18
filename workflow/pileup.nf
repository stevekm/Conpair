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
    // GATK4; https://github.com/nygenome/Conpair/issues/11
    """
    gatk Pileup \
    --reference "${params.ref_fasta}" \
    --input "${bam}" \
    --intervals "${markers}" \
    --output "${output_file}" \
    -verbose \
    -RF NotDuplicateReadFilter \
    -RF CigarContainsNoNOperator \
    -RF MatchingBasesAndQualsReadFilter \
    -RF GoodCigarReadFilter \
    -RF MappingQualityAvailableReadFilter \
    -RF MappingQualityNotZeroReadFilter
    """
    // GATK3; https://hub.docker.com/r/broadinstitute/gatk3/tags
    // """
    // java \
    // -Xmx12g \
    // -jar /usr/GenomeAnalysisTK.jar \
    // -T Pileup \
    // -R "${params.ref_fasta}" \
    // -I "${bam}" \
    // -L "${markers}" \
    // -o "${output_file}" \
    // -verbose \
    // -rf DuplicateRead \
    // -rf BadCigar \
    // --filter_reads_with_N_cigar \
    // --filter_mismatching_base_and_quals
    // """
}

// NOTE: some relevant errors from GATK3;
//  ##### ERROR MESSAGE: SAM/BAM file SAMFileReader{....bam} is malformed: read starts with deletion. Cigar: 1D74M. Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar
// ##### ERROR MESSAGE: SAM/BAM file SAMFileReader{...duplex.bam} appears to be using the wrong encoding for quality scores: we encountered an extremely high quality score of 87; please see the GATK --help documentation for options related to this error

// from GATK4 ; https://gatk.broadinstitute.org/hc/en-us/articles/360037591631-Pileup
// These Read Filters are automatically applied to the data by the Engine before processing by Pileup.
//     NotSecondaryAlignmentReadFilter
//     PassesVendorQualityCheckReadFilter
//     MappedReadFilter
//     NotDuplicateReadFilter
//     WellformedReadFilter
// https://gatk.broadinstitute.org/hc/en-us/articles/9570266920219--Tool-Documentation-Index#ReadFilters

// NOTE: samtools pileup format;
// https://samtools.sourceforge.net/pileup.shtml
// https://www.htslib.org/doc/samtools-mpileup.html


// example bad record with missing fields in GATK4;
//  1  5
//  2  154065383
//  3  A
//  4
//  5
//  6
//  7  0
//  8
