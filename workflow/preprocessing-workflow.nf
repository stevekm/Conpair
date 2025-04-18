// Workflow for pre-processing of .bam files to generate pileups and genotype likelihood Python pickle files
nextflow.enable.dsl=2

include { pileup } from './pileup.nf'
include { likelihoods } from './likelihoods.nf'

log.info("----------------")
log.info("workflow params:")
log.info("${params}")
log.info("----------------")

workflow {
    input_files = Channel.fromFilePairs("${params.input_dir}/*{.bam,.bam.bai}") // files named foo.bam and foo.bam.bai
        // "${params.input_dir}/*{.bam,.bai}" <- use this is the files are named foo.bam and foo.bai ; the bam and bai come out in different order
        // def bai = items[0]
        // def bam = items[1]
        .map { label, items ->
            def bai = items[1]
            def bam = items[0]
            return( [ bam, bai, file("${params.markers_bed}") ] )
        }

    pileup(input_files)
    pileup.out.pileup | map { [ it, file("${params.markers_txt}") ] } | set { pileups_markers }
    likelihoods(pileups_markers)
}
