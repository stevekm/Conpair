nextflow.enable.dsl=2

include { pileup } from './pileup.nf'
include { likelihoods } from './likelihoods.nf'

log.info("----------------")
log.info("workflow params:")
log.info("${params}")
log.info("----------------")

workflow {
    input_files = Channel.fromFilePairs("${params.input_dir}/*{.bam,.bai}")
        .map { label, items ->
            def bai = items[0]
            def bam = items[1]
            return( [ bam, bai, file("${params.markers_bed}") ] )
        }

    pileup(input_files)
    pileup.out.pileup | map { [ it, file("${params.markers_txt}") ] } | set { pileups_markers }
    likelihoods(pileups_markers)
}
