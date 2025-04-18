process likelihoods {
    // make a Python Pickle of the genotype likelihoods for faster loading in Conpair
    // NOTE: be careful which version of Python is used here; "python" might be Python 3 on newer systems by default
    publishDir "${params.output_dir}/likelihoods", mode: 'copy'

    input:
    tuple path(pileup), path(markers)

    output:
    path "${output_file}"

    script:
    output_file = "${pileup}".replaceFirst(/.pileup$/, ".pickle")
    """
    make_genotype_likelihoods.py \
    --pileup "${pileup}" \
    --markers "${markers}"
    """
}
