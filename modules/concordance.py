#!/usr/bin/env python2.7
import sys
import os
import optparse
import math
from collections import defaultdict
from ContaminationMarker import genotype_likelihoods_for_markers
import pickle

def concordance(
    tumor_pileup,
    normal_pileup,
    markers_data,
    min_mapping_quality = 10,
    normal_homozygous_markers_only = False,
    min_cov = 10,
    min_base_quality = 20
    ):
    """
    Calculate the concordance between a tumor and a normal sample. Both tumor and normal can be loaded from either a standard GATK pileup, or from a pre-saved genotypes likelihoods Python .pickle file.

    Parameters
    ----------
    tumor_pileup: str
        path to tumor pileup file to use for concordance. Can be a standard GATK .pileup file, or a pre-saved genotypes likelihoods .pickle file
    normal_pileup: str
        path to normal pileup file to use for concordance. Can be a standard GATK .pileup file, or a pre-saved genotypes likelihoods .pickle file
    markers_data:
        data load for markers set from a call to `ContaminationMarker.get_markers`
    min_mapping_quality: int
        the min mapping quality to use
    normal_homozygous_markers_only: bool
        use only homozygous markers in the Normal sample
    min_cov: int
        the minimum coverage value to use
    min_base_quality: int
        the minimum base quality to use


    Returns
    -------
    (float, int, int)
        returns values for concordance, num_markers_used, num_total_markers based on the given pair

    Notes
    -----
    Trying to calculate concordance between two samples with no shared markers can cause a divide by zero ZeroDivisionError error; this is currently handled in the `run.py` script.
    TODO: figure out a good way to handle the ZeroDivisionError errors
    """
    if normal_pileup.endswith('.pickle'):
        with open(normal_pileup,"rb") as fin:
            Normal_genotype_likelihoods = pickle.load(fin)
    else:
        Normal_genotype_likelihoods = genotype_likelihoods_for_markers(markers_data, normal_pileup, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)

    if tumor_pileup.endswith('.pickle'):
        with open(tumor_pileup,"rb") as fin:
            Tumor_genotype_likelihoods = pickle.load(fin)
    else:
        Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(markers_data, tumor_pileup, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)

    concordant = 0
    discordant = 0
    for m in markers_data:
        NL = Normal_genotype_likelihoods[m]
        TL = Tumor_genotype_likelihoods[m]
        if NL is None or TL is None:
            continue
        if NL['coverage'] < min_cov or TL['coverage'] < min_cov:
            continue
        if normal_homozygous_markers_only:
            if NL['likelihoods'].index(max(NL['likelihoods'])) == 1:
                continue
        if NL['likelihoods'].index(max(NL['likelihoods'])) == TL['likelihoods'].index(max(TL['likelihoods'])):
            concordant += 1
        else:
            discordant += 1

    # TODO: find a better way to handle this
    # if concordant+discordant == 0:
    #     print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(min_cov))
    #     sys.exit(0)

    concordance = float(concordant)/float(concordant+discordant)
    num_markers_used = concordant + discordant
    num_total_markers = len(markers_data)
    return(concordance, num_markers_used, num_total_markers)
