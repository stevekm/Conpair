#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.15
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

import sys
import os
import optparse
import math
from collections import defaultdict
from ContaminationMarker import genotype_likelihoods_for_markers

def main(
    tumor_pileup,
    normal_pileup,
    markers_data,
    min_mapping_quality = 10,
    normal_homozygous_markers_only = False,
    min_cov = 10,
    min_base_quality = 20
    ):
    """
    Main control function for the script

    #     parser.add_option('-T', '--tumor_pileup', help='TUMOR PILEUP FILE [mandatory field]', action='store')
    #     parser.add_option('-N', '--normal_pileup', help='NORMAL PILEUP FILE [mandatory field]', action='store')
    #     parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
    #     parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [default: 10]', default=10, type='int', action='store')
    #     parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY [default: 20]', default=20, type='int', action='store')
    #     parser.add_option('-H', '--normal_homozygous_markers_only', help='USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)', default=False, action='store_true')
    """
    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(markers_data, normal_pileup, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)
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
