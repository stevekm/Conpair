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
from ContaminationMarker import get_markers, genotype_likelihoods_for_markers

# need to import the module from the other dir
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)

# need to find a default set of targets to use
default_marker_file = os.path.join(PARENT_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

def main(
    tumor_pileup,
    normal_pileup,
    min_mapping_quality = 10,
    normal_homozygous_markers_only = False,
    markers = default_marker_file,
    min_cov = 10,
    min_base_quality = 20
    ):
    """
    Main control function for the script

    #     parser.add_option('-T', '--tumor_pileup', help='TUMOR PILEUP FILE [mandatory field]', action='store')
    #     parser.add_option('-N', '--normal_pileup', help='NORMAL PILEUP FILE [mandatory field]', action='store')
    #     parser.add_option('-M', '--markers', help='MARKER FILE [Conpair-GRCh37-default]', action='store')
    #     parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
    #     parser.add_option('-O', '--outfile', help='TXT OUTPUT FILE [stdout by default]', default="-", type='string', action='store')
    #     parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [default: 10]', default=10, type='int', action='store')
    #     parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY [default: 20]', default=20, type='int', action='store')
    #     parser.add_option('-H', '--normal_homozygous_markers_only', help='USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)', default=False, action='store_true')
    """
    Markers_data = get_markers(markers)

    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers_data, normal_pileup, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)
    Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(Markers_data, tumor_pileup, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)

    concordant = 0
    discordant = 0
    for m in Markers_data:
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
    num_total_markers = len(Markers_data)
    return(concordance, num_markers_used, num_total_markers)
