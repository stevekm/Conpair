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

# need to import the module from the other dir
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from modules.ContaminationMarker import get_markers, genotype_likelihoods_for_markers
sys.path.pop(0)

def main(opts):
    """
    Main control function for the script
    """
    # get the parameters from the options passed
    min_mapping_quality = opts.get('min_mapping_quality', 10)
    normal_homozygous_markers_only = opts.get('normal_homozygous_markers_only', False)
    markers = opts.get('markers', None)
    outfile = opts.get('outfile', '-')
    min_cov = opts.get('min_cov', 10)
    min_base_quality = opts.get('min_base_quality', 20)
    tumor_pileup = opts.get('tumor_pileup', None)
    normal_pileup = opts.get('normal_pileup', None)

    CONPAIR_DIR = os.environ['CONPAIR_DIR']
    MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

    if markers:
        MARKER_FILE = markers

    Markers = get_markers(MARKER_FILE)
    COVERAGE_THRESHOLD = min_cov
    MMQ = min_mapping_quality
    MBQ = min_base_quality
    AA_BB_only = normal_homozygous_markers_only

    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, normal_pileup, min_map_quality=MMQ, min_base_quality=MBQ)
    Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, tumor_pileup, min_map_quality=MMQ, min_base_quality=MBQ)

    concordant = 0
    discordant = 0
    for m in Markers:
        NL = Normal_genotype_likelihoods[m]
        TL = Tumor_genotype_likelihoods[m]
        if NL is None or TL is None:
            continue
        if NL['coverage'] < COVERAGE_THRESHOLD or TL['coverage'] < COVERAGE_THRESHOLD:
            continue
        if AA_BB_only:
            if NL['likelihoods'].index(max(NL['likelihoods'])) == 1:
                continue
        if NL['likelihoods'].index(max(NL['likelihoods'])) == TL['likelihoods'].index(max(TL['likelihoods'])):
            concordant += 1
        else:
            discordant += 1

    # TODO: find a better way to handle this
    # if concordant+discordant == 0:
    #     print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(COVERAGE_THRESHOLD))
    #     sys.exit(0)

    concordance = float(concordant)/float(concordant+discordant)
    return(concordance)

def parse():
    """
    Parse command line options
    """
    desc = """Program to verify tumor-normal sample concordance"""
    parser = optparse.OptionParser(version='%prog version 0.15 3/August/2016', description=desc)
    parser.add_option('-T', '--tumor_pileup', help='TUMOR PILEUP FILE [mandatory field]', action='store')
    parser.add_option('-N', '--normal_pileup', help='NORMAL PILEUP FILE [mandatory field]', action='store')
    parser.add_option('-M', '--markers', help='MARKER FILE [Conpair-GRCh37-default]', action='store')
    parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
    parser.add_option('-O', '--outfile', help='TXT OUTPUT FILE [stdout by default]', default="-", type='string', action='store')
    parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [default: 10]', default=10, type='int', action='store')
    parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY [default: 20]', default=20, type='int', action='store')
    parser.add_option('-H', '--normal_homozygous_markers_only', help='USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)', default=False, action='store_true')

    (opts, args) = parser.parse_args()

    if not opts.tumor_pileup or not opts.normal_pileup:
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(opts.tumor_pileup):
        print('ERROR: Input tumor file {0} cannot be found.'.format(opts.tumor_pileup))
        sys.exit(1)

    if not os.path.exists(opts.normal_pileup):
        print('ERROR: Input normal file {0} cannot be found.'.format(opts.normal_pileup))
        sys.exit(1)

    if opts.markers:
        if not os.path.exists(opts.markers):
            print('ERROR: Marker file {0} cannot be find.'.format(opts.markers))
            sys.exit(2)

    print(main(vars(opts)))

if __name__ == '__main__':
    parse()
