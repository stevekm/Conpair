#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Load the genotype likelihoods for all the samples in the list and save them to Python pickle for faster loading
"""
import os
import sys
import pickle
import argparse

# need to import the module from the other dir
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from modules.ContaminationMarker import get_markers, genotype_likelihoods_for_markers
sys.path.pop(0)

# need to find a default set of targets to use; Conpair-GRCh37-default
default_marker_file = os.path.join(THIS_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')



def main(**kwargs):
    """
    """
    pileup_file = kwargs.pop('pileup_file', 'normals.txt')
    output_dir = kwargs.pop('output_dir', 'output')
    markers = kwargs.pop('markers', default_marker_file)
    min_base_quality = kwargs.pop('min_base_quality', 20)
    min_mapping_quality = kwargs.pop('min_mapping_quality', 10)

    # read paths to pileup files from list
    with open(pileup_file) as fin:
        all_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    Markers = get_markers(markers)

    for pileup_file in all_pileups:
        likelihoods = genotype_likelihoods_for_markers(Markers, pileup_file, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)

        output_file = os.path.join(output_dir, os.path.basename(pileup_file)) + '.pickle'
        if os.path.exists(output_file):
            print("ERROR: file already exists: " + output_file)
            raise
        else:
            print("saving file: " + output_file)
            with open(output_file,"wb") as fout:
                pickle.dump(likelihoods, fout)

    # with open("likelihood.pickle","rb") as fin:
    #     loaded_genotype_likelihoods = pickle.load(fin)
    #
    # print(len(likelihoods), len(loaded_genotype_likelihoods))
    # print(type(likelihoods), type(loaded_genotype_likelihoods))
    # # print(likelihoods.values())
    # for key, value in likelihoods.items():
    #     print(value == likelihoods[key])

def parse():
    """
    Parse the command line options
    """
    parser = argparse.ArgumentParser(description = 'Generate cBio Portal metadata files from various input files')
    parser.add_argument('--pileup-file', dest = 'pileup_file', default = "normals.txt", help = 'File with a list filepaths to the pileups of the tumor samples to use')
    parser.add_argument('--output-dir', dest = 'output_dir', default = 'output', help = 'Output location for files')
    parser.add_argument('--markers', dest = 'markers', default = default_marker_file, help = 'Markers to use for analysis')
    parser.add_argument('--min-base-quality', dest = 'min_base_quality', type = int, default = 20, help = 'Minimum base quality to use in output')
    parser.add_argument('--min-mapping-quality', dest = 'min_mapping_quality', type = int, default = 10, help = 'Minimum mapping quality to use in output')

    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
