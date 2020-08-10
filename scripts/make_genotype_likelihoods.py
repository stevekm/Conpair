#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Load the genotype likelihoods for all the samples in the list and save them to Python pickle for faster loading
"""
import os
import sys
import pickle

# need to import the module from the other dir
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from modules.ContaminationMarker import get_markers, genotype_likelihoods_for_markers
sys.path.pop(0)

min_base_quality = 20
min_mapping_quality = 10
normal_pileup_file = "normal_pileups.txt"

with open(normal_pileup_file) as fin:
    all_normal_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

CONPAIR_DIR = os.environ['CONPAIR_DIR']
MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

Markers = get_markers(MARKER_FILE)

for pileup_file in all_normal_pileups:
    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, pileup_file, min_map_quality=min_mapping_quality, min_base_quality=min_base_quality)
    # print(Normal_genotype_likelihoods)
    break

with open("likelihood.pickle","wb") as fout:
    pickle.dump(Normal_genotype_likelihoods, fout)


with open("likelihood.pickle","rb") as fin:
    loaded_genotype_likelihoods = pickle.load(fin)

print(len(Normal_genotype_likelihoods), len(loaded_genotype_likelihoods))
print(type(Normal_genotype_likelihoods), type(loaded_genotype_likelihoods))
# print(Normal_genotype_likelihoods.values())
for key, value in Normal_genotype_likelihoods.items():
    print(value == loaded_genotype_likelihoods[key])
