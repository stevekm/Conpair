#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for the concordance module
"""
import os
import unittest
from concordance import concordance
from ContaminationMarker import get_markers
import tempfile
import shutil

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
PILEUP_DIR = os.path.join(PARENT_DIR, "data", "example", "pileup")
marker_file = os.path.join(PARENT_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

class TestCocordance(unittest.TestCase):
    def test_concordance1(self):
        """
        Test the concordance of the included pileup tumor and normal files
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        normal_pileup = os.path.join(PILEUP_DIR, 'NA12878_normal40x.gatk.pileup.txt')
        markers_data = get_markers(marker_file)

        concordance_val, num_markers_used, num_total_markers = concordance(tumor_pileup = tumor_pileup, normal_pileup = normal_pileup, markers_data = markers_data)

        self.assertEqual(concordance_val, 0.9993209289691701)
        self.assertEqual(num_markers_used, 7363)
        self.assertEqual(num_total_markers, 7387)

    def test_concordance2(self):
        """
        Test that the included tumor pileup has 100% concordance to itself
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        markers_data = get_markers(marker_file)

        concordance_val, num_markers_used, num_total_markers = concordance(tumor_pileup = tumor_pileup, normal_pileup = tumor_pileup, markers_data = markers_data)

        self.assertEqual(concordance_val, 1.0)
        self.assertEqual(num_markers_used, 7374)
        self.assertEqual(num_total_markers, 7387)

    def test_concordance_subset(self):
        """
        Check the concordance values of the included example tumor and normal for various numbers of subset entries from the pileup
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        normal_pileup = os.path.join(PILEUP_DIR, 'NA12878_normal40x.gatk.pileup.txt')
        markers_data = get_markers(marker_file)

        tmpdirpath = tempfile.mkdtemp()
        subset_tumor_5000_file = os.path.join(tmpdirpath, "tumor5000.pileup")
        subset_normal_5000_file = os.path.join(tmpdirpath, "normal5000.pileup")
        with open(tumor_pileup) as fin, open(subset_tumor_5000_file, 'w') as fout:
            for i, line in enumerate(fin):
                if i < 5000:
                    fout.write(line)
                else:
                    break
        with open(normal_pileup) as fin, open(subset_normal_5000_file, 'w') as fout:
            for i, line in enumerate(fin):
                if i < 5000:
                    fout.write(line)
                else:
                    break

        concordance_val, num_markers_used, num_total_markers = concordance(tumor_pileup = subset_tumor_5000_file, normal_pileup = subset_normal_5000_file, markers_data = markers_data)

        self.assertEqual(concordance_val, 0.9993995196156925)
        self.assertEqual(num_markers_used, 4996)
        self.assertEqual(num_total_markers, 7387)

        shutil.rmtree(tmpdirpath)

    def test_load_normal_from_pickle(self):
        """
        Test that the Normal pileup likelihoods can be loaded from .pickle file
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        normal_pileup = os.path.join(PILEUP_DIR, 'NA12878_normal40x.gatk.pileup.txt.pickle')
        markers_data = get_markers(marker_file)

        concordance_val, num_markers_used, num_total_markers = concordance(tumor_pileup = tumor_pileup, normal_pileup = normal_pileup, markers_data = markers_data)

        self.assertEqual(concordance_val, 0.9993209289691701)
        self.assertEqual(num_markers_used, 7363)
        self.assertEqual(num_total_markers, 7387)






if __name__ == "__main__":
    unittest.main()
