#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for the concordance module
"""
import os
import unittest
import concordance
from ContaminationMarker import get_markers

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
PILEUP_DIR = os.path.join(PARENT_DIR, "data", "example", "pileup")
marker_file = os.path.join(PARENT_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

class TestCocordance(unittest.TestCase):
    def test_concordance1(self):
        """
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        normal_pileup = os.path.join(PILEUP_DIR, 'NA12878_normal40x.gatk.pileup.txt')
        markers_data = get_markers(marker_file)

        concordance_val, num_markers_used, num_total_markers = concordance.main(tumor_pileup = tumor_pileup, normal_pileup = normal_pileup, markers_data = markers_data)

        self.assertEqual(concordance_val, 0.9993209289691701)
        self.assertEqual(num_markers_used, 7363)
        self.assertEqual(num_total_markers, 7387)

    def test_concordanc2(self):
        """
        """
        tumor_pileup = os.path.join(PILEUP_DIR, 'NA12878_tumor80x.gatk.pileup.txt')
        markers_data = get_markers(marker_file)

        concordance_val, num_markers_used, num_total_markers = concordance.main(tumor_pileup = tumor_pileup, normal_pileup = tumor_pileup, markers_data = markers_data)

        self.assertEqual(concordance_val, 1.0)
        self.assertEqual(num_markers_used, 7374)
        self.assertEqual(num_total_markers, 7387)

if __name__ == "__main__":
    unittest.main()
