#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for the concordance module
"""
import os
import unittest
from loader import load_comparisons

# need to get the path to some sample files
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
PILEUP_DIR = os.path.join(PARENT_DIR, "data", "example", "pileup")

class TestLoadComparisons(unittest.TestCase):
    maxDiff = None
    def test_load_single_files1(self):
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt")
        normal_file = os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt")
        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normal = normal_file)
        expected_pairs = [(tumor_file,
            normal_file,
            "NA12878_tumor80x",
            "NA12878_normal40x"
            )]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 1)

    def test_load_glob1(self):
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt")
        normal_file = os.path.join(PILEUP_DIR, "*normal*.pileup*.txt")
        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normal = normal_file)
        expected_pairs = [
            (tumor_file, os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt"), "NA12878_tumor80x", "NA12878_normal40x"),
            (tumor_file, os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt"), "NA12878_tumor80x", "NA12878_normal40x"),
            ]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 2)

    def test_load_glob2(self):
        tumor_file = os.path.join(PILEUP_DIR, "*tumor*.pileup.txt")
        normal_file = os.path.join(PILEUP_DIR, "*normal*.pileup*.txt")
        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normal = normal_file)
        expected_pairs = [
            (
                os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt"),
                os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt"),
                "NA12878_tumor80x",
                "NA12878_normal40x"),
            (
                os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt"),
                os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt"),
                "NA12878_tumor80x",
                "NA12878_normal40x"),
            (
                os.path.join(PILEUP_DIR, "NA12878-2_tumor80x.gatk.pileup.txt"),
                os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt"),
                "NA12878-2_tumor80x",
                "NA12878_normal40x"),
            (
                os.path.join(PILEUP_DIR, "NA12878-2_tumor80x.gatk.pileup.txt"),
                os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt"),
                "NA12878-2_tumor80x",
                "NA12878_normal40x"),
            ]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 2)
        self.assertEqual(num_normals_loaded, 2)




if __name__ == "__main__":
    unittest.main()

