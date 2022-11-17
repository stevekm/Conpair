#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for the concordance module
"""
import os
import json
import unittest
import shutil
from tempfile import mkdtemp
from loader import load_comparisons

# need to get the path to some sample files
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
PILEUP_DIR = os.path.join(PARENT_DIR, "data", "example", "pileup")

class TestLoadComparisons(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.tmpdir = mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_single_files1(self):
        """
        Load single files from their paths
        """
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
        """
        Load multiple files from glob patterns
        """
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
        """
        Load even more files from their glob patterns
        """
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

    def test_load_filelist1(self):
        """
        Load files from a .txt list of paths
        """
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt")
        normal_list_file = os.path.join(self.tmpdir, "normals.txt")
        with open(normal_list_file, "w") as fout:
            fout.write(os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt") + "\n")

        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normals_list = normal_list_file)
        expected_pairs = [
            (tumor_file, os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt"), "NA12878_tumor80x", "NA12878_normal40x"),
            ]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 1)

    def test_load_filelist2(self):
        """
        Load even more files from a .txt list of paths
        """
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt")
        normal_list_file = os.path.join(self.tmpdir, "normals.txt")
        with open(normal_list_file, "w") as fout:
            fout.write(os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt") + "\n")
            fout.write(os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt") + "\n")

        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normals_list = normal_list_file)
        expected_pairs = [
            (tumor_file, os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.10lines.txt"), "NA12878_tumor80x", "NA12878_normal40x"),
            (tumor_file, os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt"), "NA12878_tumor80x", "NA12878_normal40x"),
            ]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 2)

    def test_load_manifests1(self):
        """
        Load files from path but discover adjacent .json manifest files to load sample ID's from
        """
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt") # NA12878_tumor80x.gatk.pileup.txt.json
        normal_file = os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt")
        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normal = normal_file, use_manifests = True)
        expected_pairs = [(tumor_file,
            normal_file,
            "NA12878",
            "NA12878_normal40x"
            )]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 1)

    def test_load_manifest_dir1(self):
        """
        Load sample ID's from manifests in a different location
        """
        tumor_file = os.path.join(PILEUP_DIR, "NA12878_tumor80x.gatk.pileup.txt")
        normal_file = os.path.join(PILEUP_DIR, "NA12878_normal40x.gatk.pileup.txt")

        tumor_manifest_filepath = os.path.join(self.tmpdir, "NA12878_tumor80x.gatk.pileup.txt.json")
        data = {'id': "NA12878_foo"}
        with open(tumor_manifest_filepath, "w") as fout:
            json.dump(data, fout)

        labeled_pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(tumor = tumor_file, normal = normal_file, use_manifests = True, manifest_dir = self.tmpdir)
        expected_pairs = [(tumor_file,
            normal_file,
            "NA12878_foo",
            "NA12878_normal40x"
            )]
        self.assertEqual(labeled_pairs, expected_pairs)
        self.assertEqual(num_tumors_loaded, 1)
        self.assertEqual(num_normals_loaded, 1)









if __name__ == "__main__":
    unittest.main()

