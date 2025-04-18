#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import csv
import glob
import datetime
import itertools
import argparse
from multiprocessing import Pool
from modules.ContaminationMarker import get_markers
from modules.concordance import concordance
from modules.loader import load_comparisons

# get the path to the included default margers; Conpair-GRCh37-default
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
default_marker_file = os.path.join(THIS_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

timestart = datetime.datetime.now()

def run_parallel_concordance(
    pairs,
    markers_data,
    num_threads,
    min_mapping_quality,
    normal_homozygous_markers_only,
    min_cov,
    min_base_quality):
    """
    Run all the parallel instances of concordance comparisons and yield the results
    """
    # start multiprocessing pool
    pool = Pool(int(num_threads))

    # generate the aysnc result instances
    results = []
    for tumor_pileup, normal_pileup, tumor_name, normal_name in pairs:
        result = pool.apply_async(concordance, args=(
            tumor_pileup,
            normal_pileup,
            markers_data,
            min_mapping_quality,
            normal_homozygous_markers_only,
            min_cov,
            min_base_quality
            )
        )
        result_tup = (tumor_pileup, normal_pileup, tumor_name, normal_name, result)
        results.append(result_tup)

    # get each result
    for result_tup in results:
        tumor_pileup, normal_pileup, tumor_name, normal_name, result = result_tup
        try:
            concordance_val, num_markers_used, num_total_markers = result.get()
        except ZeroDivisionError:
            # if concordant+discordant == 0:
            #     print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(min_cov))
            #     sys.exit(0)
            concordance_val = None
            num_markers_used = None
            num_total_markers = None
        yield(tumor_pileup, normal_pileup, tumor_name, normal_name, concordance_val, num_markers_used, num_total_markers)

def save_benchmarks_to_file(benchmarks_file, num_threads, num_pairs, num_tumors, num_normals, action):
    """
    Append benchmark metrics to a file
    """
    timestop = datetime.datetime.now()
    time_taken = (timestop - timestart).seconds
    with open(benchmarks_file, 'a') as fout:
        line = '\t'.join([str(num_threads), str(time_taken), str(num_pairs), str(num_tumors), str(num_normals), str(action)]) + '\n'
        fout.write(line)

def run_concordance(**kwargs):
    """
    Main control function for running concordance in parallel for a list of tumors and normals
    """
    tumor = kwargs.pop('tumor', None)
    normal = kwargs.pop('normal', None)
    output_file = kwargs.pop('output_file', None) # 'concordance.tsv'
    num_normals = kwargs.pop('num_normals', 'all')
    num_tumors = kwargs.pop('num_tumors', 'all')
    normals_list = kwargs.pop('normals_list', None) # "normals.txt"
    tumors_list = kwargs.pop('tumors_list', None) # "tumors.txt"
    markers = kwargs.pop('markers', default_marker_file)
    num_threads = kwargs.pop('num_threads', 4)
    min_mapping_quality = kwargs.pop('min_mapping_quality', 10)
    normal_homozygous_markers_only = kwargs.pop('normal_homozygous_markers_only', False)
    min_cov = kwargs.pop('min_cov', 10)
    min_base_quality = kwargs.pop('min_base_quality', 20)
    save_benchmarks = kwargs.pop('save_benchmarks', False)
    benchmarks_file = kwargs.pop('benchmarks_file', 'benchmarks.tsv')
    print_filepath = kwargs.pop('print_filepath')
    use_manifests = kwargs.pop('use_manifests', False)
    manifest_dir = kwargs.pop('manifest_dir', None)


    # load all comparisons of each tumor vs each normal
    pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(
        tumor = tumor,
        normal = normal,
        normals_list = normals_list,
        tumors_list = tumors_list,
        num_tumors = num_tumors,
        num_normals = num_normals,
        use_manifests = use_manifests,
        manifest_dir = manifest_dir
        )

    num_pairs = len(pairs)

    # load the data for the markers
    markers_data = get_markers(markers)

    # open output file for writing
    if output_file is None or output_file == '-':
        fout = sys.stdout
    else:
        fout = open(output_file, "w")

    # initialize the output file writer
    conc_fieldnames = ['concordance', 'num_markers_used', 'num_total_markers', 'tumor', 'normal', 'tumor_filename', 'normal_filename']
    if print_filepath:
        conc_fieldnames.append("tumor_filepath")
        conc_fieldnames.append("normal_filepath")
    conc_writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = conc_fieldnames, lineterminator='\n')
    conc_writer.writeheader()

    # run all the comparisons in parallel and write their concordance outputs as they arrive
    for tumor_pileup, normal_pileup, tumor_name, normal_name, concordance_val, num_markers_used, num_total_markers in run_parallel_concordance(pairs, markers_data, num_threads, min_mapping_quality, normal_homozygous_markers_only, min_cov, min_base_quality):
        row = {
        'tumor_filename': os.path.basename(tumor_pileup),
        'normal_filename': os.path.basename(normal_pileup),
        'tumor': tumor_name,
        'normal': normal_name,
        'concordance': concordance_val,
        'num_markers_used': num_markers_used,
        'num_total_markers': num_total_markers
        }
        if print_filepath:
            row["tumor_filepath"] = tumor_pileup
            row["normal_filepath"] = normal_pileup
        conc_writer.writerow(row)
    fout.close()

    if save_benchmarks:
        save_benchmarks_to_file(benchmarks_file = benchmarks_file, num_threads = num_threads, num_pairs = num_pairs, num_tumors = num_tumors_loaded, num_normals = num_tumors_loaded, action = "concordance")

def parse():
    """
    Command line argument parsing when run as a script
    """
    parser = argparse.ArgumentParser(description = 'Run Conpair wrapper script')
    subparsers = parser.add_subparsers(help ='Sub-commands available')

    concordance_parser = subparsers.add_parser('concordance', help = 'Run Conpair concordance')

    concordance_parser.add_argument('tumor', nargs = '?', help = "File path or glob pattern for tumor pileup or likelihoods file")
    concordance_parser.add_argument('normal', nargs = '?', help = "File path or glob pattern for normal pileup or likelihoods file")

    concordance_parser.add_argument('--tumors-list', dest = 'tumors_list', help = 'File with a list filepaths to the pileups of the tumor samples to use') # , default = "tumors.txt"
    concordance_parser.add_argument('--normals-list', dest = 'normals_list', help = 'File with a list filepaths to the pileups of the normal samples to use') # , default = "normals.txt"

    concordance_parser.add_argument('--num-tumors', dest = 'num_tumors', default = 'all', help = 'The number of tumor samples to use from the list')
    concordance_parser.add_argument('--num-normals', dest = 'num_normals', default = 'all', help = 'The number of normal samples to use from the list')
    concordance_parser.add_argument('--markers', dest = 'markers', default = default_marker_file, help = 'Markers to use for analysis')
    concordance_parser.add_argument('-t', '--threads', dest = 'num_threads', default = 4, type = int, help = 'The number of CPU threads to use')
    concordance_parser.add_argument('--min-mapping-quality', dest = 'min_mapping_quality', default = 10, type = int, help = 'Minimum mapping quality value')
    concordance_parser.add_argument('--min-cov', dest = 'min_cov', default = 10, type = int, help = 'Minimum coverage quality value')
    concordance_parser.add_argument('--min_base_quality', dest = 'min_base_quality', default = 20, type = int, help = 'Minimum base quality value')
    concordance_parser.add_argument('--normal-homozygous-markers-only', dest = 'normal_homozygous_markers_only', default = False, action = "store_true", help = 'Use only homozygous markers')
    concordance_parser.add_argument('--output-file', dest = 'output_file', help = 'File to output concordance data to. Use "-" for stdout') # , default = "concordance.tsv"
    concordance_parser.add_argument('--save-benchmarks', dest = 'save_benchmarks', action='store_true', help = 'Append benchmarks to a file')
    concordance_parser.add_argument('--benchmarks-file', dest = 'benchmarks_file', default='benchmarks.tsv', help = 'File to append benchmarks to')
    concordance_parser.add_argument('--filepath', dest = 'print_filepath', action = "store_true", help = "Print the file path in the output")
    concordance_parser.add_argument('--manifests', dest = 'use_manifests', action = "store_true", help = "Load sample IDs from adjacent .json manifest files for each input file")
    concordance_parser.add_argument('--manifest-dir', dest = 'manifest_dir', default = None, help = "Alternate directory to load manifest files from")

    concordance_parser.set_defaults(func = run_concordance)

    args = parser.parse_args()
    args.func(**vars(args))


if __name__ == '__main__':
    parse()
