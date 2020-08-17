#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import csv
import datetime
import itertools
import argparse
from multiprocessing import Pool
from modules.ContaminationMarker import get_markers
from modules.concordance import concordance

# get the path to the included default margers; Conpair-GRCh37-default
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
default_marker_file = os.path.join(THIS_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

timestart = datetime.datetime.now()

def load_comparisons(normals_list, tumors_list, num_tumors = 'all', num_normals = 'all'):
    """
    Load the samples from file and make the comparisons of tumors vs normals
    """
    # load the paths to pileups
    with open(tumors_list) as fin:
        tumor_pileups  = [ line.strip() for line in fin if line.strip() != '' ]
    with open(normals_list) as fin:
        normal_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    # subset the lists if desired
    if num_tumors != 'all':
        tumor_pileups = [ t for t in tumor_pileups[0:int(num_tumors)] ]
    if num_normals != 'all':
        normal_pileups = [ t for t in normal_pileups[0:int(num_normals)] ]

    num_tumors_loaded = len(tumor_pileups)
    num_normals_loaded = len(normal_pileups)

    # get all combinations of both lists of pileups
    pairs = list(itertools.product(tumor_pileups, normal_pileups))

    # add sample ID labels to each pair
    labeled_pairs = []
    for pair in pairs:
        tumor_pileup, normal_pileup = pair
        tumor_name = os.path.basename(tumor_pileup).split('.')[0]
        normal_name = os.path.basename(normal_pileup).split('.')[0]
        labeled_pairs.append((tumor_pileup, normal_pileup, tumor_name, normal_name))
    return(labeled_pairs, num_tumors_loaded, num_normals_loaded)


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
            concordance_val = None
            num_markers_used = None
            num_total_markers = None
        yield(tumor_pileup, normal_pileup, tumor_name, normal_name, concordance_val, num_markers_used, num_total_markers)

def save_benchmarks_to_file(benchmarks_file, num_threads, num_pairs, num_tumors, num_normals, action):
    """
    """
    timestop = datetime.datetime.now()
    time_taken = (timestop - timestart).seconds
    with open(benchmarks_file, 'a') as fout:
        line = '\t'.join([str(num_threads), str(time_taken), str(num_pairs), str(num_tumors), str(num_normals), str(action)]) + '\n'
        fout.write(line)

def run_concordance(**kwargs):
    """
    """
    output_file = kwargs.pop('output_file', 'concordance.tsv')
    num_normals = kwargs.pop('num_normals', 'all')
    num_tumors = kwargs.pop('num_tumors', 'all')
    normals_list = kwargs.pop('normals_list', "normals.txt")
    tumors_list = kwargs.pop('tumors_list', "tumors.txt")
    markers = kwargs.pop('markers', default_marker_file)
    num_threads = kwargs.pop('num_threads', 4)
    min_mapping_quality = kwargs.pop('min_mapping_quality', 10)
    normal_homozygous_markers_only = kwargs.pop('normal_homozygous_markers_only', False)
    min_cov = kwargs.pop('min_cov', 10)
    min_base_quality = kwargs.pop('min_base_quality', 20)
    save_benchmarks = kwargs.pop('save_benchmarks', False)
    benchmarks_file = kwargs.pop('benchmarks_file', 'benchmarks.tsv')

    # load all comparisons
    pairs, num_tumors_loaded, num_normals_loaded = load_comparisons(normals_list, tumors_list, num_tumors, num_normals)
    num_pairs = len(pairs)

    # load the data for the markers
    markers_data = get_markers(markers)

    # open output file for writing
    if output_file == '-':
        fout = sys.stdout
    else:
        fout = open(output_file, "w")

    conc_fieldnames = ['concordance', 'num_markers_used', 'num_total_markers', 'tumor', 'normal', 'tumor_pileup', 'normal_pileup']
    conc_writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = conc_fieldnames, lineterminator='\n')
    conc_writer.writeheader()

    for tumor_pileup, normal_pileup, tumor_name, normal_name, concordance_val, num_markers_used, num_total_markers in run_parallel_concordance(pairs, markers_data, num_threads, min_mapping_quality, normal_homozygous_markers_only, min_cov, min_base_quality):
        conc_writer.writerow({
        'tumor_pileup': tumor_pileup,
        'normal_pileup': normal_pileup,
        'tumor': tumor_name,
        'normal': normal_name,
        'concordance': concordance_val,
        'num_markers_used': num_markers_used,
        'num_total_markers': num_total_markers
        })
    fout.close()

    if save_benchmarks:
        save_benchmarks_to_file(benchmarks_file = benchmarks_file, num_threads = num_threads, num_pairs = num_pairs, num_tumors = num_tumors, num_normals = num_normals, action = "concordance")




def parse():
    """
    Command line argument parsing when run as a script
    """
    parser = argparse.ArgumentParser(description = 'Run Conpair wrapper script')
    subparsers = parser.add_subparsers(help ='Sub-commands available')

    concordance_parser = subparsers.add_parser('concordance', help = 'Run Conpair concordance')
    concordance_parser.add_argument('--tumors-list', dest = 'tumors_list', default = "tumors.txt", help = 'File with a list filepaths to the pileups of the tumor samples to use')
    concordance_parser.add_argument('--normals-list', dest = 'normals_list', default = "normals.txt", help = 'File with a list filepaths to the pileups of the normal samples to use')
    concordance_parser.add_argument('--num-tumors', dest = 'num_tumors', default = 'all', help = 'The number of tumor samples to use from the list')
    concordance_parser.add_argument('--num-normals', dest = 'num_normals', default = 'all', help = 'The number of normal samples to use from the list')
    concordance_parser.add_argument('--markers', dest = 'markers', default = default_marker_file, help = 'Markers to use for analysis')
    concordance_parser.add_argument('-t', '--threads', dest = 'num_threads', default = 4, type = int, help = 'The number of CPU threads to use')
    concordance_parser.add_argument('--min-mapping-quality', dest = 'min_mapping_quality', default = 10, type = int, help = 'Minimum mapping quality value')
    concordance_parser.add_argument('--min-cov', dest = 'min_cov', default = 10, type = int, help = 'Minimum coverage quality value')
    concordance_parser.add_argument('--min_base_quality', dest = 'min_base_quality', default = 20, type = int, help = 'Minimum base quality value')
    concordance_parser.add_argument('--normal-homozygous-markers-only', dest = 'normal_homozygous_markers_only', default = False, action = "store_true", help = 'Use only homozygous markers')
    concordance_parser.add_argument('--output-file', dest = 'output_file', default = "concordance.tsv", help = 'File to output concordance data to')
    concordance_parser.add_argument('--save-benchmarks', dest = 'save_benchmarks', action='store_true', help = 'Append benchmarks to a file')
    concordance_parser.add_argument('--benchmarks-file', dest = 'benchmarks_file', default='benchmarks.tsv', help = 'File to append benchmarks to')

    concordance_parser.set_defaults(func = run_concordance)

    args = parser.parse_args()
    args.func(**vars(args))


if __name__ == '__main__':
    parse()
