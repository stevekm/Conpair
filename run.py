#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import the Conpair functionality from its scripts and modules and run it in parallel for all samples in the lists
"""
import csv
import os
import sys
import datetime
import itertools
from multiprocessing import Pool
import json
import argparse
import modules.concordance
from modules.ContaminationMarker import get_markers, genotype_likelihoods_for_markers

# need to find a default set of targets to use; Conpair-GRCh37-default
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
default_marker_file = os.path.join(THIS_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

# default parameters to run concordance with
concordance_params = {
'min_mapping_quality': 10,
'normal_homozygous_markers_only': False,
'min_cov': 10,
'min_base_quality': 20
}

def run_conpair(tumor_pileup, normal_pileup, actions_list, concordance_params, markers_data):
    """
    Run the Conpair wrapper scripts on the tumor normal pair and collect the results
    """
    tumor_name = os.path.basename(tumor_pileup).split('.')[0]
    normal_name = os.path.basename(normal_pileup).split('.')[0]
    timestart = datetime.datetime.now()

    # dict to hold results and metrics
    result = {
        'pair': {
            'tumor': tumor_name,
            'normal': normal_name,
            'tumor_pileup': tumor_pileup,
            'normal_pileup': normal_pileup
        }
    }
    if 'concordance' in actions_list:
        try:
            concordance, num_markers_used, num_total_markers = modules.concordance.main(
                tumor_pileup = tumor_pileup,
                normal_pileup = normal_pileup,
                min_mapping_quality = concordance_params['min_mapping_quality'],
                normal_homozygous_markers_only = concordance_params['normal_homozygous_markers_only'],
                markers_data = markers_data,
                min_cov = concordance_params['min_cov'],
                min_base_quality = concordance_params['min_base_quality']
            )
        except ZeroDivisionError:
            print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements; tumor_pileup: {}, normal_pileup: {}'.format(tumor_pileup, normal_pileup))
            concordance = None
            num_markers_used = None
            num_total_markers = None
        result['concordance'] = {}
        result['concordance']['time'] = (datetime.datetime.now() - timestart).seconds
        result['concordance']['concordance'] = concordance
        result['concordance']['num_markers_used'] = num_markers_used
        result['concordance']['num_total_markers'] = num_total_markers

    result['time'] = (datetime.datetime.now() - timestart).seconds
    return(result)

def print_to_console(timestop, num_pairs, time_taken):
    """
    """
    # print to console
    print("[{timestamp}] {num_pairs} pairs processed in {time_taken}s".format(
    timestamp = timestop,
    num_pairs = num_pairs,
    time_taken = time_taken
    ))

def save_JSON_output(output, timestart, time_taken, num_pairs, num_tumors, num_normals, num_threads, actions_list, log_file = None):
    """
    """
    if not log_file:
        log_file = "output.{}pairs.{}t.{}n.{}thread.{}s.{}.json".format(
            num_pairs,
            num_tumors,
            num_normals,
            num_threads,
            time_taken,
            timestart.strftime('%Y-%m-%d_%H-%M-%S')
        )
    log = {
    'output': output,
    'time': time_taken,
    'num_pairs': num_pairs,
    'num_tumors': num_tumors,
    'num_normals': num_normals,
    'threads': num_threads,
    'actions': actions_list
    }
    with open(log_file, "w") as fout:
        json.dump(log, fout, indent = 4)

def save_benchmarks(num_threads, time_taken, num_pairs, num_tumors, num_normals, actions_str, write_mode = "a", output_file = "benchmarks.tsv"):
    """
    """
    with open(output_file, write_mode) as fout:
        line = '\t'.join([str(num_threads), str(time_taken), str(num_pairs), str(num_tumors), str(num_normals), actions_str]) + '\n'
        fout.write(line)

def save_table_output(output):
    """
    """
    has_concordance = any([ 'concordance' in item for item in output ])
    has_contamination = any([ 'contamination' in item for item in output ])

    if has_concordance:
        conc_out = open('concordance.tsv', "w")
        conc_fieldnames = ['concordance', 'num_markers_used', 'num_total_markers', 'tumor', 'normal', 'tumor_pileup', 'normal_pileup']
        conc_writer = csv.DictWriter(conc_out, delimiter = '\t', fieldnames = conc_fieldnames)
        conc_writer.writeheader()
    if has_contamination:
        cont_out = open('concordance.tsv', "w")
        cont_fieldnames = ['tumor', 'normal', 'tumor_pileup', 'normal_pileup']
        cont_writer = csv.DictWriter(cont_out, delimiter = '\t', fieldnames = cont_fieldnames)
        cont_writer.writeheader()

    for item in output:
        if 'concordance' in item:
            row = {
            'tumor': item['pair']['tumor'],
            'tumor_pileup': item['pair']['tumor_pileup'],
            'normal': item['pair']['normal'],
            'normal_pileup': item['pair']['normal_pileup'],
            'concordance': item['concordance']['concordance'],
            'num_total_markers': item['concordance']['num_total_markers'],
            'num_markers_used': item['concordance']['num_markers_used']
            }
            conc_writer.writerow(row)

    if has_concordance:
        conc_out.close()
    if has_contamination:
        cont_out.close()

def main(**kwargs):
    """
    Main control function for the script
    """
    num_threads = kwargs.pop('num_threads', 4)
    num_tumors = kwargs.pop('num_tumors', 1)
    num_normals = kwargs.pop('num_normals', 1)
    actions = kwargs.pop('actions', "concordance,contamination")
    tumor_file = kwargs.pop('tumor_file', "tumors.txt")
    normal_file = kwargs.pop('normal_file', "normals.txt")
    markers = kwargs.pop('markers', default_marker_file)

    actions_list = actions.split(',')

    timestart = datetime.datetime.now()

    # load the data for the markers
    markers_data = get_markers(markers)

    # load the paths to pileups
    with open(tumor_file) as fin:
        all_tumor_pileups  = [ line.strip() for line in fin if line.strip() != '' ]
    with open(normal_file) as fin:
        all_normal_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    tumor_pileups = all_tumor_pileups[0:num_tumors]
    normal_pileups = all_normal_pileups[0:num_normals]

    # get all combinations of both lists of pileups
    pairs = list(itertools.product(tumor_pileups, normal_pileups))

    # run the Copair scripts in parallele in a multi-thread process pool
    pool = Pool(num_threads)
    results = [ pool.apply_async(run_conpair, args=(tumor_pileup, normal_pileup, actions_list, concordance_params, markers_data)) for tumor_pileup, normal_pileup in pairs ]
    output = [ p.get() for p in results ]

    timestop = datetime.datetime.now()
    time_taken = (timestop - timestart).seconds
    num_pairs = len(output)
    num_tumors = len(tumor_pileups)
    num_normals = len(normal_pileups)

    print_to_console(timestop, num_pairs, time_taken)
    save_JSON_output(output, timestart, time_taken, num_pairs, num_tumors, num_normals, num_threads, actions_list)
    save_benchmarks(num_threads, time_taken, num_pairs, num_tumors, num_normals, '.'.join(actions_list))
    save_table_output(output)

def parse():
    """
    Parse the command line options
    """
    parser = argparse.ArgumentParser(description = 'Generate cBio Portal metadata files from various input files')
    parser.add_argument('-t', '--threads', dest = 'num_threads', default = 4, type = int, help = 'The number of CPU threads to use')
    parser.add_argument('--tumors', dest = 'num_tumors', default = 1, type = int, help = 'The number of tumor samples to use from the list')
    parser.add_argument('--normals', dest = 'num_normals', default = 1, type = int,  help = 'The number of normal samples to use from the list')
    parser.add_argument('--actions', dest = 'actions', default = 'concordance,contamination', help = 'A comma separate list of actions to run')
    parser.add_argument('--tumor-file', dest = 'tumor_file', default = "tumors.txt", help = 'File with a list filepaths to the pileups of the tumor samples to use')
    parser.add_argument('--normal-file', dest = 'normal_file', default = "normals.txt", help = 'File with a list filepaths to the pileups of the normal samples to use')
    parser.add_argument('--markers', dest = 'markers', default = default_marker_file, help = 'Markers to use for analysis')

    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    parse()
