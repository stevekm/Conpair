#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import the Conpair functionality from its scripts and modules and run it in parallel for all samples in the lists
"""
import os
import sys
import datetime
import itertools
from multiprocessing import Pool
import json

# import the script as a module from the other directory
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
SCRIPT_DIR = os.path.join(THIS_DIR, "scripts")
sys.path.insert(0, SCRIPT_DIR)
import verify_concordance2
sys.path.pop(0)

def run_conpair(tumor_pileup, normal_pileup, actions_list):
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
        concordance = verify_concordance2.main({'tumor_pileup':tumor_pileup, 'normal_pileup':normal_pileup})
        result['concordance'] = {}
        result['concordance']['time'] = (datetime.datetime.now() - timestart).seconds
        result['concordance']['concordance'] = concordance

    result['time'] = (datetime.datetime.now() - timestart).seconds
    return(result)

def main():
    """
    run2.py <num_threads> <num_tumors> <num_normals> concordance,contamination
    """
    # check if cli args were passed
    args = sys.argv[1:]
    if len(args) > 0:
        num_threads = int(args[0])
        num_tumors = int(args[1])
        num_normals = int(args[2])
        actions = args[3] # concordance,contamination
    else:
        num_threads = 8
        num_tumors = 1
        num_normals = 1
        actions = "concordance,contamination"

    actions_list = actions.split(',')

    timestart = datetime.datetime.now()

    # load the paths to pileups
    tumor_pileup_file = "tumor_pileups.txt"
    normal_pileup_file = "normal_pileups.txt"
    with open(tumor_pileup_file) as fin:
        all_tumor_pileups  = [ line.strip() for line in fin if line.strip() != '' ]
    with open(normal_pileup_file) as fin:
        all_normal_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    tumor_pileups = all_tumor_pileups[0:num_tumors]
    normal_pileups = all_normal_pileups[0:num_normals]

    # get all combinations of both lists of pileups
    pairs = list(itertools.product(tumor_pileups, normal_pileups))

    # run the Copair scripts in parallele in a multi-thread process pool
    pool = Pool(num_threads)
    results = [ pool.apply_async(run_conpair, args=(tumor_pileup, normal_pileup, actions_list)) for tumor_pileup, normal_pileup in pairs ]
    output = [ p.get() for p in results ]

    timestop = datetime.datetime.now()
    time_taken = (timestop - timestart).seconds
    num_pairs = len(output)
    num_tumors = len(tumor_pileups)
    num_normals = len(normal_pileups)
    actions_str = '.'.join(actions_list)

    # print to console
    print("[{timestamp}] {num_pairs} pairs processed in {time_taken}s".format(
    timestamp = timestop,
    num_pairs = num_pairs,
    time_taken = time_taken
    ))

    # save the output log
    log_file = "output2.{}pairs.{}t.{}n.{}thread.{}s.{}.json".format(
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

    # save benchmarks
    with open("benchmarks2.tsv", "a") as fout:
        line = '\t'.join([str(num_threads), str(time_taken), str(num_pairs), str(num_tumors), str(num_normals), actions_str]) + '\n'
        fout.write(line)


if __name__ == '__main__':
    main()
