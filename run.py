#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the Conpair scripts in parallel for all combinations of tumor and normal pileups
"""
import sys
import os
import subprocess as sp
import itertools
from multiprocessing import Pool
import datetime
import re
import json

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
SCRIPT_DIR = os.path.join(THIS_DIR, "scripts")
verify_concordance_script = os.path.join(SCRIPT_DIR, 'verify_concordance.py')
estimate_tumor_normal_contamination_script = os.path.join(SCRIPT_DIR, 'estimate_tumor_normal_contamination.py')

def run_command(args):
    """
    Helper function to run a shell command easier
    Parameters
    ----------
    args: list
        a list of shell args to execute
    """
    process = sp.Popen(args, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True)
    proc_stdout, proc_stderr = process.communicate()
    returncode = process.returncode
    proc_stdout = proc_stdout.strip()
    proc_stderr = proc_stderr.strip()
    return(returncode, proc_stdout, proc_stderr)

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
        command = [ 'python2.7',  verify_concordance_script, '--tumor_pileup', tumor_pileup, '--normal_pileup', normal_pileup ]
        returncode, proc_stdout, proc_stderr = run_command(command)

        result['concordance'] = {}
        result['concordance']['returncode'] = returncode
        result['concordance']['stdout'] = proc_stdout
        result['concordance']['stderr'] = proc_stderr
        result['concordance']['time'] = str((datetime.datetime.now() - timestart).seconds)

        concordance = None
        if len(proc_stdout.split()) > 0:
            try:
                concordance = float(proc_stdout.split()[0].strip())
            except:
                pass
        result['concordance']['concordance'] = concordance

    if 'contamination' in actions_list:
        command = [ 'python2.7',  estimate_tumor_normal_contamination_script, '--tumor_pileup', tumor_pileup, '--normal_pileup', normal_pileup ]
        returncode, proc_stdout, proc_stderr = run_command(command)

        result['contamination'] = {}
        result['contamination']['returncode'] = returncode
        result['contamination']['stdout'] = proc_stdout
        result['contamination']['stderr'] = proc_stderr
        result['contamination']['time'] = str((datetime.datetime.now() - timestart).seconds)

        normal_contamination = None
        tumor_contamination = None
        matches = re.findall(r'[0-9]+.[0-9]+', proc_stdout)
        if len(matches) > 0:
            try:
                normal_contamination = float(matches[0])
                tumor_contamination = float(matches[1])
            except:
                pass

        result['contamination']['normal_contamination'] = normal_contamination
        result['contamination']['tumor_contamination'] = tumor_contamination

    result['time'] = str((datetime.datetime.now() - timestart).seconds)

    return(result)

def main():
    """
    run.py <num_threads> <num_tumors> <num_normals> concordance,contamination
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

    # get the first one for testing
    # run the Copair scripts in parallele in a multi-thread process pool
    pool = Pool(num_threads)
    results = [ pool.apply_async(run_conpair, args=(tumor_pileup, normal_pileup, actions_list)) for tumor_pileup, normal_pileup in pairs ] # [0:num_pairs]
    output = [ p.get() for p in results ]

    timestop = datetime.datetime.now()
    time_taken = str((timestop - timestart).seconds)
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

    # save benchmarks
    with open("benchmarks.tsv", "a") as fout:
        line = '\t'.join([str(num_threads), str(time_taken), str(num_pairs), str(num_tumors), str(num_normals), actions_str]) + '\n'
        fout.write(line)




if __name__ == '__main__':
    main()
