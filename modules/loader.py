#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
module for loading the tumor and normal files to run Conpair against
This should support loading from
- filepath
- filepath glob
- directory
- file list .txt (one filepath per line)
- sample manifest ; file containing sample ID's for each filepath
"""
import os
import glob
import json
import itertools

def load_comparisons(
    normals_list = None, # str: text file with file paths, one per line
    tumors_list = None, # str: text file with file paths, one per line
    tumor = None, # str: glob pattern or path to tumor file(s)
    normal = None, # str: glob pattern or path to tumor file(s)
    num_tumors = 'all', # str: string representing int of number of files to load if we only want a subset of the input list for testing
    num_normals = 'all', # str: string representing int of number of files to load if we only want a subset of the input list for testing
    use_manifests = False, # bool: for each file loaded, if an adjacent .json file exists, load it and look for an 'id' field with an alternative sample ID to use
    manifest_dir = None # str: alternative location to look for manifest files
    ): # -> Tuple[ List[str, str, str, str], int, int ]
    """
    Load the samples from file and make the comparisons of tumors vs normals
    """
    # load the paths to pileups
    tumor_pileups = None
    normal_pileups = None

    # evaluate if paths / glob patterns passed
    if tumor:
        tumor_pileups = glob.glob(tumor)
    if normal:
        normal_pileups = glob.glob(normal)

    # try to load from files if nothing has been loaded yet
    if tumor_pileups is None and tumors_list:
        with open(tumors_list) as fin:
            tumor_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    if normal_pileups is None and normals_list:
        with open(normals_list) as fin:
            normal_pileups  = [ line.strip() for line in fin if line.strip() != '' ]

    help_message = "Please run again with positional args for 'tumor' and 'normal', or with args for '--tumors-list' and '--normals-list'."
    if tumor_pileups is None:
        raise Exception("No tumor files loaded. " + help_message)
    if normal_pileups is None:
        raise Exception("No normal files loaded. " + help_message)

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
        tumor_name = get_sample_name(tumor_pileup, use_manifests = use_manifests, manifest_dir = manifest_dir)
        normal_name = get_sample_name(normal_pileup, use_manifests = use_manifests, manifest_dir = manifest_dir)
        labeled_pairs.append((tumor_pileup, normal_pileup, tumor_name, normal_name))
    return(labeled_pairs, num_tumors_loaded, num_normals_loaded)

def get_sample_name(
    filepath, # str: filepath to the input file to evaluate
    use_manifests = True, # bool: for each file loaded, if an adjacent .json file exists, load it and look for an 'id' field with an alternative sample ID to use
    manifest_dir = None # str: alternative location to look for manifest files
    ): # -> str : sample ID for the file
    """
    Method for getting a sample name from a file
    if a valid manifest JSON file is available, return its 'id' field value
    otherwise return default sample name from the file's basename

    a valid manifest file has the same name as the input filepath but with .json appended to the end,
    and contains a single dict with string field 'id'
    """
    if use_manifests:
        # see if a manifest .json file exists
        dirpath = os.path.dirname(filepath)
        filename = os.path.basename(filepath)

        # in case we want to look in another dir for manifests
        if manifest_dir:
            dirpath = manifest_dir
        manifest_filepath = os.path.join(dirpath, filename + ".json")

        # if the file exists, load it and look for 'id' field
        # use try/except here to attempt to minimize the number of operations we execute
        # because this method will potentially get called thousands of times on large datasets
        try:
            with open(manifest_filepath) as fin:
                return(str(json.load(fin)['id']))
        except KeyError: # 'id' not in JSON dict
            return(make_default_sample_name(filepath))
        except TypeError: # JSON contents is not a dict
            return(make_default_sample_name(filepath))
        except IOError: # [Errno 2] No such file or directory
            return(make_default_sample_name(filepath))
    return(make_default_sample_name(filepath))

def make_default_sample_name(filepath):
    # get the sample name from the basename of the file, strip off anything after "."
    return(os.path.basename(filepath).split('.')[0])
