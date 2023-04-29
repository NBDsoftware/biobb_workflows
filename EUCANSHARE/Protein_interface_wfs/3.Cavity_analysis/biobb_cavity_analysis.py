##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import re
import time
import glob
import argparse
import csv
import json
import shutil
from pathlib import Path

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter


def add_suffix(original_path, suffix):
    '''
    Adds suffix to original_path before file extension. For example:

    Inputs
    ------

        original_path  (Path class):  path to original file including filename
        suffix                (str):  suffix string 

    Output
    ------

        new_path   (str): new path

    original_path = Path('/home/user/ClusteringWF/output/step2_extract_models/output.pdb')
    suffix = '0'

    new_path = '/home/user/ClusteringWF/output/step2_extract_models/output_0.pdb'
    '''
    # Identify original file name, extension and path
    filename = original_path.stem
    fileExtension = original_path.suffix
    directory = original_path.parent

    # Create new file name
    new_fileName = filename + '_' + suffix + fileExtension

    # Create new path
    new_path = os.path.join(str(directory), new_fileName)

    return new_path

def add_suffix_to_paths(all_paths, suffix, *keywords):
    '''
    Goes over paths corresponding to keywords in all_paths dictionary, adds 
    suffix to each filename.
    
    Inputs
    ------

        all_paths  (dict): dictionary with old paths of step
        suffix      (str): string corresponding to suffix
        keywords    (str): keywords of output/input paths that will be modified

                    /path/to/file/filename.txt  ---> /path/to/file/filename_suffix.txt

    Output
    ------


        (implicit) all_paths  (dict): dictionary with paths corresponding to "keywords" modified
    '''

    # For all keys passed in keywords, modify path 
    for key in keywords:

        # Original path
        original_path = Path(all_paths[key])

        # Add suffix to path
        newpath = add_suffix(original_path, suffix)

        # Update paths dictionary
        all_paths.update({key : newpath})

    return

def main_wf(configuration_path):
    '''
    Main cavity analysis workflow. This workflow analyzes the cavities of the input structures structures, filters the cavities 
    according to pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------
    
        configuration_path (str): path to input.yml 

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Create a folder for the extracted poses
    fu.create_dir(global_prop["step0_extract_poses"]["path"])

    # STEP 0: extract poses from zip file

    # Extract the poses from the zip file
    poses_paths = fu.unzip_list(global_paths["step0_extract_poses"]["input_zip_path"], global_prop["step0_extract_poses"]["path"])

    # STEP 1: Cavity analysis
    global_log.info("step1_cavity_analysis: Compute cavities using fpocket")

    # Analyze cavities of each pdb 
    for pose_path in poses_paths:

        # Copy and modify the path names according to pdb name
        paths_fpocket = global_paths["step1_cavity_analysis"].copy()
        add_suffix_to_paths(paths_fpocket, Path(pose_path).stem, "output_pockets_zip", "output_summary")

        # Update input pdb path
        paths_fpocket.update({'input_pdb_path': pose_path})

        # Run cavity analysis
        fpocket_run(**paths_fpocket, properties=global_prop["step1_cavity_analysis"])

    # STEP 4: Filtering cavities
    global_log.info("step2_filter_cavities: Filter found cavities")

    # Filter the cavities found for each pdb
    for pose_path in poses_paths:

        # Copy and modify the path names according to pdb name
        paths_filter = global_paths["step2_filter_cavities"].copy()
        add_suffix_to_paths(paths_filter, Path(pose_path).stem, "input_pockets_zip", "input_summary", "output_filter_pockets_zip")

        # Filter cavities
        fpocket_filter(**paths_filter, properties=global_prop["step2_filter_cavities"])

    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)

    args = parser.parse_args()

    main_wf(configuration_path = args.config_path)

