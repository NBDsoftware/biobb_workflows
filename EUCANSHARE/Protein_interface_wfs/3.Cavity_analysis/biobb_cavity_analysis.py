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

def create_summary(pockets_path, default_summary_path, global_log):
    '''
    Print in log file all available pockets for each model found in the input folder for the pocket selection step. 
    Include also information for each pocket.
    
    Inputs
    ------

        pockets_path         (str): Path to pockets folder 
        default_summary_path (str): Default path to cavity summary with information for each pocket
        global_log        (Logger): Object of class Logger (dumps info to global log file of this workflow)
    
    Output
    ------

        Info dumped to log file:

        global_summary (dict):  dictionary with models and pockets ordered by modelname. See example:
        
            {
                modelname1 : {         
                    pockets : [1 , 2], 
                    pocket1 : {info pocket 1}, 
                    pocket2 : {info pocket 2}
                    }, 
                modelname2 : {
                    ...
                } 
            }
    '''
    
    global_log.info("    Available models and pockets after filtering:")

    # Pattern that can be used to extract ID from string with pocket name: (str) pocket3 -> (int) 3
    pocketID_pattern = r'pocket(\d+)'

    # Pattern that can be used to extract model ID from string with model name
    # This ID reflects the ordering of the population, not the original ID of the clustering step
    modelID_pattern = r'cluster_(\d+)'

    # Dictionary where all available model_summary dictionaries will be saved
    global_summary = {}

    # Convert string to Path class
    pockets_path = Path(pockets_path)

    # Pattern for pocket filtering log file names
    logNamePattern = pockets_path.name + "_log*.out"

    # Find all log files matching pattern
    logList = sorted(pockets_path.rglob(logNamePattern))

    # Iterate through all log files -> find model name and available pockets if any
    for log in logList:

        # Find model name NOTE: hardcoded name of step and file, should be changed
        modelName = find_matching_line(pattern=r'step1_cavity_analysis/all_pockets_(\S+).zip', filepath=str(log))

        # Find pockets  NOTE: is there a more robust way of finding models and pockets after filtering? 
        pocket_lines = find_matching_lines(pattern=r'pocket\d+$', filepath=str(log))

        # If some model and pocket/s are found in log
        if None not in (modelName, pocket_lines):

            # Dictionary where pockets and information of each pocket will be saved
            model_summary = {}

            # Path to this model's summary with information for all pocket found
            summary_path = add_suffix(original_path=Path(default_summary_path), suffix=modelName)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # list with pocket IDs for this model
            pocket_IDs = []

            # For each filtered pocket: pocket1, pocket4 ...
            for pocket_line in pocket_lines:

                # Strip whitespace if any
                pocket_line = pocket_line.strip()

                # Save entry with pocket information
                model_summary.update({pocket_line : all_pockets_summary[pocket_line]})

                # Search for the pocket ID in pocket name
                match = re.search(pocketID_pattern, pocket_line)

                # Append pocket ID to list
                pocket_IDs.append(int(match[1]))

            # Save pocket IDs
            model_summary.update({'pockets' : pocket_IDs})

            # Update global_summary
            global_summary.update({modelName : model_summary})

    # Print in log file with readable format
    pretty_summary = json.dumps(global_summary, indent=2)

    global_log.info(pretty_summary)

    return 

def find_matching_line(pattern, filepath):
    '''
    Finds first line in file containing a given pattern

    Inputs
    ------

        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------
    
        line (str): line matching the pattern or None if no line matches the pattern
    '''

    # Open log file
    file = open(filepath, 'r')
    
    # Read all lines
    lines = file.readlines()

    # Search matching string in each line
    for line in lines:

        match = re.search(pattern, line)

        if match is not None:

            # Close file
            file.close()

            return match[1]
    
    file.close()

    return None

def find_matching_lines(pattern, filepath):
    '''
    Finds all lines in file containing a given pattern

    Inputs
    ------

        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------

        lines (list(str)): lines matching the pattern or None if no line matches the pattern
    '''

    # Open log file
    file = open(filepath, 'r')
    
    # Read all lines
    lines = file.readlines()

    # List to save matching lines
    matchingLines = []

    # Search matching string in each line
    for line in lines:

        match = re.search(pattern, line)

        if match is not None:

            matchingLines.append(match[0])

    # Close file
    file.close()

    # If no lines contained the pattern, return None
    if len(matchingLines) == 0:
        matchingLines = None
    
    return matchingLines

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

    # Create and print available models with their pockets and populations after filtering
    create_summary(pockets_path = global_prop["step2_filter_cavities"]["path"], 
                   default_summary_path = global_paths['step1_cavity_analysis']['output_summary'],
                   global_log = global_log)
    
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

