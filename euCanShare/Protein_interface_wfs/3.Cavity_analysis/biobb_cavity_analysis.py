##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import re
import time
import argparse
import json
import yaml
from pathlib import Path

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_structure_utils.utils.extract_molecule import extract_molecule

from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def create_summary(poses_name_list, global_paths, output_path, output_summary_path = None):
    '''
    Print in log file all available pockets after filtering for each centroid. 
    Include also information for each pocket and population for each centroid if available.
    
    Inputs
    ------

        poses_name_list         (list): List with names of centroids
        global_paths            (dict): global paths dictionary
        output_path              (str): path to output folder
        output_summary_path      (str): path to output summary file
    
    Output
    ------

        Info dumped to yaml summary file
    '''

    # Pattern that can be used to extract pocket ID from string with pocket name: (str) pocket3 -> (int) 3
    pocketID_pattern = r'pocket(\d+)'

    # Find step names
    filter_cavities_folder = 'step3_filter_cavities'
    cavity_analysis_folder = 'step2_cavity_analysis'

    # Find summary file name
    pockets_summary_filename = Path(global_paths[cavity_analysis_folder]['output_summary']).name

    # Dictionary where all available cluster_summary dictionaries will be saved
    global_summary = {}

    # For each centroid
    for pose_name in poses_name_list:

        # Name of filtered pockets log file
        filtering_log_name =  f"{pose_name}_{filter_cavities_folder}_log.out"

        # Path to filtered pockets log file
        filtering_log_path = os.path.join(output_path, pose_name, filter_cavities_folder, filtering_log_name)

        # Find pockets in log file
        pocket_names = find_matching_lines(pattern=r'pocket\d+$', filepath=filtering_log_path)

        # If any pockets are found 
        if pocket_names is not None:

            # Dictionary with information for this cluster
            cluster_summary = {}

            # Path to this cluster's summary with information for all pocket found
            summary_path = os.path.join(output_path, pose_name, cavity_analysis_folder, pockets_summary_filename)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # list with pocket IDs for this cluster
            pocket_IDs = []

            # For each filtered pocket: pocket1, pocket4 ...
            for pocket_name in pocket_names:

                # Strip whitespaces
                pocket_name = pocket_name.strip()

                # Save entry with pocket information
                cluster_summary.update({pocket_name : all_pockets_summary[pocket_name]})

                # Search for the pocket ID in pocket name
                match = re.search(pocketID_pattern, pocket_name)

                # Append pocket ID to list
                pocket_IDs.append(int(match[1]))

            # Save pocket IDs
            cluster_summary.update({'pockets' : pocket_IDs})

            # Update global_summary
            global_summary.update({pose_name : cluster_summary})

    # Save global summary in YAML file
    if output_summary_path is None:
        output_summary_path = os.path.join(output_path, 'summary.yml')
    with open(output_summary_path, 'w') as outfile:
        yaml.dump(global_summary, outfile)
    outfile.close()

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

def main_wf(configuration_path, input_zip_path, output_path, output_summary_path):
    '''
    Main cavity analysis workflow. This workflow analyzes the cavities of the input structures structures, filters the cavities 
    according to pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------
    
        configuration_path  (str): path to input.yml 
        input_zip_path      (str): path to input zip file
        output_path         (str): path to output folder
        output_summary_path (str): path to output summary file


    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Enforce output_path if provided
    if output_path is not None:
        output_path = fu.get_working_dir_path(output_path, restart = conf.properties.get('restart', 'False'))
        conf.properties['working_dir_path'] = output_path
    else:
        output_path = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Create a folder for the extracted poses
    fu.create_dir(global_prop["step0_extract_poses"]["path"])

    # Enforce input_zip_path if provided
    if input_zip_path is not None:
        global_paths["step0_extract_poses"]["input_zip_path"] = input_zip_path
         
    # STEP 0: extract poses from zip file
    poses_path_list = fu.unzip_list(global_paths["step0_extract_poses"]["input_zip_path"], global_prop["step0_extract_poses"]["path"])
    poses_name_list = [Path(path).stem for path in poses_path_list]

    # Analyze cavities of each pdb 
    for path, name in zip(poses_path_list, poses_name_list):

        pose_prop = conf.get_prop_dic(prefix=name)
        pose_paths = conf.get_paths_dic(prefix=name)

        # Update input structure path
        pose_paths['step1_extractMolecule']['input_structure_path'] = path

        # STEP 1: Extract molecule
        global_log.info("step1_extract_molecule: Extract molecule from input structure")
        extract_molecule(**pose_paths['step1_extractMolecule'], properties=pose_prop["step1_extractMolecule"])

        # STEP 2: Cavity analysis
        global_log.info("step2_cavity_analysis: Compute cavities using fpocket")
        fpocket_run(**pose_paths['step2_cavity_analysis'], properties=pose_prop["step2_cavity_analysis"])

        # STEP 3: Filtering cavities
        global_log.info("step3_filter_cavities: Filter found cavities")
        fpocket_filter(**pose_paths['step3_filter_cavities'], properties=pose_prop["step3_filter_cavities"])

    # Create summary with available pockets per cluster 
    global_log.info("    Creating YAML summary file...")
    create_summary(poses_name_list, global_paths, output_path, output_summary_path)

    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % output_path)
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

    parser.add_argument('--input_zip', dest='input_zip_path',
                        help="Input zip file with all the poses (default: input_zip_path in step 0 of configuration file)",
                        required=False)
    
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)

    parser.add_argument('--output_summary', dest='output_summary_path',
                        help="Output summary path (default: /output_path/summary.yml)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path,
            input_zip_path = args.input_zip_path,
            output_path = args.output_path,
            output_summary_path = args.output_summary_path)

