##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import re
import glob
import time
import argparse
import shutil
from pathlib import Path

from biobb_docking_htvs import docking_htvs
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def addModelPrefix(original_path, name):
    '''
    Adds name to original_path. For example:

    original_path = Path('/home/user/ClusteringWF/output/step2_extract_models/output.pdb')
    name = model1

    new_path = '/home/user/ClusteringWF/output/step2_extract_models/model1_output.pdb'

    Inputs
    ------

        original_path  (Path class from pathlib):  path to original file including filename
        i  (int):  integer 
    '''

    filename = original_path.name
    directory = original_path.parent

    new_name = name + "_" + filename

    new_path = os.path.join(directory, new_name)

    return new_path

def moveFileIfItExists(origin_src, dest_src):
    '''
    Checks the existence of 'origin_src' file and moves it to 'dest_src' if it exists

    Inputs
    ------

        origin_src  (str):  path to original file including filename
        dest_src    (str):  path to destination file including filename
    '''

    if (os.path.exists(origin_src)):
        shutil.move(origin_src, dest_src)
    
    return

def validateStep(*output_paths):
    '''
    Check all output files exist and are not empty
    
    Inputs
    ------
        *output_paths (str): variable number of paths to output file/s

    Output
    ------
        validation_result (bool): result of validation
    '''

    # Initialize value 
    validation_result = True

    # Check existence of files
    for path in output_paths:
        validation_result = validation_result and os.path.exists(path)

    # Check files are not empty if they exist
    if (validation_result):

        for path in output_paths:
            file_not_empty = os.stat(path).st_size > 0
            validation_result = validation_result and file_not_empty

    return validation_result

def printAvailablePockets(pockets_path, global_log):
    '''
    Print in log file all available pockets for each model found in the input folder for the pocket selection step
    
    Inputs
    ------
        pockets_path  (str): Path to pockets folder including file name
        global_log (Logger): Object of class Logger (dumps info to global log file of this workflow)
    Output
    ------
        Info dumped to log file
        Dictionary with models and pockets  {modelname1 : [1,2,3], modelname2 : [2,4]}
    '''
    # Dictionary where pockets will be saved
    pockets_dict = {}

    # Convert string to Path class
    pocketsLogPath = Path(pockets_path)

    # Path of pocket filtering step
    stepPath = pocketsLogPath.parent

    # Pattern for pocket filtering log files names
    logNamePattern = stepPath.name + "_log*.out"

    # Find all log files matching pattern
    logList = stepPath.rglob(logNamePattern)

    # Iterate through all log files -> print the model name + available pockets in the global log
    for log in logList:

        # Find model name NOTE: hardcoded name of step
        modelName = findMatchingLine(pattern=r'step1_cavity_analysis/(\S+)_all_pockets.zip', filepath=str(log))

        # Find pockets 
        pockets = findMatchingLines(pattern=r'pocket\d+$', filepath=str(log))

        # Print to global log
        if modelName is not None:

            global_log.info("   Model {}".format(modelName))

            if pockets is not None:

                pocket_indices = []

                for pocket in pockets:

                    pattern = r'pocket(\d+)'

                    match = re.search(pattern, pocket)
                    
                    global_log.info("        {}".format(pocket))

                    pocket_indices.append(int(match[1]))

                pockets_dict.update({modelName : pocket_indices})

    return pockets_dict

def findMatchingLine(pattern, filepath):
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

def findMatchingLines(pattern, filepath):
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

def findNumberInString(string):
    '''
    Finds and returns first integer number in string. If no integer is found a 0 is returned.

    Inputs
    ------
        string (str): string to search int in
    
    Output
    ------
        number (int): integer found in string or 0 if no integer is found.
    '''

    # Pattern corresponding to a digit of one or more characters
    digitPattern = r'\d+'

    # Search for the first match in string
    number = re.search(digitPattern, string)
    
    if (number == None):
        # If there is no match, return a 0
        number = 0
    else:
        # If there is a match, convert match object to string
        number = number[0]

    return int(number)


def main(args):
    
    start_time = time.time()

    # Set default value for 'to_do' arg
    if args.to_do is None:
        args.to_do = 'all'

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.main_config_path )
    wdir = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=wdir, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

# STEP 1: Cavity analysis 

    # Write next action to global log
    global_log.info("step1_cavity_analysis: Compute protein cavities for each structure using fpocket")

    # Action: Search for cavities in each pdb using fpocket
    # Properties and paths of step
    prop_fpocket = global_prop["step1_cavity_analysis"]
    paths_fpocket = global_paths["step1_cavity_analysis"]

    # Default paths from input.yml
    default_zip = Path(paths_fpocket['output_pockets_zip'])
    default_summary = Path(paths_fpocket['output_summary'])

    # Generic name for pdb file in centroids folder
    all_pdb_files = os.path.join(args.input_path, "*.pdb")

    # Find number of pdbs in input folder
    for model in glob.glob(all_pdb_files):

        name = Path(model).stem

        zip_path = addModelPrefix(default_zip, name)
        summary_path = addModelPrefix(default_summary, name)

        # Update path dictionary
        paths_fpocket.update({'input_pdb_path':model})
        paths_fpocket.update({'output_pockets_zip':zip_path})
        paths_fpocket.update({'output_summary':summary_path})
        
        fpocket_run(**paths_fpocket, properties=prop_fpocket) 

        # Validate step
        if not validateStep(paths_fpocket["output_pockets_zip"], paths_fpocket["output_summary"]):
            global_log.info("    ERROR: no output from fpocket.")
            return 0
    
# STEP 2: Filtering cavities

    # Write next action to global log
    global_log.info("step2_filter_cavities: Filter found cavities")

    # Properties and paths of step
    prop_filter = global_prop["step2_filter_cavities"]
    paths_filter = global_paths["step2_filter_cavities"]

    # Default paths from input.yml
    default_zip_input = Path(paths_filter['input_pockets_zip'])
    default_summary_input = Path(paths_filter['input_summary'])
    default_zip_output = Path(paths_filter['output_filter_pockets_zip'])

    # For each pdb model
    for model in glob.glob(all_pdb_files):

        name = Path(model).stem

        # Modify the path names according to model index
        zip_input = addModelPrefix(default_zip_input, name)
        summary_input = addModelPrefix(default_summary_input, name)
        output_path = addModelPrefix(default_zip_output, name)

        # Update path dictionary
        paths_filter.update({'input_pockets_zip':zip_input})
        paths_filter.update({'input_summary':summary_input})
        paths_filter.update({'output_filter_pockets_zip':output_path})

        fpocket_filter(**paths_filter, properties=prop_filter)

# STEP 3: Use active ligands to screen the cavities

    # NOTE it would be nice to add a region filter...
    # Find available models with pockets after filtering
    pockets_dict = printAvailablePockets(paths_filter["output_filter_pockets_zip"], global_log)

    global_log.info("step3_screen_models:")

    # For each model
    for model in pockets_dict.keys():

        # Print model name
        global_log.info("    For model {}".format(model))

        # For each pocket
        for pocket_ID in pockets_dict[model]:

            # Find path to filtered pockets
            input_pockets_path = addModelPrefix(default_zip_output, model)

            # Find path to PDB model
            input_structure_path = os.path.join(args.input_path, model + ".pdb")

            # Dock all ligands and save results
            output_dir, bestAffinities, _ = docking_htvs(configuration_path = args.htvs_config_path, 
                                                                ligand_lib_path = args.ligand_lib, 
                                                                last_step = "all", 
                                                                input_pockets_path = input_pockets_path, 
                                                                pocket_ID = pocket_ID,
                                                                input_structure_path = input_structure_path)

            # New path to save docking_vs results
            new_output_dir = os.path.join(wdir, model + "_" + str(pocket_ID) + "_VS")

            # Move output_dir so it's not overwritten
            shutil.move(output_dir, new_output_dir)

            # Average affinity
            average_affinity = sum(bestAffinities)/len(bestAffinities)

            # Print results for this pocket in global log
            global_log.info("        Pocket {} yields average affinity of {}".format(pocket_ID, average_affinity))

            # Print detailed affinities per compound
            
    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % args.main_config_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('--main-config', dest='main_config_path',
                        help="Configuration file (YAML)", 
                        required=True)

    parser.add_argument('--htvs-config', dest='htvs_config_path',
                        help="Configuration file (YAML)", 
                        required=True)
    
    parser.add_argument('--input', dest='input_path',
                        help="Path to folder with pdb models (For example PDB centroids from clustering an MD trajectory)", 
                        required=True)

    # Execute workflow until 'to_do' step -> all executes all steps (all is default)
    parser.add_argument('--until', dest='to_do', 
                        help="(Opt, default: all) Extent of the pipeline to execute (cavity, all)", 
                        required=False)

    parser.add_argument('--lig-lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (Ligand PDB code, SMILES or Drug Bank ID) per line.",
                        required=False)

    args = parser.parse_args()

    main(args)

