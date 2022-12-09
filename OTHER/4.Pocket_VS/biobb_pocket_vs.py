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

from biobb_docking_htvs import main_wf as docking_htvs
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter
from biobb_structure_utils.utils.extract_residues import extract_residues
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_structure_utils.utils.extract_model import extract_model

def addPrefixToPath(original_path, prefix):
    '''
    Adds prefix to original_path file name. For example:

    original_path = Path('/home/user/ClusteringWF/output/step2_extract_models/output.pdb')
    prefix = model1

    new_path = '/home/user/ClusteringWF/output/step2_extract_models/model1_output.pdb'

    Inputs
    ------

        original_path  (Path class):  path to original file including filename
        prefix                (str):  prefix string
    
    Outputs
    ------- 

        new_path (str): path to new file including filename with prefix
    '''

    filename = original_path.name
    directory = original_path.parent

    new_name = prefix + "_" + filename

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

    global_log.info("    Available models and pockets after filtering:")

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

            global_log.info("    Model {}".format(modelName))

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

def extractPocketResidues(props, paths, pocket_residues, model):
    '''
    Extracts residues from PDB, returns path to the pdb file with the residues

    Inputs
    ------
        props               (dict): properties of step
        paths               (dict): paths of step
        pocket_residues      (str): string with residues ID separated by white spaces
        model                (str): pdb model name
    
    Output
    ------
        pocket_residues_path (str): path to file with residues

    '''

    pocket_residues = pocket_residues.split()

    pocket_residues = [int(res) for res in pocket_residues]

    props.update({'residues' : pocket_residues})


    original_output_path = Path(paths['output_residues_path'])
    original_input_path = Path(paths['input_structure_path'])

    pocket_residues_path = os.path.join(str(original_output_path.parent), model + "_" + original_output_path.name)
    input_structure_path = os.path.join(str(original_input_path.parent), model + "_" + original_input_path.name)

    paths.update({'output_residues_path' : pocket_residues_path})
    paths.update({'input_structure_path' : input_structure_path})
    
    extract_residues(**paths, properties=props)

    return pocket_residues_path

def renumberStructure(props, paths, input_structure_path, model):
    '''
    Renumber atoms and residues indexes

    Inputs
    ------
        props               (dict): properties of step
        paths               (dict): paths of step
        input_structure_path (str): input pdb path
        model                (str): pdb model name
    
    Output
    ------
        output_structure_path (str): path to renumbered pdb file 

    '''

    original_struct_path = Path(paths['output_structure_path'])
    original_map_path = Path(paths['output_mapping_json_path'])

    output_structure_path = os.path.join(str(original_struct_path.parent), model + "_" + original_struct_path.name)
    output_mapping_json_path = os.path.join(str(original_map_path.parent), model + "_" + original_map_path.name)

    paths.update({'input_structure_path' : input_structure_path})
    paths.update({'output_structure_path' : output_structure_path})
    paths.update({'output_mapping_json_path' : output_mapping_json_path})

    renumber_structure(**paths, properties=props)

    return output_structure_path

def extractModel(props, paths, model):
    '''
    Extract model

    Inputs
    ------
        props               (dict): properties of step
        paths               (dict): paths of step
        model                (str): pdb model name
    
    Output
    ------
        output_structure_path (str): path to pdb file with extracted model

    '''

    original_output_path = Path(paths['output_structure_path'])
    original_input_path = Path(paths['input_structure_path'])

    output_structure_path = os.path.join(str(original_output_path.parent), model + "_" + original_output_path.name)
    input_structure_path = os.path.join(str(original_input_path.parent), model + "_" + original_input_path.name)

    paths.update({'input_structure_path' : input_structure_path})
    paths.update({'output_structure_path' : output_structure_path})

    extract_model(**paths, properties=props)

    return output_structure_path


def main(args):
    
    start_time = time.time()

    # Set default value for 'last_step' arg
    if args.last_step is None:
        args.last_step = 'all'

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
    all_pdb_files = sorted(glob.glob(os.path.join(args.input_path, "*.pdb")))

    # Find number of pdbs in input folder
    for model in all_pdb_files:

        name = Path(model).stem

        zip_path = addPrefixToPath(default_zip, name)
        summary_path = addPrefixToPath(default_summary, name)

        # Update path dictionary
        paths_fpocket.update({'input_pdb_path':model})
        paths_fpocket.update({'output_pockets_zip':zip_path})
        paths_fpocket.update({'output_summary':summary_path})
        
        fpocket_run(**paths_fpocket, properties=prop_fpocket) 

        # Validate step
        if not validateStep(paths_fpocket["output_pockets_zip"], paths_fpocket["output_summary"]):
            global_log.info("    ERROR: no output from fpocket.")
            return 0

    if args.last_step == 'cavity':
        global_log.info("Cavity analysis completed.")
        return 

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
    for model in all_pdb_files:

        name = Path(model).stem

        # Modify the path names according to model index
        zip_input = addPrefixToPath(default_zip_input, name)
        summary_input = addPrefixToPath(default_summary_input, name)
        output_path = addPrefixToPath(default_zip_output, name)

        # Update path dictionary
        paths_filter.update({'input_pockets_zip':zip_input})
        paths_filter.update({'input_summary':summary_input})
        paths_filter.update({'output_filter_pockets_zip':output_path})

        fpocket_filter(**paths_filter, properties=prop_filter)

    if args.last_step == 'filtering':
        global_log.info("Cavity analysis completed.")
        return 

    # Find available models with pockets after filtering
    pockets_dict = printAvailablePockets(paths_filter["output_filter_pockets_zip"], global_log)

    # NOTE: cavity results should better inform autodock on each model 
    # right now if there is a good pocket somewhere (doesn't matter where) we proceed to next step

    # NOTE: parallelize over ligands

# STEP 3: Use active ligands to screen all the filtered cavities or pocket defined by residues

    global_log.info("step3_screen_models: screen pockets docking the probe ligands \n")

    if args.pocket_res:
        global_log.info("    Residues defining pocket are given: {} \n".format(args.pocket_res))
    else:
        global_log.info("    Residues defining pocket are not given, all pockets will be screened \n")

    model_history = []
    avg_affinity_history = []

    # For each model
    for model in pockets_dict.keys():

        # Print model name
        global_log.info(" ")
        global_log.info("  For model **** {} **** \n".format(model))

        # Find path to PDB model
        input_structure_path = os.path.join(args.input_path, model + ".pdb")
        
        # Find path to filtered pockets
        input_pockets_path = addPrefixToPath(default_zip_output, model)

        # Extract residues defining the pocket if list is provided -> dock ligands to that region
        if args.pocket_res:

            # Renumber structure
            global_log.info("   step3A_renumber_structure: Renumbering structure")

            prop_renumStruct = global_prop['step3A_renumber_structure'].copy()
            paths_renumStruct = global_paths['step3A_renumber_structure'].copy()

            output_structure_path = renumberStructure(props = prop_renumStruct,
                                                    paths = paths_renumStruct,
                                                    input_structure_path = input_structure_path,
                                                    model = model)

            # Extract model from structure
            global_log.info("   step3B_extract_model: Extracting model")

            prop_model = global_prop['step3B_extract_model'].copy()
            paths_model = global_paths['step3B_extract_model'].copy()

            output_structure_path = extractModel(props = prop_model,
                                                paths = paths_model,
                                                model = model)

            # Extract residues defining pocket
            global_log.info("   step3C_extract_residues: Extracting residues")

            prop_pocketRes = global_prop['step3C_extract_residues'].copy()
            paths_pocketRes = global_paths['step3C_extract_residues'].copy()

            pocket_residues_path = extractPocketResidues(props = prop_pocketRes,
                                                        paths = paths_pocketRes, 
                                                        pocket_residues = args.pocket_res, 
                                                        model = model)

            # Extract residues defining pocket
            global_log.info("   step4: Docking probe ligands")

            # Dock all ligands using residues to define box
            # NOTE: change this, now is a list of tuples (affinity, ligand identifier)
            output_dir, bestAffinities = docking_htvs(configuration_path = args.htvs_config_path, 
                                                                    ligand_lib_path = args.ligand_lib, 
                                                                    last_step = "all", 
                                                                    pocket_residues_path = pocket_residues_path,
                                                                    input_structure_path = output_structure_path)
            # New path to save docking_vs results
            new_output_dir = os.path.join(wdir, model + "_VS")

            # Move output_dir so it's not overwritten
            shutil.move(output_dir, new_output_dir)

            # Compute average affinity over docked ligands
            average_affinity = 0

            # Go over all tuples (affinity, ligand identifier)
            for ligand_info_tuple in bestAffinities:
                
                # Sum all the affinities 
                average_affinity += ligand_info_tuple[0]
            
            # Divide over total
            average_affinity = average_affinity/len(bestAffinities)

            # Print results for this pocket in global log
            global_log.info("     Yields average affinity of * {} * ".format(round(average_affinity, 3)))

            model_history.append(model)
            avg_affinity_history.append(average_affinity)

            # Print detailed affinities per compound
            for ligand_info_tuple in bestAffinities:

                global_log.info("        Ligand {} yields affinity of {}".format(ligand_info_tuple[1], round(ligand_info_tuple[0],3)))

        # If residues defining pocket are not provided -> dock ligands to all pockets found 
        else:

            # For each pocket
            for pocket_ID in pockets_dict[model]:

                # Dock all ligands and save results
                output_dir, bestAffinities = docking_htvs(configuration_path = args.htvs_config_path, 
                                                            ligand_lib_path = args.ligand_lib, 
                                                            last_step = "all", 
                                                            input_pockets_path = input_pockets_path, 
                                                            pocket_ID = pocket_ID,
                                                            input_structure_path = input_structure_path)

                # New path to save docking_vs results
                new_output_dir = os.path.join(wdir, model + "_" + str(pocket_ID) + "_VS")

                # Move output_dir so it's not overwritten
                shutil.move(output_dir, new_output_dir)

                # Compute average affinity over docked ligands
                average_affinity = 0

                # Go over all tuples (affinity, ligand identifier)
                for ligand_info_tuple in bestAffinities:
                    
                    # Sum all the affinities 
                    average_affinity += ligand_info_tuple[0]
                
                # Divide over total
                average_affinity = average_affinity/len(bestAffinities)

                # Print results for this pocket in global log
                global_log.info("        Pocket {} yields average affinity of {}".format(pocket_ID, round(average_affinity, 3)))

                # Print detailed affinities per compound
                for ligand_info_tuple in bestAffinities:

                    global_log.info("        Ligand {} yields affinity of {}".format(ligand_info_tuple[1], round(ligand_info_tuple[0],3)))
    
    if args.pocket_res:

        # Print ranking of models
        avg_affinity_history, model_history = zip(*sorted(zip(avg_affinity_history, model_history)))

        global_log.info(" ")
        global_log.info("*** FINAL RANKING ***")

        for model_index in range(len(model_history)):

            global_log.info("  Model {} gives average affinity of: {}".format(model_history[model_index], round(avg_affinity_history[model_index], 3)))


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

    return

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

    # Execute workflow until 'last_step' step -> all executes all steps (all is default)
    parser.add_argument('--until', dest='last_step', 
                        help="(Opt, default: all) Extent of the pipeline to execute (cavity, filtering, all)", 
                        required=False)

    parser.add_argument('--lig-lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (Ligand PDB code, SMILES or Drug Bank ID) per line.",
                        required=False)

    parser.add_argument('--pocket-res', dest='pocket_res',
                        help="List of residue IDs defining the pocket. If provided, the ligands will be docked in a box covering this region, ignoring pockets found.",
                        required=False)

    args = parser.parse_args()

    main(args)

