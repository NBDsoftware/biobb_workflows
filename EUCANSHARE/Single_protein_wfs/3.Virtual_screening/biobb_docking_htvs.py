#!/usr/bin/env python3

# Importing all the needed libraries
import os
import re
import json
import time
import argparse
from glob import glob
from pathlib import Path
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_vs.utils.box import box
from biobb_chemistry.babelm.babel_convert import babel_convert
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens
from biobb_vs.vina.autodock_vina_run import autodock_vina_run

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

def find_matching_str(pattern, filepath):
    '''
    Finds str in file corresponding to a given pattern

    Inputs
    ------
        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------
        match_str (str): string matching the pattern or None if no piece matches the pattern
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

            return match.group(1)

    file.close()

    return None

def print_pockets(pockets_path, global_log):
    '''
    Print in log file all available pockets for each model found in the input folder for the pocket selection step
    
    Inputs
    ------
        pockets_path  (str): Path to pockets folder including file name
        global_log (Logger): Object of class Logger (dumps info to global log file of this workflow)

    Output
    ------
        Info dumped to log file
    '''

    global_log.error("    fpocket_select failed to select a pocket from input file.")
    global_log.error("    Check the selected pocket exists:")

    # Convert string to Path class
    pocketsLogPath = Path(pockets_path)

    # Path of pocket filtering step
    stepPath = pocketsLogPath.parent

    # Pattern for pocket filtering log files names
    logNamePattern = stepPath.name + "_log*.out"

    # Find all log files matching pattern
    logList = stepPath.rglob(logNamePattern)

    # Iterate through all log files -> print the centroid name + available pockets in the global log
    for log in logList:

        # Find pockets 
        pockets = find_matching_lines(pattern=r'pocket\d+$', filepath=str(log))

        # Find name of centroid
        centroidName = find_matching_str(pattern=r'/all_pockets_(\S+).zip', filepath=str(log))

        if None not in (pockets, centroidName):

            # Print to global log
            global_log.info("   Model: {}".format(centroidName))

            for pocket in pockets:
                global_log.info("        {}".format(pocket))

    return

def find_affinity(log_path):
    '''
    Find best binding affinity among different poses, returns None is no affinity is found
    
    Inputs
    ------
        log_path (str): paths of log file with results

    Outputs
    -------
        affinity (float): best affinity among poses  
    '''

    # Find affinity of best pose from log
    affinity = find_matching_str(pattern=r'\s+1\s+(\S+)\s+', filepath=log_path)

    if affinity is not None:
        affinity = float(affinity)
    
    return affinity

def create_ranking(ranking_path, num_ligands, affinities_list):
    '''
    Create file with ranking of ligands according to affinity

    Inputs
    ------
        affinities_list  (list): list with tuples -> (affinity, ligand identifier, ligand_ID)
        ranking_path      (str): path to output file
        num_ligands       (int): number of ligands to save in ranking file
    '''

    # Sort list according to affinity
    affinities_list = sorted(affinities_list)

    # Find number of ligands to save in global_summary
    ranking_length = min(num_ligands, len(affinities_list)) 

    # Exclude the rest
    affinities_list = affinities_list[:ranking_length]

    # Find parent folder
    step_folder = str(Path(ranking_path).parent)

    # Create folder if it does not exist
    if not os.path.exists(step_folder):
        os.makedirs(step_folder)
        
    # Create file
    with open(ranking_path, 'w') as file:

        # Write header
        file.write("Rank Affinity Ligand_ID \n")

        # Write ranking
        for i, affinity_tuple in enumerate(affinities_list):

            affinity, ligand_identifier, ligand_ID = affinity_tuple

            file.write("{}\t{}\t{}\n".format(i+1, affinity, ligand_ID))

    return 

def validate_step(*output_paths):
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
    
    # Erase files if they are empty
    if (validation_result == False):
        remove_files(*output_paths)

    return validation_result

def read_ligand_lib(ligand_lib_path):
    '''
    Read all ligand identifiers from ligand library file. 
    The expected format is one of the following:

    Format 1:

    ligand1_id
    ligand2_id
    .
    .
    .

    Format 2:

    ligand1_id  name_ligand1
    ligand2_id  name_ligand2
    .
    .
    .

    Where ligand_id is a SMILES and name_ligand is a string with the ligand name
    
    Inputs
    ------
        ligand_lib_path (str): path pointing to ligand library, including file name

    Output
    ------
        ligand_smiles   (list(str)): list of ligand SMILES
        ligand_names (list(str)): list with ligand names if any 
    '''

    ligand_smiles = []
    ligand_names = []

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Process every line
        for line in ligand_lib:

            line = line.split()

            # Append ligand ID to list
            ligand_smiles.append(line[0])

            # If there is a name, append it to ligand_names
            if len(line)>1:
                ligand_names.append(line[1])
            else:
                ligand_names.append(None)

    return ligand_smiles, ligand_names

def remove_files(*paths):
    """
    Remove temporal files if they exist

    Inputs
    ------

    paths   (str): variable number of paths to files including filename
    """

    tmp_files = []

    for path in paths:
        matching_files = glob(path)
        tmp_files.extend(matching_files)

    for file in tmp_files:
        os.remove(file)

    return

def add_suffix(original_path, suffix):
    '''
    Adds suffix to original_path before file extension. For example:

    Inputs
    ------

        suffix                (str):  suffix string 
        original_path  (Path class):  path to original file including filename

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

def add_ligand_info_to_paths(all_paths, ligand_name, ligand_index, *keywords):
    '''
    Modifies all_paths according to ligand_index and ligand_name if any. 
    
    Inputs
    ------

        all_paths      (dict): dictionary with all paths
        ligand_name     (str): name of ligand
        ligand_index    (int): ligand index
        keywords        (str): keywords of output/input paths that will be modified

    For example:

    path in all_paths = '/home/user/DockingWF/step4_drugbank/drug.sdf'
    ligand_name = None
    ligand_index = 4

    new_path = '/home/user/DockingWF/step4_drugbank/drug_4.sdf'
    '''

    # Prepare name of ligand
    if ligand_name is None:
        ligand_name = ""
    else:
        ligand_name = "_" + ligand_name
    
    # Create ligand identifier from ligand index and ligand name
    ligand_identifier = f"{ligand_index}{ligand_name}"

    # For all keys passed in keywords, modify path 
    add_suffix_to_paths(all_paths, ligand_identifier, *keywords)

def write_smiles(smiles, smiles_path):
    '''
    Writes a SMILES code into a file in step_path. 

    Inputs
    ------
        smiles              (str):  SMILES code
        smiles_path         (str):  smiles file path
    
    Output
    ------

        smiles_path (str): path to file containing SMILES code
    '''

    # Find step path 
    step_path = str(Path(smiles_path).parent)

    # If step directory has not been created yet -> create dir 
    if not os.path.exists(step_path):
        os.makedirs(step_path)

    # Save SMILES in tmp file inside step_path, overwrite if exists
    smiles_tmp_file = open(smiles_path, 'w')
    smiles_tmp_file.write(smiles)
    smiles_tmp_file.close()

def main_wf(configuration_path, ligand_lib_path):
    '''
    Main HTVS workflow. This workflow takes a ligand library, a pocket (defined by the output of a cavity analysis or some residues) 
    and a receptor to screen the pocket of the receptor using the ligand library (with Autodock).

    Inputs
    ------

        configuration_path   (str): path to input.yml 
        ligand_lib_path      (str): path to ligand library with SMILES

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)
    
    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    
    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()
        
    # STEP 1: Pocket selection from filtered list 
    global_log.info("step1_fpocket_select: Extract pocket cavity")
    fpocket_select(**global_paths["step1_fpocket_select"], properties=global_prop["step1_fpocket_select"])

    # Print available pockets if selection was not valid
    if not validate_step(global_paths["step1_fpocket_select"]["output_pocket_pdb"], global_paths["step1_fpocket_select"]["output_pocket_pqr"]):
        print_pockets(global_paths["step1_fpocket_select"]["input_pockets_zip"], global_log)
        return global_paths, global_prop, None

    # STEP 2: Generate box around selected cavity or residues
    global_log.info("step2_box: Generating cavity box")
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])

    # STEP 3: Prepare target protein for docking 
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 

    # STEP 4-5-6-7. For each ligand: prepare ligand, run docking, prepare poses

    # List where results will be saved
    affinities_list = []

    # Load drug list - not suitable for large libraries
    ligand_smiles, ligand_names = read_ligand_lib(ligand_lib_path)

    for ligand_index in range(len(ligand_smiles)):

        # ** Add ligand index to log file names NOTE: should we make a copy of global dicts before adapting to pycompps?
        global_prop["step4_babel_prepare_lig"].update({'prefix' : str(ligand_index)})
        global_prop["step5_autodock_vina_run"].update({'prefix' : str(ligand_index)})
        global_prop["step6_babel_prepare_pose"].update({'prefix' : str(ligand_index)})

        # ** Copy original paths to modify them (so that different processes do not use the same files)
        step4_paths = global_paths["step4_babel_prepare_lig"].copy()
        step5_paths = global_paths["step5_autodock_vina_run"].copy()
        step6_paths = global_paths["step6_babel_prepare_pose"].copy()

        # ** Add ligand name or index to paths
        add_ligand_info_to_paths(step4_paths, ligand_names[ligand_index], ligand_index, 'input_path', 'output_path')
        add_ligand_info_to_paths(step5_paths, ligand_names[ligand_index], ligand_index, 'input_ligand_pdbqt_path', 'output_pdbqt_path', 'output_log_path')
        add_ligand_info_to_paths(step6_paths, ligand_names[ligand_index], ligand_index, 'input_path', 'output_path')
        
        # STEP 4: Convert from smiles format to pdbqt format
        global_log.info("step4_babel_prepare_lig: Prepare ligand for docking")
        write_smiles(smiles = ligand_smiles[ligand_index], smiles_path = step4_paths["input_path"])
        try:
            babel_convert(**step4_paths, properties = global_prop["step4_babel_prepare_lig"])
            lastStep_successful = validate_step(step4_paths['output_path'])
        except:
            global_log.info(f"step4_babel_prepare_lig: Open Babel failed to convert ligand {ligand_smiles[ligand_index]} to pdbqt format")
            lastStep_successful = False

        # STEP 5: Autodock vina
        if lastStep_successful:
            try:
                global_log.info("step5_autodock_vina_run: Docking the ligand")            
                autodock_vina_run(**step5_paths, properties=global_prop["step5_autodock_vina_run"])
                lastStep_successful = validate_step(step5_paths['output_log_path'], 
                                                    step5_paths['output_pdbqt_path'])
                
                # Find best affinity among different poses
                affinity = find_affinity(log_path = step5_paths['output_log_path'])
                if affinity:
                    affinities_list.append((affinity, ligand_names[ligand_index], ligand_smiles[ligand_index]))
            except:
                global_log.info(f"step5_autodock_vina_run: Autodock Vina failed to dock ligand {ligand_smiles[ligand_index]}")
                lastStep_successful = False
                
        # STEP 6: Convert poses to PDB
        if lastStep_successful:
            try:
                global_log.info("step6_babel_prepare_pose: Converting ligand pose to PDB format")    
                babel_convert(**step6_paths, properties=global_prop["step6_babel_prepare_pose"])
            except:
                global_log.info(f"step6_babel_prepare_pose: Open Babel failed to convert pose for ligand {ligand_smiles[ligand_index]} to PDB format")
                lastStep_successful = False

        # Find paths to temporal log files from steps
        step4_stdout_path = os.path.join(global_prop["step4_babel_prepare_lig"]["path"], f'{ligand_index}_{global_prop["step4_babel_prepare_lig"]["step"]}_log.out')
        step4_stderr_path = os.path.join(global_prop["step4_babel_prepare_lig"]["path"], f'{ligand_index}_{global_prop["step4_babel_prepare_lig"]["step"]}_log.err')
        step5_stdout_path = os.path.join(global_prop["step5_autodock_vina_run"]["path"], f'{ligand_index}_{global_prop["step5_autodock_vina_run"]["step"]}_log.out')
        step5_stderr_path = os.path.join(global_prop["step5_autodock_vina_run"]["path"], f'{ligand_index}_{global_prop["step5_autodock_vina_run"]["step"]}_log.err')

        # Remove all temporal i/o files if they exist (all except the final pose)
        remove_files(step4_paths["input_path"],
                     step4_paths['output_path'],
                     step5_paths["output_pdbqt_path"],
                     step5_paths["output_log_path"],
                     step4_stdout_path, step4_stderr_path,
                     step5_stdout_path, step5_stderr_path)

    # STEP 8: Find top ligands (according to lowest affinity)
    global_log.info("step7_show_top_ligands: save identifiers of top ligands ranked by lowest affinity")  
    create_ranking(global_paths['step7_show_top_ligands']['output_csv_path'],
                   global_prop['step7_show_top_ligands']['number_top_ligands'],
                   affinities_list)
    
    # Timing information
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('  Ligand library: %s' % ligand_lib_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('--ligand_lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (SMILES format) per line.",
                        required=True)

    args = parser.parse_args()

    main_wf(args.config_path, args.ligand_lib)
