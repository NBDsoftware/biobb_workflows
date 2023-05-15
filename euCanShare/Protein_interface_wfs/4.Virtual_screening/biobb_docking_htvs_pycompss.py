#!/usr/bin/env python3

# Importing all the needed libraries
import os
import re
import time
import shutil
import argparse
from pathlib import Path

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_adapters.pycompss.biobb_vs.utils.box import box
from biobb_adapters.pycompss.biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_adapters.pycompss.biobb_vs.vina.autodock_vina_run import autodock_vina_run
from biobb_adapters.pycompss.biobb_chemistry.babelm.babel_convert import babel_convert
from biobb_adapters.pycompss.biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens

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

def get_affinity(log_path):
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

def create_ranking(affinities, global_paths, global_prop, output_path):
    '''
    Create file with ranking of ligands according to affinity

    Inputs
    ------
        affinities      (list): list with tuples -> (affinity, ligand identifier, ligand_ID)
        paths           (dict): dictionary with paths of the step
        properties      (dict): dictionary with properties of the step
    '''
    # Step names    
    pose_step_name = 'step6_babel_prepare_pose'
    ranking_step_name = 'step7_show_top_ligands'

    # Read paths
    pose_file_name = Path(global_paths[ranking_step_name]['input_pose_path']).name
    ranking_path = global_paths[ranking_step_name]['output_csv_path']

    # Read properties
    num_top_ligands = global_prop[ranking_step_name]['num_top_ligands']
    keep_poses = global_prop[ranking_step_name]['keep_poses']

    # Sort list according to affinity
    affinities = sorted(affinities)

    # Find number of ligands to save in ranking
    ranking_length = min(num_top_ligands, len(affinities)) 

    # Exclude the rest
    affinities = affinities[:ranking_length]

    # Find step path
    step_path = global_prop[ranking_step_name]['path']

    # Create folder if it does not exist
    if not os.path.exists(step_path):
        os.makedirs(step_path)
        
    # Create ranking 
    with open(ranking_path, 'w') as file:

        # Write header
        file.write("Rank Affinity Ligand_name Ligand_ID \n")

        # For each ligand
        for rank, affinity_tuple in enumerate(affinities):

            affinity, ligand_name, ligand_ID = affinity_tuple

            # Write ranking line
            file.write(f"{rank+1}\t{affinity}\t{ligand_name}\t{ligand_ID}\n")

            # If keep_poses is True, copy pose to step folder
            if keep_poses:

                # Find pose path
                pose_path = os.path.join(output_path, ligand_name, pose_step_name, pose_file_name)

                # New pose path
                new_pose_path = os.path.join(step_path, f"{ligand_name}_poses.pdb")

                # Copy pose 
                shutil.copyfile(pose_path, new_pose_path)

    return 

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
        ligand_lib_path (str): path to ligand library file

    Output
    ------
        ligand_smiles   (list(str)): list of ligand SMILES
        ligand_names    (list(str)): list with ligand names
    '''

    ligand_smiles = []
    ligand_names = []

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Process every line
        for index, line in enumerate(ligand_lib):

            line = line.split()

            # Append ligand ID to list
            ligand_smiles.append(line[0])

            # If there is a name, append it to ligand_names
            if len(line)>1:
                ligand_names.append(line[1])
            else:
                ligand_names.append(str(index))

    return ligand_smiles, ligand_names

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

def get_affinities(ligand_smiles, ligand_names, global_paths, output_path):

    """
    Reads autodock log files to find best affinity for each ligand

    Inputs
    ------

        ligand_smiles   (list): list of SMILES codes
        ligand_names    (list): list of ligand names
        global_prop     (dict): global properties dictionary
        global_paths    (dict): global paths dictionary
    
    Output
    ------

        affinities      (list): list of tuples (affinity, ligand_name, ligand_smiles)
    """

    # AutoDock step name
    autodock_step_name = 'step5_autodock_vina_run'

    # Find AutoDock log name
    log_name = Path(global_paths[autodock_step_name]['output_log_path']).name

    # List where best affinity for each ligand will be stored
    affinities = []
    
    for smiles, name in zip(ligand_smiles, ligand_names):

        # Find AutoDock log path
        log_path = os.path.join(output_path, name, autodock_step_name, log_name)

        # Find best affinity among different poses
        affinity = get_affinity(log_path = log_path)

        if affinity:
            affinities.append((affinity, name, smiles))
    
    return affinities

def clean_output(ligand_names, output_path):
    """
    Removes all ligand subdirectories in the output directory
    """

    # Remove all ligand subdirectories - NOTE: this could be potentially dangerous
    for name in ligand_names:
        ligand_path = os.path.join(output_path, name)
        shutil.rmtree(ligand_path)
    
def main_wf(configuration_path, ligand_lib_path, structure_path, input_pockets_zip, pocket, output_path, num_top_ligands):
    '''
    Main HTVS workflow. This workflow takes a ligand library, a pocket (defined by the output of a cavity analysis or some residues) 
    and a receptor to screen the pocket of the receptor using the ligand library (with AutoDock).

    Inputs
    ------

        configuration_path   (str): path to input.yml 
        ligand_lib_path      (str): path to ligand library with SMILES
        structure_path       (str): path to receptor structure
        input_pockets_zip    (str): path to input pockets zip file
        pocket               (str): pocket name
        output_path          (str): path to output directory
        num_top_ligands      (int): number of top ligands to be saved

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties
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
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Enforce input_pockets_zip if provided
    if input_pockets_zip is not None:
        global_paths['step1_fpocket_select']['input_pockets_zip'] = input_pockets_zip
    
    # Enforce pocket if provided
    if pocket is not None:
        global_prop['step1_fpocket_select']['pocket'] = pocket

    # STEP 1: Pocket selection from filtered list 
    global_log.info("step1_fpocket_select: Extract pocket cavity")
    fpocket_select(**global_paths["step1_fpocket_select"], properties=global_prop["step1_fpocket_select"])

    # STEP 2: Generate box around selected cavity or residues
    global_log.info("step2_box: Generating cavity box")
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])

    # Enforce structure_path if provided
    if structure_path is not None:
        global_paths['step3_str_check_add_hydrogens']['input_structure_path'] = structure_path
        
    # STEP 3: Prepare target protein for docking 
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 

    # STEP 4-5-6-7. For each ligand: prepare ligand, run docking, prepare poses

    # Load drug list - NOTE: the whole library is loaded into memory
    ligand_smiles, ligand_names = read_ligand_lib(ligand_lib_path)

    for smiles, name in zip(ligand_smiles, ligand_names):

        # Add ligand name to properties and paths
        ligand_prop = conf.get_prop_dic(prefix=name)
        ligand_paths = conf.get_paths_dic(prefix=name)

        # Write smiles to file
        write_smiles(smiles = smiles, smiles_path = ligand_paths['step4_babel_prepare_lig']['input_path'])

        # STEP 4: Convert from smiles to pdbqt
        global_log.info("step4_babel_prepare_lig: Prepare ligand for docking")
        babel_convert(**ligand_paths['step4_babel_prepare_lig'], properties = ligand_prop["step4_babel_prepare_lig"])

        # Update common paths - NOTE: different processes will try to access the same files here ...
        ligand_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path'] = global_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path']
        ligand_paths['step5_autodock_vina_run']['input_box_path'] = global_paths['step5_autodock_vina_run']['input_box_path']

        # STEP 5: Docking with Autodock Vina
        global_log.info("step5_autodock_vina_run: Docking the ligand")            
        autodock_vina_run(**ligand_paths['step5_autodock_vina_run'], properties=ligand_prop["step5_autodock_vina_run"])
                
        # STEP 6: Convert poses to PDB
        global_log.info("step6_babel_prepare_pose: Converting ligand pose to PDB format")    
        babel_convert(**ligand_paths['step6_babel_prepare_pose'], properties=ligand_prop["step6_babel_prepare_pose"])

    # Find the best affinity for each ligand
    affinities = get_affinities(ligand_smiles, ligand_names, global_paths, output_path)

    # Enforce num_top_ligands if specified
    if num_top_ligands is not None:
        global_prop['step7_show_top_ligands']['num_top_ligands'] = int(num_top_ligands)

    # STEP 7: Find top ligands 
    global_log.info("step7_show_top_ligands: create ranking and save poses for top ligands")  
    create_ranking(affinities, global_paths, global_prop, output_path)
    
    # Clean up the output folder
    clean_output(ligand_names, output_path)

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
    
    parser.add_argument('--structure_path', dest='structure_path',
                        help="Path to file with target structure (PDB format)",
                        required=False)
    
    parser.add_argument('--input_pockets_zip', dest='input_pockets_zip',
                        help="Path to file with pockets in a zip file",
                        required=False)

    parser.add_argument('--pocket', dest='pocket',
                        help="Pocket number to be used (default: 1)",
                        required=False)
    
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    parser.add_argument('--num_top_ligands', dest='num_top_ligands',
                        help="Number of top ligands to be saved (default: corresponding value in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            ligand_lib_path = args.ligand_lib,
            structure_path = args.structure_path,
            input_pockets_zip = args.input_pockets_zip,
            pocket = args.pocket,
            output_path = args.output_path,
            num_top_ligands = args.num_top_ligands)
