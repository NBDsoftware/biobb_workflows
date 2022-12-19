##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import re
import json
import glob
import time
import argparse
import shutil
from pathlib import Path

# Import pre-defined workflows
from biobb_docking_htvs import main_wf as docking_htvs    
#from biobb_docking_htvs_mp import main_wf as docking_htvs
from biobb_clustering_cavity_analysis import main_wf as cavity_analysis

# Import biobb
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_structure_utils.utils.extract_residues import extract_residues
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_structure_utils.utils.extract_model import extract_model

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

def findPocketResidues(pocket_residues_list, input_structure_path, modelName, global_log, global_prop, global_paths):
    '''
    Fix the input structure (renumbering and extracting model if there is more than one) and extract relevant residues
    that define the pocket to be screened

    Inputs
    ------

        pocket_residues_list (list): list of residue ID numbers that define the pocket
        input_structure_path  (str): path to input structure
        modelName             (str): model name 
        global_log                 : global logger of pocket VS wf
        global_prop                : global properties of pocket VS wf
        global_paths               : global paths of pocket VS wf

    Output
    ------

        pocket_residues_path   (str): path to pdb with relevant residues
        output_structure_path  (str): path to pdb from which the residues are extracted (input pdb after renumbering and model extraction)

    '''

    # We have found a list of residues
    global_log.info("    Residues defining pocket are given: {} \n".format(pocket_residues_list))

    # Renumber structure
    global_log.info("   step3A_renumber_structure: Renumbering structure")

    prop_renumStruct = global_prop['step3A_renumber_structure'].copy()
    paths_renumStruct = global_paths['step3A_renumber_structure'].copy()

    output_structure_path = renumberStructure(props = prop_renumStruct,
                                            paths = paths_renumStruct,
                                            input_structure_path = input_structure_path,
                                            model = modelName)

    # Extract model from structure
    global_log.info("   step3B_extract_model: Extracting model")

    prop_model = global_prop['step3B_extract_model'].copy()
    paths_model = global_paths['step3B_extract_model'].copy()

    output_structure_path = extractModel(props = prop_model,
                                        paths = paths_model,
                                        model = modelName)

    # Extract residues defining pocket
    global_log.info("   step3C_extract_residues: Extracting residues")

    prop_pocketRes = global_prop['step3C_extract_residues'].copy()
    paths_pocketRes = global_paths['step3C_extract_residues'].copy()

    pocket_residues_path = extractPocketResidues(props = prop_pocketRes,
                                                paths = paths_pocketRes, 
                                                pocket_residues = pocket_residues_list, 
                                                model = modelName)

    return pocket_residues_path, output_structure_path

def analyzeDockingResults(wdir, modelName, docking_output_dir, docking_summary, global_log, avg_affinity_history, pocket_ID = None):
    '''
    Analyze the results from a docking screening: move output folder so they are not overwritten, compute average affinity, 
    print detailed results for each ligand

    Inputs
    ------
    
        wdir                  (str): current working dir path
        modelName             (str): model name
        docking_output_dir    (str): docking workflow output path
        docking_summary      (dict): summary with docking results 
        global_log         (Logger): global logger of pocket VS wf
        avg_affinity_history (list): 
        pocket_ID             (int):

    Output
    ------

    '''

    # New path to save docking_vs results
    new_output_dir = os.path.join(wdir, modelName + "_VS")

    # Move output_dir so it's not overwritten
    shutil.move(docking_output_dir, new_output_dir)

    # Compute average affinity over docked ligands
    average_affinity = 0

    # Ligand identifiers (index + ligand ID)
    ligand_identifiers = docking_summary.keys()

    # Go over all tuples (affinity, ligand identifier)
    for ligand in ligand_identifiers:
        
        # Sum all the affinities 
        average_affinity += docking_summary[ligand]['affinity']
    
    # Divide over total to find average
    average_affinity = average_affinity/len(ligand_identifiers)

    # We add the average affinity to the history
    avg_affinity_history.append(average_affinity)

    if pocket_ID is None:
        # Print results for this model in global log
        global_log.info("     Yields average affinity of * {} * ".format(round(average_affinity, 3)))
    else:
        # Print results for this pocket in global log
        global_log.info("       Pocket {} yields average affinity of {}".format(pocket_ID, round(average_affinity, 3)))

    # Print detailed affinities per compound
    for ligand in ligand_identifiers:

        global_log.info("        Ligand {} yields affinity of {}".format(ligand, round(docking_summary[ligand]['affinity'],3)))

    return

def pocketsVirtualScreening(configuration_path, input_structures_path, ligand_lib_path, pocket_residues_list, 
                            pockets_path, cavities_summary, global_conf, global_log):

    '''
    Virtual Screening of pockets using ligand library and (optionally) information regarding a specific region of interest

    Inputs
    ------

        configuration_path      (str): htvs configuration path
        input_structures_path   (str): path to input structures
        ligand_lib_path         (str): path to ligand library
        pocket_residues_list   (list): list of pocket residues
        pockets_path            (str): path to filtered pockets
        cavities_summary       (dict): dictionary with models and pockets
        global_conf                  : global configuration loaded as a configuration object
        global_log                   : global log object
    
    Output
    ------

        avg_affinity_history (list): average affinity among all ligands for each pocket
        pocket_list          (list): corresponding pocket/model

    '''

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = global_conf.get_prop_dic(global_log=global_log)
    global_paths = global_conf.get_paths_dic()

    # Find global log from properties
    wdir = global_conf.get_working_dir_path()

    # Average affinity results
    avg_affinity_history = []

    # List of paths to filtered pockets
    all_filtered_pockets = sorted(glob.glob(os.path.join(pockets_path, "*.zip")))

    # List of all model names 
    all_models = sorted(cavities_summary.keys())

    # For each model
    for model_index in range(len(all_models)):

        # Model Name
        modelName = all_models[model_index]

        # Path to PDB
        input_structure_path = os.path.join(input_structures_path, modelName + ".pdb")
        
        # Path to filtered pockets 
        pockets_path = all_filtered_pockets[model_index]

        # Print model name
        global_log.info(" ")
        global_log.info("  For model **** {} **** \n".format(modelName))

        # If list is provided 
        # -> dock ligands to that region
        if pocket_residues_list:
            
            # Extract residues defining pocket
            pocket_residues_path, output_structure_path = findPocketResidues(pocket_residues_list, 
                                                                            input_structure_path,  
                                                                            modelName, 
                                                                            global_log, 
                                                                            global_prop, 
                                                                            global_paths)

            # Dock probe ligands
            global_log.info("   step4: Docking probe ligands")

            # Dock all ligands using residues to define box
            dvs_paths, dvs_prop, dvs_summary = docking_htvs(configuration_path = configuration_path, 
                                                      ligand_lib_path = ligand_lib_path, 
                                                      last_step = "all", 
                                                      pocket_residues_path = pocket_residues_path,
                                                      input_structure_path = output_structure_path)

            # Move output folder, compute average affinity, print detailed results for each ligand
            analyzeDockingResults(wdir, modelName, dvs_prop['step1_fpocket_select']['working_dir_path'], dvs_summary, global_log, avg_affinity_history)

        # If residues defining pocket are not provided 
        # -> dock ligands to all pockets found 
        # NOTE: this should be improved, maybe only best, or those close to a certain residue 
        else:

            # For each pocket
            for pocket_ID in cavities_summary[modelName]['pockets']:

                # Dock all ligands and save results
                dvs_paths, dvs_prop, dvs_summary = docking_htvs(configuration_path = configuration_path, 
                                                            ligand_lib_path = ligand_lib_path, 
                                                            last_step = "all", 
                                                            input_pockets_path = pockets_path, 
                                                            pocket_ID = pocket_ID,
                                                            input_structure_path = input_structure_path)

                # Move output folder, compute average affinity, print detailed results for each ligand
                analyzeDockingResults(wdir, modelName, dvs_prop['step1_fpocket_select']['working_dir_path'], dvs_summary, global_log, avg_affinity_history)
    
    if pocket_residues_list:

        # Print ranking of models
        avg_affinity_history, all_models = zip(*sorted(zip(avg_affinity_history, all_models)))

        global_log.info(" ")
        global_log.info("*** FINAL RANKING ***")

        for model_index in range(len(all_models)):

            global_log.info("  Model {} gives average affinity of: {}".format(all_models[model_index], round(avg_affinity_history[model_index], 3)))

    return avg_affinity_history, all_models

def main_wf(pocket_vs_config, cavity_analysis_config, htvs_config, 
         input_structures_path, ligand_lib_path, pocket_residues_list = None, last_step = None):
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and analyzes the cavities of the most representative
    structures. Then filters the cavities according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        pocket_vs_config       (str): path to main input.yml with configuration for pocket VS wf
        cavity_analysis_config (str): path to cavity_input.yml with configuration for cavity analysis wf
        htvs_config            (str): path to htvs_input.yml with configuration for htvs wf
        input_structures_path  (str): path to the folder with the input structures in pdb format (*.pdb)
        ligand_lib_path        (str): path to ligand library
        pocket_residues_list  (list): list with residues IDs defining pocket
        last_step              (str): last step of the workflow to execute (cavity, filtering, all)

    Outputs
    -------

        /output folder
        working_dir_path (str): path to working directory 
    '''

    # NOTE: calculate COM of pocket, calculate COM of residues 
    # -> only screen the pocket defined by the residues if we have found a pocket there that passed the filtering

    # NOTE: if one ligand gives a very bad result maybe we can skip the rest? How could we implement this

    start_time = time.time()

    # Set default value for 'last_step' arg
    if last_step is None:
        last_step = 'all'

    # Receiving the main input configuration file (YAML)
    conf = settings.ConfReader(pocket_vs_config)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

# STEP 1-2: Cavity analysis and filtering 

    # Write next action to global log
    global_log.info("Cavity analysis and filtering: Compute protein cavities for each structure using fpocket and filter them")

    # Use cavity analysis workflow
    cavities_paths, cavities_prop, cavities_summary = cavity_analysis(configuration_path = cavity_analysis_config, 
                                                                        clustering_path = input_structures_path, 
                                                                        last_step = last_step)

    # If this is last step, exit
    if last_step == 'cavity' or last_step == 'filtering':
        global_log.info("Cavity analysis completed.")
        return 

# STEP 3: Use active ligands to screen the filtered cavities or pocket defined by residues
    
    global_log.info("Pocket screening: screen pockets docking the probe ligands \n")

    # Find path with pockets from cavity analysis
    pockets_path = cavities_prop['step4_filter_cavities']['path']

    avg_affinity_history, all_models = pocketsVirtualScreening(configuration_path = htvs_config, 
                                                    input_structures_path = input_structures_path, 
                                                    ligand_lib_path = ligand_lib_path, 
                                                    pocket_residues_list = pocket_residues_list, 
                                                    pockets_path = pockets_path, 
                                                    cavities_summary = cavities_summary, 
                                                    global_conf = conf,
                                                    global_log = global_log)

    # New path to save cavity analysis results
    new_output_dir = os.path.join(conf.get_working_dir_path(), "cavity_output")

    # Move cavity output to output/ folder
    shutil.move(cavities_prop['step4_filter_cavities']['working_dir_path'], new_output_dir)

    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % pocket_vs_config)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('--main-config', dest='main_config_path',
                        help="Configuration file for pocket VS (YAML)", 
                        required=True)

    parser.add_argument('--cavity-config', dest='cavity_config_path',
                        help="Configuration file for cavity analysis (YAML)", 
                        required=True)

    parser.add_argument('--htvs-config', dest='htvs_config_path',
                        help="Configuration file for HTVS (YAML)", 
                        required=True)
    
    parser.add_argument('--input', dest='input_path',
                        help="Path to folder with pdb models (For example PDB centroids from clustering an MD trajectory)", 
                        required=True)

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

    main_wf(pocket_vs_config = args.main_config_path, 
            cavity_analysis_config = args.cavity_config_path, 
            htvs_config = args.htvs_config_path, 
            input_structures_path = args.input_path, 
            last_step = args.last_step, 
            ligand_lib_path = args.ligand_lib,
            pocket_residues_list = args.pocket_res)

