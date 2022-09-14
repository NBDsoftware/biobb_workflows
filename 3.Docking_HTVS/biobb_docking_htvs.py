#!/usr/bin/env python3

# Importing all the needed libraries
import os
import time
import argparse
from pathlib import Path, PurePath
from biobb_io.api.drugbank import drugbank
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_vs.utils.box import box
from biobb_chemistry.babelm.babel_convert import babel_convert
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens
from biobb_vs.vina.autodock_vina_run import autodock_vina_run

def main(config, system=None):

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(config, system)
    
    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    
    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 


    # STEP 1: Pocket selection from filtered list 

    # Write next action to global log
    global_log.info("step1_fpocket_select: Extract pocket cavity")

    # Action: pocket selection
    fpocket_select(**global_paths["step1_fpocket_select"], properties=global_prop["step1_fpocket_select"])


    # STEP 2: Generate box around selected cavity

    # Write next action to global log
    global_log.info("step2_box: Generating cavity box")
    
    # Action: box creation
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])


    # STEP 3: Prepare target protein for docking 
    
    # Write next action to global log     
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    
    # Action: Preparation of target protein
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 


    # STEP 4: Obtain small molecules from drugbank

    # Write next action to global log
    global_log.info("step4_drugbank: Extracting small molecules from library")
    
    # Properties and paths of step
    prop_drugbank  = global_prop["step4_drugbank"]
    paths_drugbank = global_paths["step4_drugbank"]

    drug_list = prop_drugbank['drugbank_list']

    num_drugs = len(drug_list)

    for i in range(num_drugs):

        # Add 'models' keyword and value to properties
        prop_drugbank.update({'drugbank_id':drug_list[i]})

        # Action: retrieve drug list
        drugbank(**paths_drugbank, properties=prop_drugbank)


        # STEP 5: Convert from sdf format to pdbqt format

        # Write next action to global log
        global_log.info("step5_babel_convert_prep_lig: Preparing small molecule (ligand) for docking")
    
        # Action: format conversion using Open Babel
        babel_convert(**global_paths["step5_babel_convert_prep_lig"], properties=global_prop["step5_babel_convert_prep_lig"])


        # STEP 6: Autodock vina

        # Write action to global log
        global_log.info("step6_autodock_vina_run: Running the docking")

        # Paths of step
        prop_autodock  = global_prop["step6_autodock_vina_run"]
        paths_autodock = global_paths["step6_autodock_vina_run"]

        # Get the parent paths
        pdbqt_parent_path = Path(paths_autodock['output_pdbqt_path']).parent
        log_parent_path = Path(paths_autodock['output_log_path']).parent

        # Define new names according to selected drug/ligand
        pdbqt_filename = str(pdbqt_parent_path) + '/output_vina_' + str(i) + '_' + drug_list[i] + '.pdbqt'
        log_filename   = str(log_parent_path) + '/output_vina_' + str(i) + '_' + drug_list[i] + '.log'

        # Update paths in dictionary
        paths_autodock.update({'output_pdbqt_path': pdbqt_filename})
        paths_autodock.update({'output_log_path': log_filename})

        # Action: Run of autodock vina
        autodock_vina_run(**paths_autodock, properties=prop_autodock)


        # STEP 7: Convert pose to PDB

        # Write action to global log
        global_log.info("step7_babel_convert_pose_pdb: Converting ligand pose to PDB format")

        # Paths of step
        prop_babel  = global_prop["step7_babel_convert_pose_pdb"]
        paths_babel = global_paths["step7_babel_convert_pose_pdb"]

        output_parent_path = Path(paths_babel['output_path']).parent   

        # Define new names according to selected drug/ligand
        input_filename  = pdbqt_filename
        output_filename = str(output_parent_path) + '/output_model_' + str(i) + '_' + drug_list[i] + '.pdb'

        # Update paths in dictionary
        paths_babel.update({'input_path': input_filename})
        paths_babel.update({'output_path': output_filename})

        # Action: Convert pose to PDB
        babel_convert(**paths_babel, properties=prop_babel)

    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    if system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Protein-ligand Docking tutorial using BioExcel Building Blocks")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
