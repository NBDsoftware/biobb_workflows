#!/usr/bin/env python3

# Importing all the needed libraries
import os
import glob
import time
import shutil
import argparse
from pathlib import Path 
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.fix_chirality import fix_chirality
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.genion import genion
from biobb_gromacs.gromacs.mdrun import mdrun
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str

# Define a function to copy files from several paths to a single destination
# Adding a prefix to the file name
def copy_files_with_prefix(source_paths, destination_path, prefix):
    """
    Copy files from several source paths to a single destination. Add prefix to the original file names.

    Inputs
    ------

        source_paths (dict): dictionary with all source paths
        destination_path (str): destination path
        prefix (str): prefix to add to the original file names
    
    Outputs
    -------

        /output folder
        None
    """

    # Iterate over the source paths dictionary
    for file, path in source_paths.items():

        # Get the file name
        file_name = Path(path).name

        # Add the prefix to the file name
        new_file_name = prefix + "_" + file_name

        # Create new path
        new_path = str(Path(destination_path).joinpath(new_file_name))

        # Copy the file
        shutil.copy(path, new_path)
    
def main_wf(configuration_path):
    '''
    Main protein-protein docking pose refinement workflow. It takes as input a zip file containing PDB files with the protein-protein docking poses and a 
    YAML configuration file with all the parameters needed to run the workflow.

    Inputs
    ------

        configuration_path (str): path to input.yml

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
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Create a folder for the extracted poses
    fu.create_dir(global_prop["step0_extract_poses"]["path"])

    # Create a folder for the merged results
    fu.create_dir(global_prop["step18_merge_results"]["path"])

    # STEP 0: extract poses from zip file

    # Extract the poses from the zip file
    poses_paths = fu.unzip_list(global_paths["step0_extract_poses"]["input_zip_path"], global_prop["step0_extract_poses"]["path"])

    # Iterate over the poses
    for pose_path in poses_paths:

        # Copy the pose to new path
        shutil.copy(pose_path, global_paths["step0_extract_poses"]["output_pdb_path"])

        # STEP 1: model missing heavy atoms of side chains
        global_log.info("step1_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        fix_side_chain(**global_paths["step1_fixsidechain"], properties=global_prop["step1_fixsidechain"])

        # STEP 2: model SS bonds (CYS -> CYX) if necessary
        global_log.info("step2_fixssbonds: Fix SS bonds")
        fix_ssbonds(**global_paths["step2_fixssbonds"], properties=global_prop["step2_fixssbonds"])

        # STEP 3: Fix amides
        global_log.info("step3_fixamides: fix clashing amides")
        fix_amides(**global_paths["step3_fixamides"], properties=global_prop["step3_fixamides"])

        # STEP 4: Fix chirality
        global_log.info("step4_fixchirality: fix chirality of residues")
        fix_chirality(**global_paths["step4_fixchirality"], properties=global_prop["step4_fixchirality"])

        # STEP 5: renumber structure atoms and residues
        global_log.info("step5_renumberstructure: renumber structure")
        renumber_structure(**global_paths["step5_renumberstructure"], properties=global_prop["step5_renumberstructure"])

        # STEP 6: add H atoms, generate coordinate (.gro) and topology (.top) file
        global_log.info("step6_pdb2gmx: Generate the topology")
        pdb2gmx(**global_paths["step6_pdb2gmx"], properties=global_prop["step6_pdb2gmx"])

        # STEP 7: Create simulation box
        global_log.info("step7_editconf: Create the solvent box")
        editconf(**global_paths["step7_editconf"], properties=global_prop["step7_editconf"])

        # STEP 8: Add solvent molecules
        global_log.info("step8_solvate: Fill the solvent box with water molecules")
        solvate(**global_paths["step8_solvate"], properties=global_prop["step8_solvate"])

        # STEP 9: ion generation pre-processing
        global_log.info("step9_grompp_genion: Preprocess ion generation")
        grompp(**global_paths["step9_grompp_genion"], properties=global_prop["step9_grompp_genion"])

        # STEP 10: ion generation
        global_log.info("step10_genion: Ion generation")
        genion(**global_paths["step10_genion"], properties=global_prop["step10_genion"])

        # Minimization with steepest descend and constraints on the backbone

        # STEP 11: minimization pre-processing
        global_log.info("step11_grompp_min: Preprocess energy minimization 1")
        grompp(**global_paths["step11_grompp_min"], properties=global_prop["step11_grompp_min"])

        # STEP 12: minimization
        global_log.info("step12_mdrun_min: Execute energy minimization 1: steepest descend with backbone constraints")
        mdrun(**global_paths["step12_mdrun_min"], properties=global_prop["step12_mdrun_min"])

        # STEP 13: dump potential energy evolution 
        global_log.info("step13_energy_min: Compute potential energy during minimization")
        gmx_energy(**global_paths["step13_energy_min"], properties=global_prop["step13_energy_min"])

        # Minimization with steepest descend and without constraints on the backbone

        # STEP 14: minimization pre-processing
        global_log.info("step14_grompp_min: Preprocess energy minimization 2")
        grompp(**global_paths["step14_grompp_min"], properties=global_prop["step14_grompp_min"])

        # STEP 15: minimization
        global_log.info("step15_mdrun_min: Execute energy minimization 2: steepest descend without backbone constraints")
        mdrun(**global_paths["step15_mdrun_min"], properties=global_prop["step15_mdrun_min"])

        # STEP 16: dump potential energy evolution
        global_log.info("step16_energy_min: Compute potential energy during minimization")
        gmx_energy(**global_paths["step16_energy_min"], properties=global_prop["step16_energy_min"])

        # STEP 17: convert .gro to .pdb
        global_log.info("step17_trjconv: Convert structure file to pdb format")
        gmx_trjconv_str(**global_paths["step17_trjconv"], properties=global_prop["step17_trjconv"])

        # STEP 18: copy the results of this pose to the merged results folder
        global_log.info("step18_merge_results: Copy the results of this pose to the merged results folder")
        copy_files_with_prefix(global_paths["step18_merge_results"], global_prop["step18_merge_results"]["path"], Path(pose_path).stem)
    
    # Find all .pdb paths in the merged results folder
    pdb_paths = glob.glob(os.path.join(global_paths["step18_merge_results"]["path"], "*.pdb"))

    # Zip the pdb files
    fu.zip_list(zip_file = global_paths["step18_merge_results"]["output_zip_path"], file_list = pdb_paths)

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


if __name__ == "__main__":

    parser = argparse.ArgumentParser("PP docking pose refinement workflow")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    args = parser.parse_args()

    main_wf(configuration_path=args.config_path)
