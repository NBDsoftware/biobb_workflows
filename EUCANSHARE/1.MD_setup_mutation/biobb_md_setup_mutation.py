#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import os
import time
import json
import argparse
import shutil
from zipfile import ZipFile
from pathlib import Path

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_io.api.pdb import pdb
from biobb_io.api.canonical_fasta import canonical_fasta
from biobb_model.model.fix_backbone import fix_backbone
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_altlocs import fix_altlocs
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.fix_chirality import fix_chirality
from biobb_model.model.checking_log import checking_log
from biobb_model.model.mutate import mutate
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.genion import genion
from biobb_gromacs.gromacs.mdrun import mdrun
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_structure_utils.utils.structure_check import structure_check
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure

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

def removeFiles(*file_path):
    '''
    Removes files in '*file_path' if they exist

    Inputs
    ------

        file_path  (str): variable number of paths to files including filename
    '''
    for path in file_path:
        if os.path.exists(path):
            os.remove(path)
    
    return

def prep_output_file(input_file_path, output_file_path):
    '''
    Copy input_file_path to output_file_path

    Inputs
    ------

        input_file_path  (str): '/path/to/input/input.pdb'
        output_file_path (str): '/path/to/output/output.pdb'
    '''

    # Get '/path/to/output/'
    wdir = Path(output_file_path).parent

    # If '/path/to/output/' not created then create
    if not os.path.isdir(wdir):
        os.mkdir(wdir)

    # Copy 'input.pdb' from '/path/to/input/' to '/path/to/output/output.pdb'
    shutil.copy(input_file_path, output_file_path)

def write_pdb_from_gro(output_pdb_path, input_gro_path):
    '''
    Convert .gro to .pdb file

    Inputs
    ------

        output_pdb_path  (str): '/path/to/input/input.pdb'
        input_gro_path   (str): '/path/to/output/output.pdb'
    '''

    prop = {'selection': 'Protein'}

    gmx_trjconv_str(
        input_structure_path=input_gro_path,
        input_top_path=input_gro_path,
        output_str_path=output_pdb_path,
        properties=prop
    )

def readJSONSummary(summary_path):

    '''
    Reads information from structure_check JSON summary, and returns the contents of the JSON 
    file as a dictionary to be used if necessary by the extra steps

    Inputs
    ------

        summary_path   (str): '/path/to/structureCheckSummary.json'

    Output
    ------

        pdb_defects_data (dict): dictionary containing the information for each defect
    
    '''

    # Open check PDB json summary file
    jsonFile = open(summary_path, mode = 'r')

    # Save JSON data as dictionary
    pdb_defects_data = json.load(jsonFile)

    return pdb_defects_data

def findBestAltLocs(altloc_data):
    '''
    Finds best alternative location for each residue

    Inputs
    ------

        altloc_data (dict): dictionary with information regarding altlocs
    '''

    best_altlocs = []

    for residue in altloc_data.keys():

        res_descriptors = residue.split()

        res_id = res_descriptors[1]

        atomic_data = altloc_data[residue]

        atoms = list(atomic_data.keys())

        occupancy_data = atomic_data[atoms[0]]

        max_occupancy = 0
        chosen_label = "A"

        for altloc in occupancy_data:

            if altloc["occupancy"] > max_occupancy:

                max_occupancy = altloc["occupancy"]
                chosen_label = altloc["loc_label"]
    
        best_altlocs.append(res_id + ":" + chosen_label)

    return best_altlocs

def concatenateTrajs(num_trajs, paths_image_traj, props_image_traj, paths_trjcat, props_trjcat):
    '''
    Creates a zip file with all trajectories and then concatenates them

    Inputs
    ------

        num_trajs         (int): number of trajectories
        paths_image_traj (dict): paths of trajectory imaging step
        props_image_traj (dict): properties of trajectory imaging step
        paths_trjcat     (dict): paths of concatenating step
        props_trjcat     (dict): properties of concatenating step

    '''

    # Generic trajectory name
    output_traj_path = paths_image_traj["output_traj_path"]
    traj_name = Path(output_traj_path).name

    # Path of step imaging trajectories with traj0/ traj1/ traj2/
    img_step_path = props_image_traj["path"]

    # Path for zip file expected by trjcat step
    zip_path = paths_trjcat["input_trj_zip_path"]

    # trjcat step folder
    trjcat_step_path = props_trjcat["path"]

    # If the step folder is not created, create it
    if not os.path.exists(trjcat_step_path):
        os.mkdir(trjcat_step_path)

    # Create a ZipFile object
    zip_file_obj = ZipFile(zip_path, 'w')

    # Write to zip file
    for traj_index in range(num_trajs):

        # Traj folder name
        traj_folder = "traj" + str(traj_index)

        # Generic trajectory folder path
        traj_folder_path = os.path.join(img_step_path, traj_folder)

        # Generic trajectory file path
        traj_file_path = os.path.join(traj_folder_path, traj_name)

        # Add to zip file
        zip_file_obj.write(traj_file_path)
    
    # Close zip file
    zip_file_obj.close()

    # Concatenate trajectories using zip file
    trjcat(**paths_trjcat, properties=props_trjcat)

    return

def addFolderToPaths(all_paths, new_folder, *keywords):
    '''
    Goes over paths corresponding to keywords in all_paths dictionary, adds 
    folder corresponding to current traj. Creates the trajectory and step folders if needed.
    
    Inputs
    ------

        all_paths  (dict): dictionary with old paths of step
        new_folder  (str): string corresponding to new folder name inside old path
        keywords    (str): keywords of output/input paths that will be modified

                            /path/to/step  ---> /path/to/step/new_folder

    Output
    ------

        all_paths  (dict): dictionary with paths with those corresponding to "keywords" modified
    '''

    # For all keys passed in keywords, modify path according to trajname
    for key in keywords:

        # Step path
        stepfolder = Path(all_paths[key]).parent

        # Name of file
        filename = Path(all_paths[key]).name

        # Add traj path to step path
        trajfolder = os.path.join(stepfolder, new_folder)

        # Add file name to traj path
        newpath = os.path.join(trajfolder, filename)

        # Update paths dictionary
        all_paths.update({key : newpath})

    # If the step folder is not created, create it
    if not os.path.exists(stepfolder):
        os.mkdir(stepfolder)

    # If the traj folder is not created, create it
    if not os.path.exists(trajfolder):
        os.mkdir(trajfolder)

    return all_paths


def main_wf(configuration_path, input_pdb, last_step, mutation_list, num_trajs):
    '''
    Main MD Setup, mutation and run workflow. Can be used to retrieve a PDB, fix some defects of the structure, 
    add specific mutations, prepare the system to simulate, minimize it, equilibrate it and finally do N production runs.

    Inputs
    ------

        configuration_path (str): path to input.yml
        input_pdb          (str): either the PDB code in pdb:id format or the path the the pdb file 
        last_step          (str): last step of the workflow to execute ('pdb', 'fix', 'min', 'nvt', 'npt' or 'all')
        mutation_list      (str): Mutations as comma-separated list with the format: 'chain_id : Old_residue_code Residue_number New_residue_code'. 
                                Examples: 'A:G34T' or 'A:F38C,A:N39W,A:T40G'
        num_trajs          (int): number of trajectories to launch one after another

    Outputs
    -------

        /output folder
        structure_path  (str): path to equilibrated imaged structure
        trajectory_path (str): path to imaged trajectory

    '''
    
    start_time = time.time()

    # Set default value for 'last_step' arg
    if last_step is None:
        last_step = 'all'

    # Set default value for 'n_trajs' arg
    if num_trajs is None:
        num_trajs = 1
    else:
        num_trajs = int(num_trajs)
        
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Declaring the steps of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

# STEP 1 (A): Source PDB file

    # Properties and paths of step
    props_pdb = global_prop["step1A_pdb"]
    paths_pdb = global_paths["step1A_pdb"]

    # Initial value for PDB code - change it if user gives PDB code
    pdbCode = None
    
    if 'pdb:' in input_pdb:

        # If the user gives the pdb ID as 'pdb:id' -> download PDB 

        # Just the id from 'pdb:id'
        pdbCode = input_pdb.split(':')[1]

        # Update properties dictionary
        props_pdb.update({'pdb_code': pdbCode})
            
        # Write action to global log and execute step
        global_log.info("step1A_pdb: Downloading {} from PDB".format(pdbCode))
        pdb(**paths_pdb, properties=props_pdb)
    
    else:
        # If the user gives the pdb file -> just copy to step folder

        # Write action to global log and execute step
        global_log.info("step1A_pdb: Adding input PDB ({}) to working dir".format(input_pdb))
        prep_output_file(input_pdb, paths_pdb["output_pdb_path"])
        
        # Search for pdb code in input.yml
        if props_pdb['pdb_code']:
            pdbCode = props_pdb['pdb_code']
    
    # Validate step
    if not validateStep(paths_pdb["output_pdb_path"]):
        global_log.info("    ERROR: No PDB file was fetched. Check PDB code or input file")
        return None, None

# STEP 1 (B): extracts molecule of interest: protein

    # Properties and paths of step
    props_extract = global_prop["step1B_extractMolecule"]
    paths_extract = global_paths["step1B_extractMolecule"]

    # Write action to global log and execute step
    global_log.info("step1B_extractMolecule: extract molecule of interest (protein)")
    extract_molecule(**paths_extract, properties=props_extract)

    # Validate step
    if not validateStep(paths_extract["output_molecule_path"]):
        global_log.info("    ERROR: Extraction of molecule failed.")
        return None, None

    # Properties and paths of step
    props_structureChk = global_prop["step1C_structure_check"]
    paths_structureChk = global_paths["step1C_structure_check"]

    # Write action to global log and execute step
    global_log.info("step1C_structure_check: Checking the errors of the PDB structure and creating a JSON summary")
    structure_check(**paths_structureChk, properties=props_structureChk)

    # Read JSON summary and save it on dictionary
    pdb_defects_data = readJSONSummary(paths_structureChk["output_summary_path"])

    # Check if this should be the final step
    if last_step == 'pdb':
        global_log.info("Molecule extraction completed. Final structure saved on " + paths_extract["output_molecule_path"])
        return paths_extract["output_molecule_path"], None

# STEP 2 (A): Fix alternative locations

    # Properties and paths of step
    props_altloc = global_prop["step2A_fixaltlocs"]
    paths_altloc = global_paths["step2A_fixaltlocs"]

    if pdb_defects_data["altloc"]:

        # Find best alternative locations
        best_altlocs = findBestAltLocs(pdb_defects_data["altloc"])

        # Update properties dictionary
        props_altloc.update({'altlocs': best_altlocs})

        # Write action to global log and execute step
        global_log.info("step2A_fixaltlocs: Fix alternative locations")
        fix_altlocs(**paths_altloc, properties=props_altloc)
        
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_extract["output_molecule_path"], paths_altloc["output_pdb_path"])

    # Validate step
    if not validateStep(paths_altloc["output_pdb_path"]):
        global_log.info("    ERROR: fixing alternative locations failed.")
        return None, None

# STEP 2 (B): Add mutations if requested

    # Properties and paths of step
    props_mut = global_prop["step2B_mutations"]
    paths_mut = global_paths["step2B_mutations"]

    if mutation_list:

        # Update properties dictionary
        props_mut.update({'mutation_list': mutation_list})

        # Write action to global log and execute step
        global_log.info("step2B_mutations: Preparing mutated structure")
        mutate(**paths_mut, properties=props_mut)
        
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_altloc["output_pdb_path"], paths_mut["output_pdb_path"])

    # Validate step
    if not validateStep(paths_mut["output_pdb_path"]):
        global_log.info("    ERROR: mutation failed.")
        return None, None

# STEP 2 (C-D): model missing heavy atoms of backbone NOTE: if break is long only caps, How?

    # Properties and paths of step
    props_fxBCK = global_prop["step2D_fixbackbone"]
    paths_fxBCK = global_paths["step2D_fixbackbone"]

    # Missing bck heavy atoms + known PDB code -> Model
    if pdb_defects_data["backbone"] and pdbCode is not None:

        # Properties and paths of step
        props_fasta = global_prop["step2C_canonical_fasta"]
        paths_fasta = global_paths["step2C_canonical_fasta"]

        # Update properties
        props_fasta.update({'pdb_code' : pdbCode})

        # Use PDB code to find canonical FASTA
        canonical_fasta(**paths_fasta, properties=props_fasta)

        # Write action to global log and execute step
        global_log.info("step2C_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
        fix_backbone(**paths_fxBCK, properties=props_fxBCK) 
    
    # Missing bck heavy atoms + no PDB code -> No Model + Warning
    elif pdb_defects_data["backbone"] and pdbCode is None:

        # Write warning to global log
        global_log.info("    WARNING: Missing backbone atoms or breaks but PDB code not provided. Skipping fixbackbone step...")
        
        # Just copy pdb to step folder  
        prep_output_file(paths_mut["output_pdb_path"], paths_fxBCK["output_pdb_path"])
    
    # No missing bck heavy atoms -> No Model
    else:

        # Just copy pdb to step folder  
        prep_output_file(paths_renum["output_structure_path"], paths_fxBCK["output_pdb_path"])

    # Validate step
    if not validateStep(paths_fxBCK["output_pdb_path"]):
        global_log.info("    ERROR: Fixing backbone failed. Check input PDB")
        return None, None

# STEP 2 (E): model missing heavy atoms of side chains

    # Properties and paths of step
    props_fxSC = global_prop["step2E_fixsidechain"]
    paths_fxSC = global_paths["step2E_fixsidechain"]

    # Check if there are missing side chain heavy atoms
    if pdb_defects_data["fixside"]:

        # Write action to global log and execute step
        global_log.info("step2E_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        fix_side_chain(**paths_fxSC, properties=props_fxSC)
    
    else:

        # Just copy pdb to step folder  
        prep_output_file(paths_fxBCK["output_pdb_path"], paths_fxSC["output_pdb_path"])

    # Validate step
    if not validateStep(paths_fxSC["output_pdb_path"]):
        global_log.info("    ERROR: Fixing side chains failed. Check input PDB")
        return None, None

# STEP 2 (F): model SS bonds (CYS -> CYX) if necessary

    # Properties and paths of step
    props_ss = global_prop["step2F_fixssbonds"]
    paths_ss = global_paths["step2F_fixssbonds"]

    # Check if there are SS bonds
    if pdb_defects_data["getss"]:

        # Write action to global log and execute step
        global_log.info("step2F_fixssbonds: Fix SS bonds")
        fix_ssbonds(**paths_ss, properties=props_ss)

    else:
        
        # Just copy pdb to step folder  
        prep_output_file(paths_fxSC["output_pdb_path"], paths_ss["output_pdb_path"])

    # Validate step
    if not validateStep(paths_ss["output_pdb_path"]):
        global_log.info("    ERROR: Fixing SS bonds failed.")
        return None, None

# STEP 2 (G): Fix amides

    # Properties and paths of step
    props_famds = global_prop["step2G_fixamides"]
    paths_famds = global_paths["step2G_fixamides"]

    if pdb_defects_data["amide"]:

        # Write action to global log and execute step
        global_log.info("step2G_fixamides: fix clashing amides")
        fix_amides(**paths_famds, properties=props_famds)

    else:

        # Just copy pdb to step folder  
        prep_output_file(paths_ss["output_structure_path"], paths_famds["output_pdb_path"])


    # Validate step
    if not validateStep(paths_famds["output_pdb_path"]):
        global_log.info("    ERROR: ERROR while fixing clashing amides")
        return None, None

# STEP 2 (H): Fix chirality

    # Properties and paths of step
    props_fchir = global_prop["step2H_fixchirality"]
    paths_fchir = global_paths["step2H_fixchirality"]

    # NOTE: this is not the same as chirality problems... read from log file?
    if pdb_defects_data["chiral"] or pdb_defects_data["chiral_bck"]:
        
        # Write action to global log and execute step
        global_log.info("step2H_fixchirality: fix chirality of residues")
        fix_chirality(**paths_fchir, properties=props_fchir)

    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_famds["output_pdb_path"], paths_fchir["output_pdb_path"])

    # Validate step
    if not validateStep(paths_fchir["output_pdb_path"]):
        global_log.info("    ERROR: ERROR while fixing chirality of residues")
        return None, None

# STEP 2 (H): renumber structure atoms and residues

    # Properties and paths of step
    props_renum = global_prop["step2I_renumberstructure"]
    paths_renum = global_paths["step2I_renumberstructure"]

    # Write action to global log and execute step
    global_log.info("step2I_renumberstructure: renumber structure")
    renumber_structure(**paths_renum, properties=props_renum)

    # Validate step
    if not validateStep(paths_renum["output_structure_path"]):
        global_log.info("    ERROR: ERROR while renumbering structure")
        return None, None

# STEP 2 (I): Final check of the PDB structure and creation of a report for the user

    # Properties and paths of step
    props_check = global_prop["step2J_checkPDB"]
    paths_check = global_paths["step2J_checkPDB"]

    # Write action to global log and execute step
    global_log.info("step2J_checkPDB: check the errors of a PDB structure and create a report log file")
    checking_log(**paths_check, properties=props_check)

    # Check if this should be the final step
    if last_step == 'fix':
        global_log.info("Fix completed. Final structure saved on " + paths_renum["output_structure_path"])
        return paths_renum["output_structure_path"], None
    
# STEP 3: add H atoms, generate coordinate (.gro) and topology (.top) file

    # Properties and paths of step
    props = global_prop["step3_pdb2gmx"]
    paths = global_paths["step3_pdb2gmx"]

    # Write action to global log and execute step
    global_log.info("step3_pdb2gmx: Generate the topology")
    pdb2gmx(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"], paths["output_top_zip_path"]):
        global_log.info("    ERROR: Coordinates and/or topology were not generated. Check input PDB")
        return None, None

# STEP 4: Create simulation box

    # Properties and paths of step
    props = global_prop["step4_editconf"]
    paths = global_paths["step4_editconf"]

    # Write action to global log and execute step
    global_log.info("step4_editconf: Create the solvent box")
    editconf(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"]):
        global_log.info("    ERROR: Simulation box not created. Check properties of step")
        return None, None

# STEP 5: Add solvent molecules

    # Properties and paths of step
    props = global_prop["step5_solvate"]
    paths = global_paths["step5_solvate"]

    # Write next action to global log and execute step
    global_log.info("step5_solvate: Fill the solvent box with water molecules")
    solvate(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"], paths["output_top_zip_path"]):
        global_log.info("    ERROR: solvate failed to add solvent molecules.")
        return None, None

# STEP 6: ion generation pre-processing

    # Properties and paths of step
    props = global_prop["step6_grompp_genion"]
    paths = global_paths["step6_grompp_genion"]

    # Write next action to global log and execute step
    global_log.info("step6_grompp_genion: Preprocess ion generation")
    grompp(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_tpr_path"]):
        global_log.info("    ERROR: grompp didn't generate the tpr file")
        return None, None

# STEP 7: ion generation 

    # Properties and paths of step
    props = global_prop["step7_genion"]
    paths = global_paths["step7_genion"]

    # Write next action to global log and execute step
    global_log.info("step7_genion: Ion generation")
    genion(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"], paths["output_top_zip_path"]):
        global_log.info("    ERROR: genion didn't generate output files")
        return None, None

# STEP 8: minimization pre-processing

    # Properties and paths of step
    props = global_prop["step8_grompp_min"]
    paths = global_paths["step8_grompp_min"]

    # Write next action to global log and execute step
    global_log.info("step8_grompp_min: Preprocess energy minimization")
    grompp(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_tpr_path"]):
        global_log.info("    ERROR: grompp didn't generate the tpr file")
        return None, None

# STEP 9: minimization

    # Properties and paths of step
    props = global_prop["step9_mdrun_min"]
    paths = global_paths["step9_mdrun_min"]

    # Write next action to global log and execute step
    global_log.info("step9_mdrun_min: Execute energy minimization")
    mdrun(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"]):
        global_log.info("    ERROR: output .gro file was not generated")
        return None, None
    
# STEP 10: dump potential energy evolution NOTE: create plots automatically

    # Write next action to global log and execute step
    global_log.info("step10_energy_min: Compute potential energy during minimization")
    gmx_energy(**global_paths["step10_energy_min"], properties=global_prop["step10_energy_min"])

    # Check if this should be the final step
    if last_step == 'min':
        global_log.info("Minimization completed. Final structure saved on " + paths["output_gro_path"])
        return paths["output_gro_path"], None

# STEP 11: NVT equilibration pre-processing

    # Properties and paths of step
    props = global_prop["step11_grompp_nvt"]
    paths = global_paths["step11_grompp_nvt"]

    # Write next action to global log and execute step
    global_log.info("step11_grompp_nvt: Preprocess system temperature equilibration")
    grompp(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_tpr_path"]):
        global_log.info("    ERROR: grompp didn't generate the tpr file")
        return None, None

# STEP 12: NVT equilibration

    # Properties and paths of step
    props = global_prop["step12_mdrun_nvt"]
    paths = global_paths["step12_mdrun_nvt"]

    # Write next action to global log and execute step
    global_log.info("step12_mdrun_nvt: Execute system temperature equilibration")
    mdrun(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"]):
        global_log.info("    ERROR: output .gro file was not generated")
        return None, None

# STEP 13: dump temperature evolution

    # Write next action to global log and execute step
    global_log.info("step13_energy_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step13_energy_nvt"], properties=global_prop["step13_energy_nvt"])

    # Check if this should be the final step
    if last_step == 'nvt':
        global_log.info("NVT Equilibration completed. Final structure saved on " + paths["output_gro_path"])
        return paths["output_gro_path"], None

# STEP 14: NPT equilibration pre-processing

    # Properties and paths of step
    props = global_prop["step14_grompp_npt"]
    paths = global_paths["step14_grompp_npt"]

    # Write next action to global log and execute step
    global_log.info("step14_grompp_npt: Preprocess system pressure equilibration")
    grompp(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_tpr_path"]):
        global_log.info("    ERROR: grompp didn't generate the tpr file")
        return None, None

# STEP 15: NPT equilibration

    # Properties and paths of step
    props = global_prop["step15_mdrun_npt"]
    paths = global_paths["step15_mdrun_npt"]

    # Write next action to global log and execute step
    global_log.info("step15_mdrun_npt: Execute system pressure equilibration")
    mdrun(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"]):
        global_log.info("    ERROR: output .gro file was not generated")
        return None, None

# STEP 16: dump density and pressure evolution

    # Write next action to global log and execute step
    global_log.info("step16_energy_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step16_energy_npt"], properties=global_prop["step16_energy_npt"])

    # Check if this should be the final step
    if last_step == 'npt':
        global_log.info("NPT Equilibration completed. Final structure saved on " + paths["output_gro_path"])
        return paths["output_gro_path"], None

# STEP 17: free NPT production run pre-processing

    # Properties and paths of step
    props = global_prop["step17_grompp_md"]
    paths = global_paths["step17_grompp_md"]

    # Write next action to global log and execute step
    global_log.info("step17_grompp_md: Preprocess free dynamics")
    grompp(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_tpr_path"]):
        global_log.info("    ERROR: grompp didn't generate the tpr file")
        return None, None

# STEP 18: free NPT production run

    # Write next action to global log 
    global_log.info("step18_mdrun_md: Execute free molecular dynamics simulation")
    global_log.info("     Number of trajectories: {}".format(num_trajs))
    
    for traj_index in range(num_trajs):

        # Properties and paths of step
        props_mdrun = global_prop["step18_mdrun_md"].copy()
        paths_mdrun = global_paths["step18_mdrun_md"].copy()

        # Modify paths according to traj
        addFolderToPaths(paths_mdrun,
                        "traj" + str(traj_index),
                        "output_trr_path",
                        "output_gro_path",
                        "output_edr_path",
                        "output_log_path",
                        "output_cpt_path")
        
        # Execute step
        mdrun(**paths_mdrun, properties=props_mdrun)

        # Validate step NOTE: different processes accessing same global log, if we were to parallelize this?
        if not validateStep(paths_mdrun["output_gro_path"]):
            global_log.info("    ERROR: output .gro file was not generated")
            return None, None

    # STEP 19: dump RMSD with respect to equilibrated structure (first frame) NOTE: add computation of RMSF, PCA and projection onto PC for visualization

        # Write next action to global log 
        global_log.info("step19_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")

        # Properties and paths of step
        props_rmsFst = global_prop["step19_rmsfirst"].copy()
        paths_rmsFst = global_paths["step19_rmsfirst"].copy()

        # Modify paths according to traj
        addFolderToPaths(paths_rmsFst,
                        "traj" + str(traj_index),
                        "input_traj_path",
                        "output_xvg_path")

        # Execute step
        gmx_rms(**paths_rmsFst, properties=props_rmsFst)

    # STEP 20: dump RMSD with respect to minimized structure

        # Write next action to global log 
        global_log.info("step20_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
        
        # Properties and paths of step
        props_rmsExp = global_prop["step20_rmsexp"].copy()
        paths_rmsExp = global_paths["step20_rmsexp"].copy()

        # Modify paths according to traj
        addFolderToPaths(paths_rmsExp,
                        "traj" + str(traj_index),
                        "input_traj_path",
                        "output_xvg_path")
        
        # Execute step
        gmx_rms(**paths_rmsExp, properties=props_rmsExp)

    # STEP 21: dump Radius of gyration NOTE: add computation of RMSF, PCA and projection onto PC for visualization

        # Write next action to global log 
        global_log.info("step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")

        # Properties and paths of step
        props_rgyr = global_prop["step21_rgyr"].copy()
        paths_rgyr = global_paths["step21_rgyr"].copy()

        # Modify paths according to traj
        addFolderToPaths(paths_rgyr,
                        "traj" + str(traj_index),
                        "input_traj_path",
                        "output_xvg_path")

        # Execute step
        gmx_rgyr(**paths_rgyr, properties=props_rgyr)

    # STEP 22: image (correct for PBC) the trajectory centering the protein and dumping only protein atoms

        # Write next action to global log 
        global_log.info("step22_image: Imaging the resulting trajectory")

        # Properties and paths of step
        props_img = global_prop["step22_image"].copy()
        paths_img = global_paths["step22_image"].copy()

        # Modify paths according to traj
        addFolderToPaths(paths_img,
                        "traj" + str(traj_index),
                        "input_traj_path",
                        "output_traj_path")

        # Execute step
        gmx_image(**paths_img, properties=props_img)

# STEP 23: remove water and ions from structure obtained after equilibration, before production run

    # Write next action to global log 
    global_log.info("step23_dry: Removing water molecules and ions from the equilibrated structure")
        
    # Execute step
    gmx_trjconv_str(**global_paths["step23_dry"], properties=global_prop["step23_dry"])

# STEP 24: concatenate trajectories

    # Write next action to global log 
    global_log.info("step24_trjcat: Concatenate trajectories")

    # Properties and paths of step
    props_trjcat = global_prop["step24_trjcat"].copy()
    paths_trjcat = global_paths["step24_trjcat"].copy()

    # Properties and paths of step
    paths_mdrun = global_paths["step18_mdrun_md"].copy()

    concatenateTrajs(num_trajs = num_trajs, paths_image_traj = paths_img, 
                    props_image_traj=props_img, paths_trjcat = paths_trjcat, 
                    props_trjcat = props_trjcat)

    # Validate step
    if not validateStep(paths_trjcat["output_trj_path"]):
        global_log.info("    ERROR: Concatenation of trajectories failed.")
        return None, None

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

    structure_path = global_paths["step23_dry"]["output_str_path"]
    trajectory_path = paths_trjcat["output_trj_path"]

    return structure_path, trajectory_path

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Simple MD Protein Setup")

    parser.add_argument('--input', dest='input',
                        help="Input pdb file or pdb id (as pdb:id)", 
                        required=False)

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)

    # Execute workflow until 'last_step' -> all executes all steps (all is default)
    parser.add_argument('--until', dest='last_step', 
                        help="Extent of the pipeline to execute (pdb, fix, min, nvt, npt or all; default: all)", 
                        required=False)

    parser.add_argument('--mut_list', dest='mut_list',
                        help="Mutations as comma-separated list with the format: 'chain_id : Old_residue_code Residue_number New_residue_code'. Examples: 'A:G34T' or 'A:F38C,A:N39W,A:T40G' (default: None)", 
                        required=False)

    parser.add_argument('--n_trajs', dest='n_trajs',
                        help="Number of trajectories (default: 1)", 
                        required=False)    
    
    args = parser.parse_args()

    # Define global constants
    global HOMOLOGY_MAX_RES

    # NOTE: feature not implemented yet
    # Maximum number of residues to be created using homology modelling
    HOMOLOGY_MAX_RES = 10

    _,_ = main_wf(configuration_path = args.config_path, input_pdb = args.input, last_step = args.last_step, 
            mutation_list = args.mut_list, num_trajs = args.n_trajs)
