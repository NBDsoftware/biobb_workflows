#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import sys
import os
import re
import time
import json
import argparse
import shutil
from pathlib import Path, PurePath
from traceback import StackSummary

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_io.api.pdb import pdb
from biobb_io.api.canonical_fasta import canonical_fasta
from biobb_model.model.fix_backbone import fix_backbone
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_altlocs import fix_altlocs
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.checking_log import checking_log
from biobb_model.model.mutate import mutate
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
    
    Inputs:
        *output_paths (str): variable number of paths to output file/s

    Output:
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
        removeFiles(*output_paths)

    return validation_result

def removeFiles(*file_path):
    '''
    Removes files in '*file_path' if they exist

    Inputs
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
        input_file_path  (str): '/path/to/input/input.pdb'
        output_file_path (str): '/path/to/output/output.pdb'
    '''

    # Get '/path/to/output/'
    wdir = PurePath(output_file_path).parent

    # If '/path/to/output/' not created then create
    if not os.path.isdir(wdir):
        os.mkdir(wdir)

    # Copy 'input.pdb' from '/path/to/input/' to '/path/to/output/output.pdb'
    shutil.copy(input_file_path, output_file_path)

def write_pdb_from_gro(output_pdb_path, input_gro_path):
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

    Inputs:
        summary_path   (str): '/path/to/structureCheckSummary.json'

    Output:
        pdb_defects_data (dict): dictionary containing the information for each defect
    
    '''

    # Open check PDB json summary file
    jsonFile = open(summary_path, mode = 'r')

    # Save JSON data as dictionary
    pdb_defects_data = json.load(jsonFile)

    return pdb_defects_data

def checkStructure(step_pdb_path, global_log):
    '''
    Writes PDB check summary in JSON file and writes JSON file data into a python dictionary
    to be used by the workflow.
    
    Inputs:
        step_pdb_path       (str): path to pdb in current step
        global_log    (log class): logger to global log

    Output:
        pdb_defects_data   (dict): dictionary with pdb defects
    '''

    # Identify output folder 
    step_path = PurePath(step_pdb_path).parent

    # Create JSON summary output path
    output = os.path.join(step_path, "pdbStructureCheck.json")

    # Create paths dictionary
    paths = {'input_structure_path': step_pdb_path,
             'output_summary_path': output}
             
    # Create properties dictionary (None: return all features)
    props = {'features': None}

    # Write action to global log and execute step
    global_log.info("    Checking the errors of the PDB structure and creating a JSON summary")
    structure_check(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_summary_path"]):
        global_log.info("    ERROR: PDB Structure check failed.")
        return 0

    # Read information from JSON summary
    pdb_defects_data = readJSONSummary(paths["output_summary_path"])

    return pdb_defects_data

def findCanonicalFASTA(fixBackbone_paths, pdbCode):
    '''
    Finds canonical FASTA from a PDB code and saves it in step path

    Inputs:
        fixBackbone_paths (dict): dictionary with paths of current step
        pdbCode            (str): PDB code 
    '''

    # Define properties NOTE: (api_id?) pdb or pdbe? 
    prop = {'pdb_code': pdbCode}

    # Output path
    output_fasta_path=fixBackbone_paths["input_fasta_canonical_sequence_path"]

    # Identify output folder 
    wdir = PurePath(output_fasta_path).parent

    # If '/path/to/output/' not created then create
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    
    # Fetch canonical FASTA
    canonical_fasta(output_fasta_path=output_fasta_path, properties=prop)

    return

def findBestAltLocs(altloc_data):
    '''
    Finds best alternative location for each residue

    Inputs:
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


def run_wf(args):

    start_time = time.time()

    # Set default value for 'to_do' arg
    if 'to_do' not in args:
        args.to_do = 'free'
        
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.config_path)

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
    
    if 'pdb:' in args.input_pdb_path:

        # If the user gives the pdb ID as 'pdb:id' -> download PDB 

        # Just the id from 'pdb:id'
        pdbCode = args.input_pdb_path.split(':')[1]

        # Update properties dictionary
        props_pdb.update({'pdb_code': pdbCode})
            
        # Write action to global log and execute step
        global_log.info("step1A_pdb: Downloading {} from PDB".format(pdbCode))
        pdb(**paths_pdb, properties=props_pdb)
    
    else:
        # If the user gives the pdb file -> just copy to step folder

        # Write action to global log and execute step
        global_log.info("step1A_pdb: Adding input PDB ({}) to working dir".format(args.input_pdb_path))
        prep_output_file(args.input_pdb_path, paths_pdb["output_pdb_path"])
        
        # Search for pdb code in input.yml
        if props_pdb['pdb_code']:
            pdbCode = props_pdb['pdb_code']
    
    # Validate step
    if not validateStep(paths_pdb["output_pdb_path"]):
        global_log.info("    ERROR: No PDB file was fetched. Check PDB code or input file")
        return 0

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
        return 0

    # Check structure and save results in dictionary
    pdb_defects_data = checkStructure(paths_extract["output_molecule_path"], global_log)

# STEP 2 (A): Fix alternative locations

    # Properties and paths of step
    props_altloc = global_prop["step2A_fixaltlocs"]
    paths_altloc = global_paths["step2A_fixaltlocs"]

    if pdb_defects_data["altloc"]:

        # Find best alternative locations
        best_altlocs = findBestAltLocs(pdb_defects_data["altloc"])

        print(best_altlocs)

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
        return 0

# STEP 2 (B): Add mutations if requested

    # Properties and paths of step
    props_mut = global_prop["step2B_mutations"]
    paths_mut = global_paths["step2B_mutations"]

    if args.mut_list:
        # Update properties dictionary
        props_mut.update({'mutation_list': args.mut_list})
        # Write action to global log and execute step
        global_log.info("step2B_mutations: Preparing mutated structure")
        mutate(**paths_mut, properties=props_mut)

        # Check structure and save results in dictionary
        pdb_defects_data = checkStructure(paths_mut["output_pdb_path"], global_log)
        
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_altloc["output_pdb_path"], paths_mut["output_pdb_path"])

    # Validate step
    if not validateStep(paths_mut["output_pdb_path"]):
        global_log.info("    ERROR: mutation failed.")
        return 0

# STEP 2 (C): model missing heavy atoms of backbone

    # Properties and paths of step
    props_fxBCK = global_prop["step2C_fixbackbone"]
    paths_fxBCK = global_paths["step2C_fixbackbone"]

    # Check if there are missing side chain heavy atoms
    if pdb_defects_data["backbone"] and pdbCode is not None:

        # Use PDB code to find canonical FASTA
        findCanonicalFASTA(paths_fxBCK, pdbCode)

        # Write action to global log and execute step
        global_log.info("step2C_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
        fix_backbone(**paths_fxBCK, properties=props_fxBCK) 

        # Check structure and save results in dictionary
        pdb_defects_data = checkStructure(paths_fxBCK["output_pdb_path"], global_log)
    
    elif pdb_defects_data["backbone"] and pdbCode is None:
        # Write warning to global log
        global_log.info("    WARNING: Missing backbone atoms or breaks but PDB code not provided. Skipping fixbackbone step...")
        # Just copy pdb to step folder  
        prep_output_file(paths_mut["output_pdb_path"], paths_fxBCK["output_pdb_path"])
   
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_mut["output_pdb_path"], paths_fxBCK["output_pdb_path"])

    # Validate step
    if not validateStep(paths_fxBCK["output_pdb_path"]):
        global_log.info("    ERROR: Fixing backbone failed. Check input PDB")
        return 0

# STEP 2 (D): model missing heavy atoms of side chains

    # Properties and paths of step
    props_fxSC = global_prop["step2D_fixsidechain"]
    paths_fxSC = global_paths["step2D_fixsidechain"]

    # Check if there are missing side chain heavy atoms
    if pdb_defects_data["fixside"]:
        # Write action to global log and execute step
        global_log.info("step2D_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        fix_side_chain(**paths_fxSC, properties=props_fxSC)
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_fxBCK["output_pdb_path"], paths_fxSC["output_pdb_path"])

    # Validate step
    if not validateStep(paths_fxSC["output_pdb_path"]):
        global_log.info("    ERROR: Fixing side chains failed. Check input PDB")
        return 0

# STEP 2 (E): model SS bonds (CYS -> CYX) if necessary

    # Properties and paths of step
    props_ss = global_prop["step2E_fixssbonds"]
    paths_ss = global_paths["step2E_fixssbonds"]

    # Check if there are SS bonds
    if pdb_defects_data["getss"]:
        # Write action to global log and execute step
        global_log.info("step2E_fixssbonds: Fix SS bonds")
        fix_ssbonds(**paths_ss, properties=props_ss)
    else:
        # Just copy pdb to step folder  
        prep_output_file(paths_fxSC["output_pdb_path"], paths_ss["output_pdb_path"])

    # Validate step
    if not validateStep(paths_ss["output_pdb_path"]):
        global_log.info("    ERROR: Fixing SS bonds failed.")
        return 0

# STEP 2 (F): renumber structure atoms and residues

    # Properties and paths of step
    props_renum = global_prop["step2F_renumberstructure"]
    paths_renum = global_paths["step2F_renumberstructure"]

    # Write action to global log and execute step
    global_log.info("step2F_renumberstructure: renumber structure")
    renumber_structure(**paths_renum, properties=props_renum)

    # Validate step
    if not validateStep(paths_renum["output_structure_path"]):
        global_log.info("    ERROR: ERROR while renumbering structure")
        return 0

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
        prep_output_file(paths_renum["output_structure_path"], paths_famds["output_pdb_path"])


    # Validate step
    if not validateStep(paths_famds["output_pdb_path"]):
        global_log.info("    ERROR: ERROR while fixing clashing amides")
        return 0

# STEP 2 (H): Final check of the PDB structure and creation of a report for the user

    # Properties and paths of step
    props_check = global_prop["step2H_checkPDB"]
    paths_check = global_paths["step2H_checkPDB"]

    # Write action to global log and execute step
    global_log.info("step2H_checkPDB: check the errors of a PDB structure and create a report log file")
    checking_log(**paths_check, properties=props_check)


    # Check if this should be the final step
    if args.to_do == 'fix':
        shutil.copy(paths_famds["output_pdb_path"], args.output_pdb_path)
        global_log.info("Fix completed. Final structure saved on " + args.output_pdb_path)
        return 0

    # NOTE: ISSUE WARNINGS in log file before simulation?
    

# STEP 3: add H atoms, generate coordinate (.gro) and topology (.top) file

    # Properties and paths of step
    props = global_prop["step3_pdb2gmx"]
    paths = global_paths["step3_pdb2gmx"]

    # Write action to global log and execute step
    global_log.info("step3_pdb2gmx: Generate the topology")
    pdb2gmx(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"], paths["output_top_zip_path"]):
        global_log.info("    ERROR: Coordinates and/or topology where not generated. Check input PDB")
        return 0

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
        return 0

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
        return 0

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
        return 0

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
        return 0

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
        return 0

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
        return 0
    
# STEP 10: dump potential energy evolution

    # Write next action to global log and execute step
    global_log.info("step10_energy_min: Compute potential energy during minimization")
    gmx_energy(**global_paths["step10_energy_min"], properties=global_prop["step10_energy_min"])

    # Check if this should be the final step
    if args.to_do == 'min':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step9_mdrun_min"]["output_gro_path"])
        global_log.info("Minimization completed. Final structure saved on " + args.output_pdb_path)
        return 0

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
        return 0

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
        return 0

# STEP 13: dump temperature evolution

    # Write next action to global log and execute step
    global_log.info("step13_energy_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step13_energy_nvt"], properties=global_prop["step13_energy_nvt"])

    # Check if this should be the final step
    if args.to_do == 'nvt':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step12_mdrun_nvt"]["output_gro_path"])
        global_log.info("NVT Equilibration completed. Final structure saved on " + args.output_pdb_path)
        return 0

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
        return 0

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
        return 0   

# STEP 16: dump density and pressure evolution

    # Write next action to global log and execute step
    global_log.info("step16_energy_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step16_energy_npt"], properties=global_prop["step16_energy_npt"])

    # Check if this should be the final step
    if args.to_do == 'npt':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step15_mdrun_npt"]["output_gro_path"])
        global_log.info("NPT Equilibration completed. Final structure saved on " + args.output_pdb_path)
        return 0

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
        return 0

# STEP 18: free NPT production run

    # Properties and paths of step
    props = global_prop["step18_mdrun_md"]
    paths = global_paths["step18_mdrun_md"]

    # Write next action to global log and execute step
    global_log.info("step18_mdrun_md: Execute free molecular dynamics simulation")
    mdrun(**paths, properties=props)

    # Validate step
    if not validateStep(paths["output_gro_path"]):
        global_log.info("    ERROR: output .gro file was not generated")
        return 0

# STEP 19: dump RMSD with respect to equilibrated structure (first frame)

    # Write next action to global log and execute step
    global_log.info("step19_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")
    gmx_rms(**global_paths["step19_rmsfirst"], properties=global_prop["step19_rmsfirst"])

# STEP 20: dump RMSD with respect to minimized structure

    # Write next action to global log and execute step
    global_log.info("step20_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
    gmx_rms(**global_paths["step20_rmsexp"], properties=global_prop["step20_rmsexp"])

# STEP 21: dump Radius of gyration 

    # Write next action to global log and execute step
    global_log.info("step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
    gmx_rgyr(**global_paths["step21_rgyr"], properties=global_prop["step21_rgyr"])

# STEP 22: image (correct for PBC) the trajectory centering the protein and dumping only protein atoms

    # Write next action to global log and execute step
    global_log.info("step22_image: Imaging the resulting trajectory")
    gmx_image(**global_paths["step22_image"], properties=global_prop["step22_image"])

# STEP 23: remove water and ions from final structure after free NPT production run

    # Write next action to global log and execute step
    global_log.info("step23_dry: Removing water molecules and ions from the resulting structure")
    gmx_trjconv_str(**global_paths["step23_dry"], properties=global_prop["step23_dry"])

    write_pdb_from_gro(args.output_pdb_path, global_paths["step18_mdrun_md"]["output_gro_path"])
    global_log.info("Free Equilibration completed. Final structure saved on " + args.output_pdb_path)
    
    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % args.config_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return 0
    
def main():

    parser = argparse.ArgumentParser("Simple MD Protein Setup")

    parser.add_argument('-i', dest='input_pdb_path',
                        help="Input pdb file or id (as pdb:id)", required=True)

    parser.add_argument('-o', dest='output_pdb_path',
                        help="Output pdb file", required=True)

    # Execute workflow until 'to_do' step -> free executes all steps
    parser.add_argument(
        '--op', dest='to_do', help="Extent of the pipeline to execute (fix, min, nvt, npt, free")

    parser.add_argument('--mut_list', dest='mut_list',
                        help="Mutations list as comma-separated list with the format  'chain_id : Old_residue_code Residue_number New_residue_code'. Examples: 'A:G34T' or 'A:F38C,A:N39W,A:T40G' ")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)")
    
    args = parser.parse_args()

    run_wf(args)


if __name__ == "__main__":
    main()
