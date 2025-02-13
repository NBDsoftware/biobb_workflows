#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
from Bio.SeqIO.PdbIO import PdbSeqresIterator
from Bio.PDB import PDBParser
from Bio import SeqIO
from typing import List, Dict, Union
from pathlib import Path
import time
import random
import argparse
import os

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_io.api.canonical_fasta import canonical_fasta
from biobb_model.model.fix_backbone import fix_backbone
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_altlocs import fix_altlocs
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.fix_chirality import fix_chirality
from biobb_model.model.mutate import mutate

from biobb_amber.leap.leap_gen_top import leap_gen_top
from biobb_amber.sander.sander_mdrun import sander_mdrun
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
from biobb_amber.process.process_minout import process_minout

from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf
from biobb_structure_utils.utils.cat_pdb import cat_pdb
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.extract_heteroatoms import extract_heteroatoms
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_pdb_tools.pdb_tools.biobb_pdb_tofasta import biobb_pdb_tofasta


# Biopython helpers
def highest_occupancy_altlocs(pdb_file, global_log) -> List[str]:
    """
    Reads a PDB file and returns a list of the highest occupancy alternative locations
    for each residue that has multiple conformations (alternative locations). 

    The output is a list where each element is a string in the format: 
    "<chain_id><residue_number>:<altLoc>", representing the chain, residue number, 
    and the alternative location identifier with the highest occupancy.

    Args:
        pdb_file (str): Path to the PDB file to be parsed.

    Returns:
        List[str]: A list of strings where each string represents the residue's chain ID, 
                   residue number, and the highest occupancy alternative location identifier.
                   The format of each string is "<chain_id><residue_number>:<altLoc>".

    Example:
        If a residue with ID 339 in chain 'A' has two alternative locations 'A' and 'B',
        and 'A' has a higher occupancy, the output for this residue would be "A339:A".

        >>> highest_occupancy_altlocs('example.pdb')
        ["A339:A", "A171:B", "A768:A"]
    """
    
    parser = PDBParser(QUIET=True)
    
    # Check if the file exists
    if not os.path.exists(pdb_file):
        global_log.error(f"File {pdb_file} not found")
        return []
    
    structure = parser.get_structure('structure', pdb_file)

    altloc_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                
                altloc_dict = {}
                for atom in residue:
                    altloc = atom.get_altloc()
                    
                    if altloc and altloc != " ":  # Check if the atom has an alternative location
                        # Keep track of highest occupancy for each altloc
                        if altloc not in altloc_dict or atom.get_occupancy() > altloc_dict[altloc]['occupancy']:
                            altloc_dict[altloc] = {'occupancy': atom.get_occupancy(), 'residue': residue}
                            
                # If there is any alternative location
                if altloc_dict:
                    # Find the altloc with the highest occupancy
                    best_altloc = max(altloc_dict, key=lambda x: altloc_dict[x]['occupancy'])
                    res_id = residue.get_id()[1]
                    chain_id = chain.get_id()
                    altloc_residues.append(f"{chain_id}{res_id}:{best_altloc}")
    
    # Log result
    if altloc_residues:
        global_log.info(f"Found residues with alternative locations: {altloc_residues}")
    else:
        global_log.info("No residues with alternative locations found")
        
    return altloc_residues

def get_ligands(ligands_top_folder: Union[str, None], global_log) -> List[Dict[str, str]]:
    """
    Get a list of available ligands in the ligands topology folder. The function searches for all the .frcmod and .lib / .prep files in the folder.
    
    If the ligands folder is provided but doesn't exist or any of the ligands is missing a file an error is raised. 
    
    If the ligand folder is None, an empty list is returned.
    
    Inputs
    ------
    
        ligands_top_folder (str): Path to the folder with the ligand .frcmod and .lib / .prep files.
        global_log: Logger object for logging messages.
    
    Returns
    -------
    
        ligands: Dictionary with the ligand names, topology and coordinate file paths. Empty dict if no ligands are found. 
        
            Example 
                    ligands = {
                        'ZZ7': {
                            'force_modification': 'path/to/ZZ7.frcmod',
                            'library': 'path/to/ZZ7.lib'
                        },
                        'ZZ8': {
                            'force_modification': 'path/to/ZZ8.frcmod',
                            'preparation': 'path/to/ZZ8.prep'
                        }
                    }
    """
    
    # Initialize the list of ligands
    ligands = {}

    # If ligands folder is not provided, return an empty list
    if ligands_top_folder is None:
        return ligands
    
    # Check if the ligands folder exists
    if not os.path.exists(ligands_top_folder):
        global_log.error(f"Folder {ligands_top_folder} not found")
        return ligands
    
    # Search for ligands in the ligands folder
    for file in os.listdir(ligands_top_folder):
        
        # Check if the file is a .itp or .gro file
        if file.endswith(".frcmod") or file.endswith(".lib") or file.endswith(".prep"): 
            
            # Get the file name without extension
            ligand_id = Path(file).stem   
            
            # Get the file extension
            file_extension = Path(file).suffix
            
            # Check if the ligand name is already in the dictionary
            if ligand_id in ligands:
                if file_extension == ".frcmod":
                    ligands[ligand_id]['force_modification'] = os.path.join(ligands_top_folder, file)
                elif file_extension == ".lib":
                    ligands[ligand_id]['library'] = os.path.join(ligands_top_folder, file)
                elif file_extension == ".prep":
                    ligands[ligand_id]['preparation'] = os.path.join(ligands_top_folder, file)
            else:
                if file_extension == ".frcmod":
                    ligands[ligand_id] = {'force_modification': os.path.join(ligands_top_folder, file)}
                elif file_extension == ".lib":
                    ligands[ligand_id] = {'library': os.path.join(ligands_top_folder, file)}
                elif file_extension == ".prep":
                    ligands[ligand_id] = {'preparation': os.path.join(ligands_top_folder, file)}
    
    # Check if all ligands have both a force modification and a library / preparation file
    for ligand, files in ligands.items():
        if 'force_modification' not in files:
            global_log.error(f"Frcmod file for ligand {ligand} not found")
            
            # Remove ligand from the dictionary
            ligands.pop(ligand)
            continue
        
        if 'library' not in files and 'preparation' not in files:
            global_log.error(f"Library or Prep file for ligand {ligand} not found")
            
            # Remove ligand from the dictionary
            ligands.pop(ligand)
            continue
        
    return ligands
      
def fasta_from_pdb(input_pdb_path: str, output_fasta_path: str, global_log) -> bool:
    """
    Try to obtain the FASTA sequence using the SEQRES records in the PDB file with Biopython. If the SEQRES records are available, write the FASTA sequence to 
    the output file and return True. If the SEQRES records are not available, return False.
    
    Inputs
    ------
    
        input_pdb_path (str): Path to the input PDB file.
        output_fasta_path (str): Path to the output FASTA file.
    
    Returns
    -------
    
        bool: Whether the FASTA sequence was obtained from SEQRES records or not.
    """
    
    # Open the PDB file and use the PdbSeqresIterator to extract sequences
    with open(input_pdb_path, "r") as handle:
        sequences = list(PdbSeqresIterator(handle))  # Extract all sequences into a list
        
        if not sequences:
            global_log.warning(f"PDB doesn't contain SEQRES records")
            return False
        
        global_log.info(f"PDB does contain SEQRES records {sequences}")
        
        # Find parent folder of output_fasta_path
        parent_folder = os.path.dirname(output_fasta_path)
        
        # Create parent folder if it does not exist
        if not os.path.exists(parent_folder):
            os.makedirs(parent_folder)
        
        # Write sequences to a FASTA file
        with open(output_fasta_path, "w") as fasta_out:
            SeqIO.write(sequences, fasta_out, "fasta")
    
    return True

def get_pdb_code(pdb_file: str) -> Union[str, None]:
    """
    Retrieve the PDB code from a PDB file, if available.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    str or None:
        The PDB code if available, otherwise None.
    """
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Access the header information
    header = structure.header

    # Retrieve the PDB code (idcode)
    pdb_code = header.get('idcode', None)

    return pdb_code

def find_amber_his(pdb_file: str, global_log) -> List[str]: # NOTE: This should not be needed!
    """
    Find all histidine residues in the PDB file (AMBER naming convention).
    Histidine names after pdb4amber should be: HIE (epsilon protonated), HID (delta protonated), 
    HIP (both protons present).
    
    Parameters
    ----------
    
    pdb_file : str
        Path to the PDB file.
    
    global_log : Logger
        Logger object for logging messages.
    
    Returns 
    -------
    
    List[str]:
        A list of histidine residues in the PDB file.
    """
    
    his_names = ['HIE', 'HID', 'HIP']
    
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    # Initialize the list of histidine residues
    histidine_residues = {}
    
    # Iterate over all residues in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                
                # Find the residue name
                res_name = residue.get_resname()
                
                # Find the residue number
                res_num = residue.get_id()[1]
                
                # Check if the residue is a histidine
                if res_name in his_names:
                    histidine_residues[res_num] = res_name
                    
    # Log the result
    if histidine_residues:
        global_log.info(f"Found histidine residues: {histidine_residues}")
    else:
        global_log.info("No histidine residues found")

    return list(histidine_residues.values())

def get_pdb2gmx_his(his_residues: List[str]) -> str: # NOTE: This should not be needed!
    """
    Transform the histidine residue names from AMBER to a string
    that can be used as input for the -his flag in the pdb2gmx tool.
    
        HID -> 0, HIE -> 1, HIP -> 2 
        
    Parameters
    ----------
    
    his_residues : List[str]
        List of histidine residues in the PDB file.
    
    Returns
    -------
    
    str:
        A string with the histidine residues transformed to the pdb2gmx format.
    """

    his_dict = {'HID': '0', 'HIE': '1', 'HIP': '2'}
    
    his_pdb2gmx = [his_dict[his] for his in his_residues]
    
    return " ".join(his_pdb2gmx)

# Set additional general properties not considered by the configuration reader
def set_amber_path(global_properties: dict, binary_path: str) -> None: # NOTE: Adapt this
    """
    Set the path to the AMBER binary for all steps using AMBER.

    Inputs
    ------

        global_properties (dict): Dictionary containing the global_properties.
        binary_path (str): Path to the AMBER binary.
    """

    list_of_steps = ['step3B_structure_topology', 'step3K_editconf', 'step3L_solvate', 'step3M_grompp_genion', 'step3N_genion',
                        'step4A_grompp_min', 'step4B_mdrun_min', 'step4E_grompp_nvt', 'step4F_mdrun_nvt',
                        'step4H_grompp_npt', 'step4I_mdrun_npt', 'step5A_grompp_md', 'step5B_mdrun_md']

    for step in list_of_steps:
        global_properties[step]['binary_path'] = binary_path

def set_mpi_path(global_properties: dict, mpi_bin: str, mpi_np: int) -> None: # NOTE: Adapt this
    """
    Set the path to the MPI binary for all steps using MPI.

    Inputs
    ------

        global_properties (dict): Dictionary containing the global_properties.
        mpi_bin (str): Path to the MPI binary.
        mpi_np (int): Number of processors to be used.
    """

    list_of_steps = ['step4B_mdrun_min', 'step4F_mdrun_nvt', 'step4I_mdrun_npt', 'step5B_mdrun_md']

    for step in list_of_steps:
        global_properties[step]['mpi_bin'] = mpi_bin
        global_properties[step]['mpi_np'] = mpi_np
def set_gpu_use(global_properties: dict, gpu_use: bool) -> None: # NOTE: Adapt this
    """
    Set the use of GPU for all steps using AMBER that support it.

    Inputs
    ------

        global_properties (dict): Dictionary containing the global_properties.
        gpu_use (bool): Whether to use GPU or not.
    """

    list_of_steps = ['step4F_mdrun_nvt', 'step4I_mdrun_npt', 'step5B_mdrun_md']

    for step in list_of_steps:
        global_properties[step]['use_gpu'] = gpu_use
        
def set_global_amber_properties(global_properties: dict, amber_properties: dict, global_log) -> None: # NOTE: Adapt this
    """
    Set all the gmx global properties of this workflow, i.e. those global properties included at the beginning of the YAML configuration file that
    are general to some gmx steps.
    
    Inputs
    ------
    
        global_properties (dict): Dictionary containing the global_properties.
        amber_properties (dict): Dictionary containing the gmx properties.
        global_log (Logger): Logger object for logging messages.
    """
    
    # Enforce amber binary path for all steps using it
    if amber_properties.get('binary_path'):
        global_log.info(f"Using AMBER binary path: {amber_properties['binary_path']}")
        set_amber_path(global_properties, amber_properties['binary_path'])

    # Enforce mpi binary path for all steps using mpi
    if amber_properties.get('mpi_bin'):
        global_log.info(f"Using MPI binary path: {amber_properties['mpi_bin']}")
        set_mpi_path(global_properties, amber_properties['mpi_bin'], amber_properties.get('mpi_np'))

    # Enforce gpu use for all steps using amber that support it
    if amber_properties.get('use_gpu'):
        global_log.info(f"Using GPU for AMBER steps")
        set_gpu_use(global_properties, amber_properties['use_gpu'])
        
# Process input/output files
def zip_ligand_topologies(ligands_dict: dict, global_paths: dict, global_prop: dict) -> None:
    """
    Zip all ligand input files (library, preparation, force modification) and update the global paths.
    
    Inputs
    ------
    
        ligands_dict (dict): Dictionary with the ligand names, topology and coordinate file paths.
        global_paths (dict): Dictionary containing the global_paths.
        global_prop (dict): Dictionary containing the global_properties.
    """
    # Zip all ligand inputs
    lib_files = [ligands_dict[ligand].get('library') for ligand in ligands_dict if ligands_dict[ligand].get('library') is not None]
    prep_files = [ligands_dict[ligand].get('preparation') for ligand in ligands_dict if ligands_dict[ligand].get('preparation') is not None]
    frcmod_files = [ligands_dict[ligand].get('force_modification') for ligand in ligands_dict if ligands_dict[ligand].get('force_modification') is not None]
    
    if lib_files:
        input_lib_zip = os.path.join(global_prop["step3B_complex_topology"]["path"], "ligands_lib.zip")
        fu.zip_list(input_lib_zip, lib_files)
        global_paths["step3B_complex_topology"]["input_lib_path"] = input_lib_zip
    if prep_files:
        input_prep_zip = os.path.join(global_prop["step3B_complex_topology"]["path"], "ligands_prep.zip")
        fu.zip_list(input_prep_zip, prep_files)
        global_paths["step3B_complex_topology"]["input_prep_path"] = input_prep_zip
    if frcmod_files:
        input_frcmod_zip = os.path.join(global_prop["step3B_complex_topology"]["path"], "ligands_frcmod.zip")
        fu.zip_list(input_frcmod_zip, frcmod_files)
        global_paths["step3B_complex_topology"]["input_frcmod_path"] = input_frcmod_zip
       
            
def main_wf(configuration_path, input_pdb_path = None, pdb_code = None, pdb_chains = None, mutation_list = None, 
            ligands_top_folder = None, skip_fix_backbone = None, skip_fix_side_chain = None, 
            fix_ss = None, fix_amide_clashes = None, forcefields = ['protein.ff14SB','DNA.bsc1','gaff'], setup_only = False, 
            incrd_path = None, prmtop_path = None, equil_only = False, nsteps = None, num_parts = 1, num_replicas = 1, 
            final_analysis = None, output_path = None):
    '''
    Main setup, mutation and MD run workflow with AMBER. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system, minimize it, equilibrate it and finally do N production runs (either replicas or parts).

    Inputs
    ------

        configuration_path   (str): path to YAML configuration file
        input_pdb_path       (str): (Optional) path to input PDB file
        pdb_code             (str): (Optional) PDB code to be used to get the canonical FASTA sequence
        pdb_chains           (str): (Optional) list of chains to be extracted from the PDB file and fixed
        mutation_list        (str): (Optional) list of mutations to be introduced in the structure
        ligands_top_folder   (str): (Optional) path to the folder containing the ligand .itp and .gro files
        skip_fix_backbone   (bool): (Optional) whether to skip the fix of the backbone atoms
        skip_fix_side_chain (bool): (Optional) whether to skip the fix of the side chain atoms
        fix_ss              (bool): (Optional) wether to add disulfide bonds
        fix_amide_clashes   (bool): (Optional) wether to flip clashing amides to relieve the clashes
        forcefields          (str): (Optional) forcefield to be used in the simulation. Default: ['protein.ff14SB','DNA.bsc1','gaff']  
        setup_only          (bool): (Optional) whether to only setup the system or also run the simulations
        incrd_path           (str): (Optional) path to already-prepared input coordinates file (.inpcrd)
        prmtop_path          (str): (Optional) path to already-prepared input topology file (.prmtop)
        equil_only          (bool): (Optional) whether to only run the equilibration or also run the production simulations
        nsteps               (int): (Optional) Total number of steps of the production simulation
        num_parts            (int): (Optional) number of parts of the trajectory 
        num_replicas         (int): (Optional) number of replicas of the trajectory
        final_analysis      (bool): (Optional) whether to perform the final analysis or not
        output_path          (str): (Optional) path to output folder

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    ###########################
    # Workflow initialization #
    ###########################
    
    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Enforce output_path if provided
    if output_path is not None:
        output_path = fu.get_working_dir_path(output_path, restart = conf.properties.get('restart', 'False'))
        conf.working_dir_path = output_path
    else:
        output_path = conf.get_working_dir_path()
        
    # Remove amber-specific properties from global properties - otherwise they will be applied to all steps
    amber_properties = conf.global_properties.pop('amber', None)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    
    # Check num parts and num of replicas
    if (num_replicas is None) and (num_parts is None):
        # Default behavior: 1 replica
        num_replicas=1
        global_log.info("Number of replicas not set. Defaulting to 1")
    elif (num_replicas is not None) and (num_parts is not None):
        # Cannot set both num parts and num replicas
        global_log.error("Number of trajectories and replicas cannot be set at the same time")

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Set properties for gmx steps
    # set_global_amber_properties(global_prop, amber_properties, global_log) # NOTE: ask for the amber installation path - then complete to the different programs
    ligands_dict = get_ligands(ligands_top_folder, global_log)

    ##############################################
    # Extract atoms and prepare structure for MD #
    ##############################################
    
    # If prepared structure is not provided
    if incrd_path is None:

        # If input PDB is given as argument
        if input_pdb_path is not None:
            global_paths["step1A_extractAtoms"]["input_structure_path"] = input_pdb_path
            global_paths["step1B_extractLigands"]["input_structure_path"] = input_pdb_path
            
        # If chains are given as argument
        if pdb_chains is not None:
            global_prop["step1A_extractAtoms"]["molecule_type"] = "chains"
            global_prop["step1A_extractAtoms"]["chains"] = pdb_chains
        
        # STEP 1A: extract main structure of interest while removing water and ligands (heteroatoms)
        global_log.info("step1A_extractAtoms: extract chain of interest")
        extract_molecule(**global_paths["step1A_extractAtoms"], properties=global_prop["step1A_extractAtoms"])
        
        # Define the heteroatoms to be extracted based on the ligands in the topology folder NOTE: What happens if we try to extract a non-existing ligand? Avoid an error and issue a warning
        global_prop["step1B_extractLigands"]["heteroatoms"] = [{'name': ligand_name} for ligand_name in ligands_dict.keys()]
        
        # STEP 1B: extract ligands
        global_log.info("step1B_extractLigands: extract ligands")
        extract_heteroatoms(**global_paths["step1B_extractLigands"], properties=global_prop["step1B_extractLigands"])

        ###########################################
        # Prepare topology and coordinates for MD #
        ###########################################
        
        # STEP 2A: Fix alternative locations
        global_log.info("step2A_fixaltlocs: Fix alternative locations")
        global_prop["step2A_fixaltlocs"]["altlocs"] = highest_occupancy_altlocs(global_paths["step1A_extractAtoms"]["input_structure_path"], global_log)
        fix_altlocs(**global_paths["step2A_fixaltlocs"], properties=global_prop["step2A_fixaltlocs"])

        # STEP 2B: Add mutations if requested
        if mutation_list is not None:
            global_prop["step2B_mutations"]["mutation_list"] = ",".join(mutation_list)
        global_log.info("step2B_mutations: Preparing mutated structure")
        mutate(**global_paths["step2B_mutations"], properties=global_prop["step2B_mutations"])
        
        # Model the backbone atoms
        if not skip_fix_backbone:
            
            # STEP 2C: Try to get the FASTA sequence to model the backbone from ...
            fasta_available = False
            
            # ... an http request to the PDB
            try:
                global_log.info("step2C_canonical_fasta: Get canonical FASTA")
                file_pdb_code = get_pdb_code(global_paths["step1A_extractAtoms"]["input_structure_path"])
                
                if None not in [pdb_code, file_pdb_code]:
                    # Make sure the PDB code is the same as the one in the input PDB file
                    if pdb_code != file_pdb_code:
                        global_log.warning(f"step2C_canonical_fasta: Provided PDB code ({pdb_code}) is different from the one in the input PDB file ({file_pdb_code}).")
                        global_log.warning(f"step2C_canonical_fasta: Using the provided PDB code ({pdb_code}).")

                if pdb_code is None:
                    pdb_code = file_pdb_code
                
                if pdb_code is not None:
                    global_prop["step2C_canonical_fasta"]["pdb_code"] = pdb_code
                    canonical_fasta(**global_paths["step2C_canonical_fasta"], properties=global_prop["step2C_canonical_fasta"])
                    fasta_available = True
            except:
                global_log.warning("step2C_canonical_fasta: Could not get canonical FASTA. Check the internet connection in the machine running the workflow. Trying to get the canonical FASTA from the PDB file...")
                fasta_available = False
            
            # ... from SEQRES records in the PDB file
            if not fasta_available:
                global_log.info("step2C_pdb_tofasta: Get FASTA from SEQRES of PDB file")
                fasta_available = fasta_from_pdb(global_paths["step1A_extractAtoms"]["input_structure_path"], global_paths["step2C_pdb_tofasta"]["output_file_path"], global_log)

                # Update fix backbone input
                global_paths['step2D_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step2C_pdb_tofasta']['output_file_path']
                
            # ... from the residues in the PDB file (not canonical)
            if not fasta_available:
                global_log.info("step2C_pdb_tofasta: Get FASTA from PDB file")
                
                # Update the input file path
                global_paths['step2C_pdb_tofasta']['input_file_path'] = global_paths["step1A_extractAtoms"]["input_structure_path"]
                
                # Only existing residues in the PDB file are included
                biobb_pdb_tofasta(**global_paths["step2C_pdb_tofasta"], properties=global_prop["step2C_pdb_tofasta"])
                
                # Update fix backbone input
                global_paths['step2D_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step2C_pdb_tofasta']['output_file_path']
                fasta_available = True
                
            # STEP 2D: Model missing heavy atoms of backbone
            if fasta_available:
                global_log.info("step2D_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
                fix_backbone(**global_paths["step2D_fixbackbone"], properties=global_prop["step2D_fixbackbone"])
            else:
                global_log.warning("step2D_fixbackbone: Could not get FASTA sequence. Skipping modeling of the missing heavy atoms in the backbone.")
                global_paths['step2E_fixsidechain']['input_pdb_path'] = global_paths['step2B_mutations']['output_pdb_path']
        else:
            global_log.info("step2D_fixbackbone: Skipping modeling of the missing heavy atoms in the backbone")
            global_paths['step2E_fixsidechain']['input_pdb_path'] = global_paths['step2B_mutations']['output_pdb_path']

        # STEP 2E: model missing heavy atoms of side chains
        if not skip_fix_side_chain:
            global_log.info("step2E_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
            fix_side_chain(**global_paths["step2E_fixsidechain"], properties=global_prop["step2E_fixsidechain"])
        else:
            global_log.info("step2E_fixsidechain: Skipping modeling of the missing heavy atoms in the side chains")
            if not skip_fix_backbone:
                global_paths['step2F_fixssbonds']['input_pdb_path'] = global_paths['step2D_fixbackbone']['output_pdb_path']
            else:
                global_paths['step2F_fixssbonds']['input_pdb_path'] = global_paths['step2B_mutations']['output_pdb_path']

        # STEP 2F: model SS bonds (CYS -> CYX)
        if fix_ss:
            global_log.info("step2F_fixssbonds: Fix SS bonds")
            fix_ssbonds(**global_paths["step2F_fixssbonds"], properties=global_prop["step2F_fixssbonds"])
        else:
            global_paths['step2G_fixamides']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

        # STEP 2G: Rotate amide groups to fix clashes
        if fix_amide_clashes:
            global_log.info("step2G_fixamides: fix clashing amides")
            fix_amides(**global_paths["step2G_fixamides"], properties=global_prop["step2G_fixamides"])
        else:
            if fix_ss:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2F_fixssbonds']['output_pdb_path']
            else:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

        # STEP 2H: Fix chirality
        global_log.info("step2H_fixchirality: fix chirality of residues")
        fix_chirality(**global_paths["step2H_fixchirality"], properties=global_prop["step2H_fixchirality"])

        # STEP 2I: Concatenate all PDB files
        global_log.info("step2I_catpdb: Concatenate ligands and protein PDB files")
        cat_pdb(**global_paths["step2I_catpdb"], properties=global_prop["step2I_catpdb"])
         
        # STEP 2J: Renumber structure atoms and residues
        global_log.info("step2J_renumberstructure: renumber structure")
        renumber_structure(**global_paths["step2J_renumberstructure"], properties=global_prop["step2J_renumberstructure"])
        
        ###########################################
        # Prepare topology and coordinates for MD #
        ###########################################
        
        # NOTE: Take into account the input PDB might contain an RNA/DNA chain
        
        # STEP 3A: Clean PDB for amber
        global_log.info("step3A_pdb4amber: Clean PDB for AMBER")
        pdb4amber_run(**global_paths["step3A_pdb4amber"], properties=global_prop["step3A_pdb4amber"])
        
        # STEP 3B: Create topology and coordinate with tleap
        global_log.info("step3B_complex_topology: Create topology and coordinate files")
        global_prop["step3B_complex_topology"]["forcefield"] = forcefields
        zip_ligand_topologies(ligands_dict, global_paths, global_prop)
        leap_gen_top(**global_paths["step3B_complex_topology"], properties=global_prop["step3B_complex_topology"])
        
        # STEP 3C: Minimize all H atoms in vacuum
        global_log.info("step3C_vacuum_min_Hs: Minimize all H atoms in vacuum")
        sander_mdrun(**global_paths["step3C_vacuum_min_Hs"], properties=global_prop["step3C_vacuum_min_Hs"])
        
        # STEP 3D: Process minimization results
        global_log.info("step3D_process_min: Process minimization results")
        process_minout(**global_paths["step3D_process_min"], properties=global_prop["step3D_process_min"])
        
        # STEP 3E: Minimize protein in vacuum
        global_log.info("step3E_vacuum_min_receptor: Minimize receptor in vacuum")
        all_ligands = list(ligands_dict.keys())
        global_prop["step3E_vacuum_min_receptor"]["mdin"]["restraintmask"] = f':{",".join(all_ligands)}'
        sander_mdrun(**global_paths["step3E_vacuum_min_receptor"], properties=global_prop["step3E_vacuum_min_receptor"])
        
        # STEP 3F: Process minimization results
        global_log.info("step3F_process_min: Process minimization results")
        process_minout(**global_paths["step3F_process_min"], properties=global_prop["step3F_process_min"])
        
        # NOTE: Check the ligands contain the same H atoms as in the ligand parameterization workflow or template - so far so good
        
        # NOTE: add solvent model as arg
        
        # NOTE: Add H mass repartition - add option in command line

    else:
        
        global_log.info("Using prepared structure for MD")
        global_log.info(f"Input coordinates: {incrd_path}")
        global_log.info(f"Input topology: {prmtop_path}")
        
        # If prepared structure is provided, update the global paths
        
    
    if equil_only:
        global_log.info("Equilibration only: equil_only flag is set to True! Exiting...")
        return
    
    # Debug
    return
    
    ##########################
    # Production simulations #
    ##########################
    
    if num_replicas:
        # Folder names for replicas
        simulation_folders = [f"replica_{i}" for i in range(int(num_replicas))]
        global_log.info(f"Number of replicas: {num_replicas}")
    elif num_parts:
        # Folder names for parts
        simulation_folders = [f"parts_{i}" for i in range(int(num_parts))]
        global_log.info(f"Number of parts: {num_parts}")

    # Run each simulation (replica or part)
    traj_list = []
    for simulation in simulation_folders:
        
        traj_prop = conf.get_prop_dic(prefix=simulation)
        traj_paths = conf.get_paths_dic(prefix=simulation)

        # Set general properties for all steps
        # set_global_amber_properties(traj_prop, amber_properties, global_log)

        # NOTE: Update previous global paths needed by simulation-specific steps
        
        # Enforce nsteps if provided
        if nsteps is not None:
            traj_prop['step5A_grompp_md']['mdp']['nsteps']=int(nsteps)
        total_simulation_timesteps = traj_prop['step5A_grompp_md']['mdp']['nsteps']
            
        # Simulations are replicas
        if num_replicas:
            
            # Change seed and velocities for each replica
            traj_prop['step5A_grompp_md']['mdp']['ld-seed'] = random.randint(1, 1000000)
            traj_prop['step5A_grompp_md']['mdp']['continuation'] = 'no'
            traj_prop['step5A_grompp_md']['mdp']['gen-vel'] = 'yes'

        # Simulations are parts of a single trajectory
        if num_parts:
            
            # Divide the number of steps by the number of parts
            traj_prop['step5A_grompp_md']['mdp']['nsteps']=int(traj_prop['step5A_grompp_md']['mdp']['nsteps']/int(num_parts))
            
            # For all parts except the first one, use the previous gro and cpt files
            if simulation != simulation_folders[0]:
                traj_paths['step5A_grompp_md']['incrd_path'] = previous_gro_path
                traj_paths['step5A_grompp_md']['input_cpt_path'] = previous_cpt_path

        # NOTE: In the thermostat, couple the ligands to the receptor if several groups are used

        # NOTE: Here the production run 
        global_log.info(f"{simulation} >  step5B_mdrun_md: Execute free molecular dynamics simulation")
        
        # NOTE: Append the trajectory to the list (to merge them later if needed)
        # traj_list.append(/here/the/traj/path)
        
        # NOTE: Update the previous restart and coordinates files
        
        # NOTE: Here we can include any automated analysis that we want to do with the trajectory 
        # they can be merged later or not depending on the nature of the simulations (replicas or parts)
        
        
    # Do the final analysis with all the previous parts or replicas
    if final_analysis:
        
        # If simulations are different parts of a single trajectory
        if num_parts:
            
            # Concatenate the analysis files that can be concatenated
            concatenate_gmx_analysis(conf, simulation_folders, output_path)
        
            # STEP 23: concatenate trajectories
            global_log.info("step7A_trjcat: Concatenate trajectories")
            fu.zip_list(zip_file=global_paths["step7A_trjcat"]['input_trj_zip_path'], file_list=traj_list)
            trjcat(**global_paths["step7A_trjcat"], properties=global_prop["step7A_trjcat"])
            
            # Update properties to dry the full merged trajectory
            global_prop["step7B_dry_trj"]["start"] = 0                                                                                           # The initial time of the merged trajectory in ps
            global_prop["step7B_dry_trj"]["end"] = total_simulation_timesteps*traj_prop['step5A_grompp_md']['mdp']['dt']                         # The total time of the merged trajectory in ps
            global_prop["step7B_dry_trj"]["dt"] =  traj_prop['step5A_grompp_md']['mdp']['dt']*traj_prop['step5A_grompp_md']['mdp']['nstxout']    # The saving frequency of the trajectory # NOTE: here we are hardcoding again the kind of trajectory
            
            # STEP 24: obtain dry the merged trajectory
            global_log.info("step7B_dry_trj: Obtain dry trajectory")
            gmx_trjconv_trj(**global_paths["step7B_dry_trj"], properties=global_prop["step7B_dry_trj"])
        
            #Remove unused trajectory
            os.remove(global_paths["step7A_trjcat"]["output_trj_path"])
            
            # STEP 25: obtain dry structure
            global_log.info("step7C_dry_str: Obtain dry structure")
            gmx_trjconv_str(**global_paths["step7C_dry_str"], properties=global_prop["step7C_dry_str"])

            # STEP 26: image the trajectory
            global_log.info("step7D_image_traj: Imaging the trajectory")
            gmx_image(**global_paths['step7D_image_traj'], properties=global_prop['step7D_image_traj'])
            
            #Remove unused trajectory
            os.remove(global_paths["step7B_dry_trj"]["output_traj_path"])

            # STEP 27: fit the trajectory
            global_log.info("step7E_fit_traj: Fit the trajectory")
            gmx_image(**global_paths['step7E_fit_traj'], properties=global_prop['step7E_fit_traj'])
            
            #Remove unused trajectory
            os.remove(global_paths["step7D_image_traj"]["output_traj_path"])
            
        # If simulations are replicas
        if num_replicas:
            
            # For each replica, do the final analysis
            for simulation in simulation_folders:
                
                # Get the properties and paths for the replica
                traj_prop = conf.get_prop_dic(prefix=simulation)
                traj_paths = conf.get_paths_dic(prefix=simulation)
                
                # Set general properties for all steps
                set_global_amber_properties(traj_prop, amber_properties, global_log)
                
                # NOTE: we are hard-coding the kind of traj that we are using with these paths: output_trr_path
                # Update previous global paths needed by simulation-specific steps
                traj_paths['step7B_dry_trj']['input_traj_path'] = traj_paths['step5B_mdrun_md']['output_trr_path']
                traj_paths['step7B_dry_trj']['prmtop_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7B_dry_trj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7C_dry_str']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7C_dry_str']['prmtop_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7C_dry_str']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7D_image_traj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7E_fit_traj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                
                # STEP 24: obtain dry the trajectory
                global_log.info(f"{simulation} > step7B_dry_trj: Obtain dry trajectory")
                gmx_trjconv_trj(**traj_paths["step7B_dry_trj"], properties=traj_prop["step7B_dry_trj"])
                
                # STEP 25: obtain dry structure
                global_log.info(f"{simulation} > step7C_dry_str: Obtain dry structure")
                gmx_trjconv_str(**traj_paths["step7C_dry_str"], properties=traj_prop["step7C_dry_str"])
                
                # STEP 26: image the trajectory
                global_log.info(f"{simulation} > step7D_image_traj: Imaging the trajectory")
                gmx_image(**traj_paths['step7D_image_traj'], properties=traj_prop['step7D_image_traj'])
                
                #Remove unused trajectory
                os.remove(traj_paths["step7B_dry_trj"]["output_traj_path"])
                
                # STEP 27: fit the trajectory
                global_log.info(f"{simulation} > step7E_fit_traj: Fit the trajectory")
                gmx_image(**traj_paths['step7E_fit_traj'], properties=traj_prop['step7E_fit_traj'])
                
                #Remove unused trajectory
                os.remove(traj_paths["step7D_image_traj"]["output_traj_path"])

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


if __name__ == "__main__":

    parser = argparse.ArgumentParser("MD Simulation with AMBER")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('--input_pdb', dest='input_pdb_path',
                        help="Input PDB file. Default: input_structure_path in step 1 of configuration file.",
                        required=False)

    parser.add_argument('--pdb_code', dest='pdb_code',
                        help="PDB code to get the canonical FASTA sequence of the input PDB file. If not given the workflow will look for it in the HEADER of the PDB. Default: None",
                        required=False)

    parser.add_argument('--pdb_chains', nargs='+', dest='pdb_chains',
                        help="Protein PDB chains to be extracted from PDB file and fixed. Default: A.",
                        required=False, default=['A'])

    parser.add_argument('--mutation_list', nargs='+', dest='mutation_list',
                        help="List of mutations to be introduced in the protein (default: None, ex: A:Arg220Ala)",
                        required=False)

    parser.add_argument('--ligands_topology', dest='ligands_top_folder',
                        help="Folder with .frcmod and .prep / .lib files for the ligands that should be included in the simulation. Note that the ligand coordinates should be present in the PDB. Default: None",
                        required=False)

    parser.add_argument('--skip_fix_backbone', action='store_true', dest='skip_fix_backbone',
                        help="Skip the backbone modeling of missing atoms. Default: False",
                        required=False)
    
    parser.add_argument('--skip_fix_sc', action='store_true', dest='skip_fix_side_chain',
                        help="Skip the side chain modeling of missing atoms. Default: False",
                        required=False)
    
    parser.add_argument('--fix_ss', action='store_true',
                        help="Add disulfide bonds to the protein. Use carefully! Default: False",
                        required=False)

    parser.add_argument('--fix_amides', action='store_true', dest='fix_amide_clashes',
                        help="Flip clashing amides to relieve the clashes. Default: False",
                        required=False)
    
    parser.add_argument('--ff', dest='forcefields',
                        nargs='+', help="Forcefields to be used in the simulation. Default: ['protein.ff14SB','DNA.bsc1','gaff']",
                        required=False, default=['protein.ff14SB','DNA.bsc1','gaff'])

    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system. Default: False",
                        required=False, default=False)

    parser.add_argument('--incrd', dest='incrd_path',
                        help="Input structure coordinates file ready to minimize (.incrd). To provide an externally prepared system, use together with --prmtop (default: None)",
                        required=False)

    parser.add_argument('--prmtop', dest='prmtop_path',
                        help="Input topology file ready to minimize (.prmtop). To provide an externally prepared system, use together with --incrd (default: None)",
                        required=False)
    
    parser.add_argument('--equil_only', action='store_true',
                        help="Only run the equilibration steps. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--nsteps', dest='nsteps',
                        help="Number of steps of the production simulation",
                        required=False)
    
    parser.add_argument('--num_parts', dest='num_parts',
                        help="Number of parts to divide the simulation into. Default: 1",
                        required=False)
    
    parser.add_argument('--num_replicas', dest='num_replicas',
                        help="Number of replicas with different seeds to run the simulation. Default: 1",
                        required=False)

    parser.add_argument('--final_analysis', action='store_true', dest='final_analysis',
                        help="Run the final analysis of the trajectory/ies. Concatenation of the analysis and trajectory, trajectory drying, imaging and fitting. Default: False",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    # NOTE: add flag to determine what should remain restrained during the production run - currently everything is free

    args = parser.parse_args()

    # Check .gro structure and .zip topology are given together
    if (args.incrd_path is None and args.prmtop_path is not None) or (args.incrd_path is not None and args.prmtop_path is None):
        raise Exception("Both --input_gro and --prmtop_path must be provided together")

    # Check .pdb structure and .gro/.zip topology are not given together
    if (args.input_pdb_path is not None and args.incrd_path is not None):
        raise Exception("Both --input_pdb and --input_gro/--prmtop_path are provided. Please provide only one of them")

    main_wf(configuration_path=args.config_path, input_pdb_path=args.input_pdb_path, pdb_code=args.pdb_code, 
            pdb_chains=args.pdb_chains, mutation_list=args.mutation_list, ligands_top_folder=args.ligands_top_folder, 
            skip_fix_backbone=args.skip_fix_backbone, skip_fix_side_chain=args.skip_fix_side_chain, fix_ss=args.fix_ss, 
            fix_amide_clashes=args.fix_amide_clashes, forcefields=args.forcefields, setup_only=args.setup_only, 
            incrd_path=args.incrd_path, prmtop_path=args.prmtop_path, equil_only=args.equil_only, nsteps=args.nsteps,  
            num_parts=args.num_parts, num_replicas=args.num_replicas, final_analysis=args.final_analysis, 
            output_path=args.output_path)