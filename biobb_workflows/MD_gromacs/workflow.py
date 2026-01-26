#!/usr/bin/env python3

# Importing all the needed libraries
from typing import List, Dict, Union, Optional, Literal
from pathlib import Path
import logging
import argparse
import shutil
import time
import os

from Bio.PDB import PDBParser

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.convert_tpr import convert_tpr
from biobb_gromacs.gromacs.genion import genion
from biobb_gromacs.gromacs.mdrun import mdrun
from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_gromacs.gromacs.genrestr import genrestr
from biobb_gromacs.gromacs_extra.append_ligand import append_ligand
from biobb_gromacs.gromacs_extra.ndx2resttop import ndx2resttop
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf 
from biobb_structure_utils.utils.cat_pdb import cat_pdb

# Constants
# Titratable residues in GROMACS, HIS are already set with their resnames
# 'HIS': ('HISD', 'HISE', 'HISH', 'HIS1')
gmx_titra_resnames = {
    'LYS': ('LYSN', 'LYS'),                 # deprotonated, protonated
    'ARG': ('ARGN', 'ARG'),                 # deprotonated, protonated
    'ASP': ('ASP', 'ASPH'),                 # deprotonated, protonated
    'GLU': ('GLU', 'GLUH')                  # deprotonated, protonated
}

def check_inputs(
    input_pdb_path: Optional[str], 
    input_gro_path: Optional[str], 
    input_top_path: Optional[str], 
    ligands_top_folder: Optional[str], 
    input_tpr_path: Optional[str],
    input_cpt_path: Optional[str],
    input_ndx_path: Optional[str],
    input_plumed_path: Optional[str],
    input_plumed_folder: Optional[str],
    setup_only: Optional[bool],
    equil_only: Optional[bool],
    global_log: Optional[logging.Logger]
) -> Literal['input_pdb', 'prepared_system', 'restart_simulation']:
    """
    Check the inputs for the workflow. 
    
    1. Check if any of the compulsory input files are provided and if they exist.
    
    2. Check if the provided optional input files exist.
    
    If any of the checks fail, an error is raised.
    
    Inputs
    ------
    
        input_pdb_path (str): Path to the input PDB file.
        input_gro_path (str): Path to the input GRO file.
        input_top_path (str): Path to the input topology file.
        ligands_top_folder (str): Path to the folder with the ligands .itp files.
        input_tpr_path (str): Path to already-prepared binary input run file
        input_cpt_path (str): Path to checkpoint file.
        input_ndx_path (str): Path to the index file.
        input_plumed_path (str): Path to the plumed file.
        input_plumed_folder (str): Path to the folder with the plumed files.
        setup_only (bool): Condition for the user to request just the setup of the simulation
        equil_only (bool): Condition for the user to request just the equilibration of the simulation
        global_log: Logger object for logging messages.
    """
    
    # Types of input modes
    input_modes = {
        'input_pdb': {
            'compulsory' : [input_pdb_path],
            'optional' : [ligands_top_folder]
        },
        'prepared_system': {
            'compulsory' : [input_gro_path, input_top_path],
            'optional' : []
        },
        'restart_simulation': {
            'compulsory' : [input_tpr_path, input_cpt_path],
            'optional' : []
        }
    }
    
    # Check how many input modes are being used
    used_modes = []
    for mode, files in input_modes.items():
        # Check if the compulsory files are given
        if all(files['compulsory']):
            used_modes.append(mode)

    if len(used_modes) == 0:
        global_log.error("No valid inputs found. Provide either:\n"
                         "1) An input PDB file (input_pdb_path) to prepare the system from scratch.\n"
                         "2) An input GRO file (input_gro_path) and a topology file (input_top_path) to run MD on a prepared system.\n"
                         "3) An input TPR file (input_tpr_path) and a checkpoint file (input_cpt_path) to restart a simulation.")
        raise FileNotFoundError("No valid inputs found.")

    if len(used_modes) > 1:
        global_log.error(f"Incompatible inputs found: {used_modes}. Please provide inputs for only one mode.")
        raise ValueError("Incompatible inputs found.")

    # Check if the compulsory files exist
    for file in input_modes[used_modes[0]]['compulsory']:
        if not os.path.exists(file):
            global_log.error(f"File {file} not found.")
            raise FileNotFoundError(f"File {file} not found.")
    
    # Check the optional files exist if there are any
    for file in input_modes[used_modes[0]]['optional']:
        if file:
            if not os.path.exists(file):
                global_log.error(f"File {file} not found.")
                raise FileNotFoundError(f"File {file} not found.")
    
    # Set up only option requires input pdb mode
    if setup_only:
        if not used_modes[0] == 'input_pdb':
            global_log.error("The option setup_only = True requires a PDB file as input.")
    
    # Equil only option requires input_pdb or prepared_system
    if equil_only:
        if used_modes[0] == 'restart_simulation':
             global_log.error("The option equil_only = True requires a PDB file or .gro and .top files as input.")
    
    # Log input mode
    if used_modes[0] == 'input_pdb':
        global_log.info("Using PDB as input:")
        global_log.info(f"Input PDB file: {input_pdb_path}")
        
    elif used_modes[0] == 'prepared_system':
        global_log.info("Using prepared system as input:")
        global_log.info(f"Input GRO file: {input_gro_path}")
        global_log.info(f"Input TOP file: {input_top_path}")
        
    elif used_modes[0] == 'restart_simulation':
        global_log.info("Using restart files as input:")
        global_log.info(f"Input TPR file: {input_tpr_path}")
        global_log.info(f"Input CPT file: {input_cpt_path}")
        
    # If index file is provided, check it exists
    if input_ndx_path:
        global_log.info(f"Input NDX file: {input_ndx_path}")
        if not os.path.exists(input_ndx_path):
            global_log.error(f"File {input_ndx_path} not found.")
            raise FileNotFoundError(f"File {input_ndx_path} not found.")
    
    # If plumed file is provided, check it exists
    if input_plumed_path:
        global_log.info(f"Input PLUMED file: {input_plumed_path}")
        if not os.path.exists(input_plumed_path):
            global_log.error(f"File {input_plumed_path} not found.")
            raise FileNotFoundError(f"File {input_plumed_path} not found.")

    # If plumed folder is provided, check it exists
    if input_plumed_folder:
        global_log.info(f"Input PLUMED folder: {input_plumed_folder}")
        if not os.path.exists(input_plumed_folder):
            global_log.error(f"Folder {input_plumed_folder} not found.")
            raise FileNotFoundError(f"Folder {input_plumed_folder} not found.")
        # Check the folder is not empty
        if not os.listdir(input_plumed_folder):
            global_log.error(f"Folder {input_plumed_folder} is empty.")
            raise FileNotFoundError(f"Folder {input_plumed_folder} is empty.")

    return used_modes[0]

# Biopython helpers
def get_ligands(ligands_top_folder: Union[str, None], global_log) -> List[Dict[str, str]]:
    """
    Get a list of available ligands in the ligands topology folder. The function searches for all the 
    .itp and .gro files in the folder. If the folder is provided but doesn't exist or any .itp/.gro file 
    is missing, an error is raised.
    
    Inputs
    ------
    
        ligands_top_folder (str): Path to the folder with the ligands .itp files.
        global_log: Logger object for logging messages.
    
    Returns
    -------
    
        ligands: Dictionary with the ligand names, topology and coordinate file paths. Empty dict if no ligands are found. 
        
            Example 
                    ligands = {
                        'ZZ7': {
                            'topology': 'path/to/ZZ7.itp',
                            'coordinates': 'path/to/ZZ7.gro'
                        },
                        'ZZ8': {
                            'topology': 'path/to/ZZ8.itp',
                            'coordinates': 'path/to/ZZ8.gro'
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
        if file.endswith(".itp") or file.endswith(".gro"):
            
            # Get the file name without extension
            ligand_id = Path(file).stem   
            
            # Get the file extension
            file_extension = Path(file).suffix
            
            # Check if the ligand name is already in the dictionary
            if ligand_id in ligands:
                if file_extension == ".itp":
                    ligands[ligand_id]['topology'] = os.path.join(ligands_top_folder, file)
                elif file_extension == ".gro":
                    ligands[ligand_id]['coordinates'] = os.path.join(ligands_top_folder, file)
            else:
                if file_extension == ".itp":
                    ligands[ligand_id] = {'topology': os.path.join(ligands_top_folder, file)}
                elif file_extension == ".gro":
                    ligands[ligand_id] = {'coordinates': os.path.join(ligands_top_folder, file)}
    
    # Check if all ligands have both topology and coordinate files
    for ligand, files in ligands.items():
        if 'topology' not in files:
            global_log.error(f"Topology file for ligand {ligand} not found")
            
            # Remove ligand from the dictionary
            ligands.pop(ligand)
            continue
        
        if 'coordinates' not in files:
            global_log.error(f"Coordinate file for ligand {ligand} not found")
            
            # Remove ligand from the dictionary
            ligands.pop(ligand)
            continue
        
    return ligands
    
def get_chains_dict(pdb_file: str) -> Dict[str, Dict[str, List[int]]]:
    """ 
    Get a dictionary with the chain IDs as keys and a dictionary with the residue intervals as values.
    
    HETEROATOM residues are not considered. 
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.
    
    Returns
    -------
    Dict:
        A dictionary with the chain IDs as keys and a list with the start and end residue numbers.
        
        Example:
        {
            'A': {'residues': [1, 100]},
            'B': {'residues': [1, 200]}
        }
    """
    
    # NOTE: What happens with nucleotides
    # NOTE: what happens if there are no chains
    
    # Initialize the PDB parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein_structure', pdb_file)
    
    chains_dict = {}
    
    # Iterate over all models and chains in the structure
    for model in structure:
        for chain in model:
            # List to store residue numbers for standard amino acids/nucleotides
            standard_residues = []
            
            for residue in chain:
                # The residue ID: (hetfield, sequence_identifier, insertion_code)
                # For standard residues, hetfield is a blank space ' '
                if residue.get_id()[0] == ' ':
                    standard_residues.append(residue.get_id()[1])
            
            # Only add the chain to the dictionary if it contains standard residues
            if standard_residues:
                chain_id = chain.get_id()
                chains_dict[chain_id] = {
                    'residues': [min(standard_residues), max(standard_residues)]
                }
                
    return chains_dict
   
def read_protonation_states(pdb_file: str, resname: str, global_log) -> List[str]:
    """ 
    Read a PDB file and return the protonation states per chain of the specified residue 
    name along the protein. The protonation states are given using the pdb2gmx naming convention:
    
        0: deprotonated
        1: protonated
    
    Parameters
    ----------
    
        pdb_file : str
            Path to the PDB file.
        resname : str
            Residue name to check the protonation state. Example: "HIS"
        global_log : Logger
            Logger object for logging messages.

    Returns
    -------
    
        List[str]: 
            Protonation state of the specified residue name.
            
            Example: ["0 1 0 1", "", "0 1 0 1"] 
    """
    
    # Find the list of protonation state resnames for this residue
    protonation_resnames = gmx_titra_resnames.get(resname)
    if not protonation_resnames:
        global_log.error(f"Residue {resname} not found in the list of titratable residues.")
        return ""
    
    # Initialize
    line_count = 0
    old_pdb_resnum = 0
    old_chain_id = ""
    protonation_states = ""
    protonation_states_list = []
    
    # Parse the PDB structure manually
    with open(pdb_file, 'r') as f:
        for line in f:
            line_count += 1
            # If line contains a residue atom
            if len(line) > 26 and line.startswith("ATOM"):
                # Read residue name columns (18-20) and column 21 (to read 4-letter resnames)
                pdb_resname = line[17:21].strip()
                # Read chain ID column (22)
                pdb_chain_id = line[21:22].strip()
                # Read residue sequence number columns (23-26)
                pdb_resnum = line[22:26].strip()
                
                # Update old chain ID on first line
                if old_chain_id == "":
                    old_chain_id = pdb_chain_id
                
                # If this line is a new chain, add protonation states to the list
                if pdb_chain_id != old_chain_id:
                    # Add the protonation states to the list
                    protonation_states_list.append(protonation_states.strip())
                    # Reset the protonation states
                    protonation_states = ""
                    # Update the old chain ID
                    old_chain_id = pdb_chain_id
                
                try:
                    pdb_resnum = int(pdb_resnum)
                except ValueError:
                    raise ValueError(f"""{pdb_resnum} cannot be converted to int. Could not convert 
                                     residue sequence columns 23-26 to integer on line {line_count} 
                                     of {pdb_file}. Check the PDB file format.""")
                # If this line is a new match for the type of residue. Example: "HISD", "HISE"...
                if resname in pdb_resname and pdb_resnum > old_pdb_resnum:
                    # Find the protonation state of this residue
                    protonation_states += f"{protonation_resnames.index(pdb_resname)} "
                    # Update the old residue number
                    old_pdb_resnum = pdb_resnum
        
        # Add the last protonation states to the list
        protonation_states_list.append(protonation_states.strip())
    
    return protonation_states_list

# Process topology - temporal solution 
def process_ligand_top(input_path: str, output_path: str) -> None:
    """
    Read the input topology from the ligand. 
    Removes any [ defaults ] directive present.
    Removes any [ molecules ] directive present. 
    Copies the input topology to the output path.
    
    Inputs
    ------
    
        input_path (str): Path to the input topology file.
        output_path (str): Path to the output topology file.
    
    Returns
    -------
    
        None
    """
    
    # Read the input topology file
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    
    # Remove any defaults and molecules directives
    new_lines = []
    reading_defaults = False
    reading_molecules = False
    for line in lines:
        
        # Mark the end of the defaults section
        if reading_defaults and line.startswith("["):
            reading_defaults = False
            
        # Mark the end of the molecules section
        if reading_molecules and line.startswith("["):
            reading_molecules = False
            
        # Mark the beginning of the defaults section
        if line.startswith("[ defaults ]"):
            reading_defaults = True 
        
        # Mark the beginning of the molecules section
        if line.startswith("[ molecules ]"):
            reading_molecules = True

        # Add the line if not reading the defaults section or the molecules section
        if not reading_defaults and not reading_molecules:
            new_lines.append(line)
    
    # Write the new topology file
    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    

# YML construction
def config_contents(
    gmx_bin: str = 'gmx',
    mpi_bin: Optional[str] = None,
    mpi_np: Optional[int] = None,
    num_threads_mpi: Optional[int] = 0,
    num_threads_omp: Optional[int] = 0,
    use_gpu: bool = False,
    restart: bool = False,
    forcefield: str = 'amber99sb-ildn',
    ions_concentration: float = 0.15,
    temp: float = 300.0,
    dt: float = 2.0,
    equil_time: float = 1.0,
    equil_frames: Optional[int] = 500,
    prod_time: Optional[float] = 100.0,
    prod_frames: Optional[int] = 2000,
    seed: Optional[int] = -1,
    debug: bool = False,
    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.
    
    The YAML file contains the configuration for the protein preparation workflow.
    
    Parameters
    ----------
    
    gmx_bin : str
        Path to the GROMACS binary. Default is 'gmx'.
    mpi_bin : str
        Path to the MPI binary. Default is 'null'.
    mpi_np : int
        Number of MPI ranks to start by mpi binary (e.g. mpirun -np N). Only used when mpi_bin is provided. Use only when running multi-node with gmx_mpi.
        When using srun, this parameter is not needed as srun automatically allocates the number of ranks specified by the job scheduler.
    num_threads_mpi : int
        Number of ranks for MPI. Default is None.
    num_threads_omp : int
        Number of threads for OpenMP. Default is None.
    use_gpu : bool
        Whether to use GPU for GROMACS. Default is False.
    restart : bool
        Whether to skip steps already performed. Default is True.
    forcefield : str
        Force field to use. Default is 'amber99sb-ildn'.
    ions_concentration : float
        Concentration of salt in mols/L. Default is 0.15.
    temp : float
        Temperature in Kelvin. Default is 300.0.
    dt : float
        Time step in fs. Default is 2 fs.
    equil_time : float
       Time of each equilibration step in ns. Default: 1 ns
    equil_frames : int
        Number of frames to save during the equilibration steps. Default: 500 frames.
    prod_time : float
        Total time of the production simulation in ns. Default: 100 ns
    prod_frames : int
        Number of frames to save during the production steps. Default: 2000 frames.
    seed: int
        Seed for random number generation. Default is -1 (random seed).

    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    # Check dt is between 1 and 4 fs
    if dt < 1.0 or dt > 4.0:
        raise ValueError("dt must be between 1 and 4 fs.")
    
    # Check equil_time and prod_time are larger than 0
    if equil_time <= 0.0:
        raise ValueError("equil_time must be larger than 0 ns.")
    
    if prod_time <= 0.0:
        raise ValueError("prod_time must be larger than 0 ns.")
    
    # Calculate number of steps
    nsteps_equil = int(equil_time*10**6 // dt)
    nsteps_prod = int(prod_time*10**6 // dt)
    
    dt_in_ps = dt / 1000  # convert fs to ps

    equil_traj_freq_steps = max(nsteps_equil // equil_frames, 1)  # equil_frames frames during equilibration

    prod_traj_freq_steps = max(nsteps_prod // prod_frames, 1)      # prod_frames frames during production
    
    if mpi_bin is None:
        mpi_bin = 'null'
    
    if mpi_np is None:
        mpi_np_config = ''
    else:
        mpi_np_config = f"mpi_np: {mpi_np}"

    num_threads_mpi_config = f"num_threads_mpi: {num_threads_mpi}"
    num_threads_omp_config = f"num_threads_omp: {num_threads_omp}"
    
    if seed is None:
        seed = -1

    return f""" 
# Global properties (common for all steps)
global_properties:
  working_dir_path: output                                          # Workflow default output directory
  can_write_console_log: False                                      # Verbose writing of log information
  restart: {restart}                                                # Skip steps already performed
  remove_tmp: {not debug}                                               # Remove temporal files

##################################################################
# Section 3 (Steps A-I): Prepare topology and coordinates for MD #
##################################################################

# Generate the topology of the structure with pdb2gmx
step1_pdb2gmx:
  tool: pdb2gmx
  paths:
    input_pdb_path: path/to/input.pdb                     # Will be set by the workflow
    output_gro_path: structure.gro
    output_top_zip_path: structure_top.zip
  properties:
    binary_path: {gmx_bin}    
    force_field: {forcefield}     # Will be set by the workflow
    water_type: tip3p             # spc, spce, tip3p, tip4p, tip5p, tips3p
    ignh: False                   # Ignore hydrogens in the input structure
    merge: False                  # Merge all chains into one molecule

# Add the reference group to the index file
step1_make_ref_group:
  tool: make_ndx
  paths:
    input_structure_path: dependency/step1_pdb2gmx/output_gro_path
    output_ndx_path: chain.ndx
  properties:
    binary_path: {gmx_bin}   
    selection:  "System"         # Will be set by the workflow

# Add the restrained group to the index file
step2_make_rest_group:
  tool: make_ndx
  paths:
    input_structure_path: dependency/step1_pdb2gmx/output_gro_path
    input_ndx_path: dependency/step1_make_ref_group/output_ndx_path
    output_ndx_path: calpha.ndx
  properties:
    binary_path: {gmx_bin}   
    selection: "a CA"            # Will be set by the workflow

# Append position restraints to the topology file using the reference and restrained groups of the index file
step3_append_posres:
  tool: ndx2resttop
  paths:
    input_ndx_path: dependency/step2_make_rest_group/output_ndx_path
    input_top_zip_path: dependency/step1_pdb2gmx/output_top_zip_path
    output_top_zip_path: structure_top.zip
  properties:
    force_constants: "500 500 500"

step4_structure_pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: dependency/step1_pdb2gmx/output_gro_path
    input_structure_path: dependency/step1_pdb2gmx/output_gro_path
    output_str_path: structure.pdb
  properties:
    binary_path: {gmx_bin}   

step1_create_ligand_pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: path/to/ligand.gro        # Will be set by the workflow
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    output_str_path: ligand.pdb
  properties:   
    binary_path: {gmx_bin}   

step2_create_complex_pdb:
  tool: cat_pdb
  paths: 
    input_structure1: dependency/step4_structure_pdb/output_str_path  # Will be set by the workflow
    input_structure2: dependency/step1_create_ligand_pdb/output_str_path
    output_structure_path: complex.pdb

step3_make_ligand_ndx:
  tool: make_ndx
  paths:
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    output_ndx_path: ligand_heavy_atoms.ndx
  properties:
    binary_path: {gmx_bin}  
    selection: "0 & ! a H*"

step4_create_ligand_restraints:
  tool: genrestr
  paths:
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    input_ndx_path: dependency/step3_make_ligand_ndx/output_ndx_path
    output_itp_path: ligand_restraints.itp    # Will be set by the workflow
  properties:
    binary_path: {gmx_bin} 
    restrained_group: "System_&_!H*"
    force_constants: "500 500 500"

step5_append_ligand_topology:
  tool: append_ligand
  paths:
    input_top_zip_path: dependency/step3_append_posres/output_top_zip_path  # Will be set by the workflow
    input_itp_path: path/to/ligand.itp                                            # Will be set by the workflow
    input_posres_itp_path: dependency/step4_create_ligand_restraints/output_itp_path    # Path to ligand position restraint topology file
    output_top_zip_path: complex_top.zip

step6_editconf:
  tool: editconf
  paths:
    input_gro_path: dependency/step1_pdb2gmx/output_gro_path
    output_gro_path: editconf.gro
  properties:
    binary_path: {gmx_bin}   
    box_type: octahedron         # cubic, triclinic, octahedron, dodecahedron
    distance_to_molecule: 1.0    # Distance of the box from the outermost atom in nm

step7_solvate:
  tool: solvate
  paths:
    input_solute_gro_path: dependency/step6_editconf/output_gro_path
    input_top_zip_path: dependency/step3_append_posres/output_top_zip_path
    output_top_zip_path: solvate_top.zip
    output_gro_path: solvate.gro
  properties:
    binary_path: {gmx_bin}     
    
step8_grompp_genion:
  tool: grompp
  paths:
    input_gro_path: dependency/step7_solvate/output_gro_path
    input_top_zip_path: dependency/step7_solvate/output_top_zip_path
    output_tpr_path: gppion.tpr
  properties:
    binary_path: {gmx_bin}     
    simulation_type: minimization
    maxwarn: 10                    # NOTE: this will be ligand dependent!! :o - the warning is for each atom with redefined parameters

step9_genion:
  tool: genion
  paths:
    input_tpr_path: dependency/step8_grompp_genion/output_tpr_path
    input_top_zip_path: dependency/step7_solvate/output_top_zip_path
    output_top_zip_path: genion_top.zip
    output_gro_path: genion.gro
  properties:
    binary_path: {gmx_bin}    
    neutral: True                 # Neutralize charge of the system
    concentration: {ions_concentration}    # Concentration of ions in mols/L

#############################################################################
# Section 4 (Steps A-I): Minimize and equilibrate the initial configuration #
#############################################################################

# NOTE: Leverage hierarchy in mdp when creating new eq steps without restraints - you can overwrite the default mdp config from the type of simulation

step1_grompp_min:
  tool: grompp
  paths:
    input_gro_path: dependency/step9_genion/output_gro_path
    input_top_zip_path: dependency/step9_genion/output_top_zip_path
    output_tpr_path: gppmin.tpr
  properties:
    binary_path: {gmx_bin}     
    simulation_type: minimization
    mdp:
      integrator: steep
      nsteps: 10000
      emtol: 500
      emstep: 0.01
    
step2_mdrun_min:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step1_grompp_min/output_tpr_path
    output_trr_path: min.trr     # NOTE: Make sure we are not writing f***g trr files
    output_gro_path: min.gro
    output_edr_path: min.edr
    output_log_path: min.log
  properties:
    binary_path: {gmx_bin} 
    mpi_bin: {mpi_bin}
    {num_threads_mpi_config}
    {num_threads_omp_config}
    {mpi_np_config}
    
step3_make_ndx:
  tool: make_ndx 
  paths:
    input_structure_path: dependency/step2_mdrun_min/output_gro_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '! "Water_and_ions"'

step4_energy_min:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step2_mdrun_min/output_edr_path
    output_xvg_path: min_ene.xvg
  properties:
    binary_path: {gmx_bin}     
    terms: ["Potential"]
    xvg: xmgr 

step5_grompp_nvt:
  tool: grompp
  paths:
    input_gro_path: dependency/step2_mdrun_min/output_gro_path
    input_ndx_path: dependency/step3_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step9_genion/output_top_zip_path
    output_tpr_path: gppnvt.tpr
  properties:
    binary_path: {gmx_bin}     
    simulation_type: nvt
    mdp:
      ref-t: {temp} {temp}
      tc-grps: "!Water_and_ions Water_and_ions"
      nsteps: {nsteps_equil}
      dt: {dt_in_ps}              
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {equil_traj_freq_steps}
      
step6_mdrun_nvt:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step5_grompp_nvt/output_tpr_path
    output_xtc_path: nvt.xtc
    output_trr_path: nvt.trr
    output_gro_path: nvt.gro
    output_edr_path: nvt.edr
    output_log_path: nvt.log
    output_cpt_path: nvt.cpt
  properties:
    binary_path: {gmx_bin}     
    mpi_bin: {mpi_bin}               
    use_gpu: {use_gpu}             # Wether to use GPU support or not
    {num_threads_mpi_config}
    {num_threads_omp_config}
    {mpi_np_config}
    
step7_temp_nvt:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step6_mdrun_nvt/output_edr_path
    output_xvg_path: nvt_temp.xvg
  properties:
    binary_path: {gmx_bin}     
    terms: ["Temperature"]
    xvg: xmgr 

step8_grompp_npt:
  tool: grompp
  paths:
    input_gro_path: dependency/step6_mdrun_nvt/output_gro_path
    input_ndx_path: dependency/step3_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step9_genion/output_top_zip_path
    input_cpt_path: dependency/step6_mdrun_nvt/output_cpt_path
    output_tpr_path: gppnpt.tpr
  properties:
    binary_path: {gmx_bin}     
    simulation_type: npt
    mdp:
      pcoupltype: isotropic
      nsteps: {nsteps_equil}
      dt: {dt_in_ps}
      ref-t: {temp} {temp}
      tc-grps: "!Water_and_ions Water_and_ions"
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {equil_traj_freq_steps}

step9_mdrun_npt:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step8_grompp_npt/output_tpr_path
    output_xtc_path: npt.xtc
    output_trr_path: npt.trr
    output_gro_path: npt.gro
    output_edr_path: npt.edr
    output_log_path: npt.log
    output_cpt_path: npt.cpt
  properties: 
    binary_path: {gmx_bin}     
    mpi_bin: {mpi_bin}               
    use_gpu: {use_gpu}
    {num_threads_mpi_config}
    {num_threads_omp_config}
    {mpi_np_config}

step10_density_npt:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step9_mdrun_npt/output_edr_path
    output_xvg_path: npt_press_den.xvg
  properties:
    binary_path: {gmx_bin}     
    terms: ["Pressure", "Density"]
    xvg: xmgr 

############################################
# Section 5 (Steps A-B): MD production run #
############################################

step1_grompp_md:
  tool: grompp
  paths:
    input_gro_path: dependency/step9_mdrun_npt/output_gro_path
    input_cpt_path: dependency/step9_mdrun_npt/output_cpt_path 
    input_ndx_path: dependency/step3_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step9_genion/output_top_zip_path
    output_tpr_path: gppmd.tpr
  properties:
    binary_path: {gmx_bin}     
    simulation_type: free
    mdp:
      nsteps: {nsteps_prod}
      dt: {dt_in_ps} 
      ref-t: {temp} {temp}
      tc-grps: "!Water_and_ions Water_and_ions"
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {prod_traj_freq_steps}
      nstenergy: {prod_traj_freq_steps}
      continuation: 'no'
      gen-vel: 'yes'          
      ld-seed: {seed}
      gen-seed: {seed}

step1B_convert_tpr:
    tool: convert-tpr
    paths:
      input_tpr_path: dependency/step1_grompp_md/output_tpr_path
      output_tpr_path: new_run_file.tpr
    properties:
      extend: {prod_time*1000}   # Extend simulation by prod_time in ps
      binary_path: {gmx_bin}     

step2_mdrun_prod:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step1_grompp_md/output_tpr_path
    output_xtc_path: md.xtc
    output_trr_path: md.trr
    output_gro_path: md.gro
    output_edr_path: md.edr
    output_log_path: md.log
    output_cpt_path: md.cpt
  properties:
    binary_path: {gmx_bin}     
    mpi_bin: {mpi_bin}               
    use_gpu: {use_gpu}
    {num_threads_mpi_config}
    {num_threads_omp_config}
    {mpi_np_config}
        
#########################################
# Section 6 (Steps A-D): Basic analysis #
#########################################

step1_gro2pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: dependency/step2_mdrun_prod/output_gro_path
    input_structure_path: dependency/step2_mdrun_prod/output_gro_path
    output_str_path: topology.pdb
  properties:
    binary_path: {gmx_bin}   
    
step2_rmsd_equilibrated:
  tool: gmx_rms
  paths:
    input_structure_path: dependency/step1_gro2pdb/output_str_path
    input_traj_path: dependency/step2_mdrun_prod/output_xtc_path
    output_xvg_path: md_rmsdfirst.xvg
  properties:
    binary_path: {gmx_bin}     
    selection: Backbone
    xvg: xmgr 

step3_rmsd_experimental:
  tool: gmx_rms
  paths:
    input_structure_path: dependency/step1_gro2pdb/output_str_path
    input_traj_path: dependency/step2_mdrun_prod/output_xtc_path
    output_xvg_path: md_rmsdexp.xvg
  properties:
    binary_path: {gmx_bin}     
    selection: Backbone
    xvg: xmgr 

step4_rgyr:
  tool: gmx_rgyr
  paths:
    input_structure_path: dependency/step1_gro2pdb/output_str_path
    input_traj_path: dependency/step2_mdrun_prod/output_xtc_path
    output_xvg_path: md_rgyr.xvg
  properties:
    binary_path: {gmx_bin}     
    selection: Backbone
    xvg: xmgr 
    
step5_rmsf:
  tool: cpptraj_rmsf
  paths:
    input_top_path: dependency/step1_gro2pdb/output_str_path
    input_traj_path: dependency/step2_mdrun_prod/output_xtc_path
    output_cpptraj_path: md_rmsf.xmgr   # .dat, .agr, .xmgr, .gnu
  properties:
    start: 1
    end: -1
    steps: 1
    mask: "!@H=" # by default cpptraj already strips solvent atoms

#####################################################
# Section 7 (Steps A-D): Trajectory post-processing #
#####################################################

step6_dry_str:
  tool: gmx_trjconv_str
  paths: 
    input_structure_path: dependency/step1_gro2pdb/output_str_path
    input_top_path: dependency/step1_grompp_md/output_tpr_path
    output_str_path: dry_structure.gro
  properties:
    binary_path: {gmx_bin}     
    selection: Protein
    center: True
    pbc: mol
    ur: compact

step7_dry_traj:
  tool: gmx_trjconv_trj
  paths: 
    input_traj_path: dependency/step2_mdrun_prod/output_xtc_path
    input_top_path: dependency/step1_grompp_md/output_tpr_path
    output_traj_path: dry_traj.xtc
  properties:
    binary_path: {gmx_bin}     
    selection: Protein

step8_center:
  tool: gmx_image 
  paths:
    input_traj_path: dependency/step7_dry_traj/output_traj_path
    input_top_path: dependency/step1_grompp_md/output_tpr_path
    output_traj_path: center_traj.xtc
  properties:
    binary_path: {gmx_bin}     
    center_selection: Protein
    output_selection: Protein
    center: True
    ur: compact
    pbc: none

step9_image_traj:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step8_center/output_traj_path
    input_top_path: dependency/step1_grompp_md/output_tpr_path
    output_traj_path: imaged_traj.xtc
  properties:
    binary_path: {gmx_bin}     
    output_selection: Protein
    cluster_selection: Protein
    center_selection: Protein  # NOTE: why is this used??
    center: False
    ur: compact
    pbc: mol

step10_fit_traj:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step9_image_traj/output_traj_path
    input_top_path: dependency/step1_grompp_md/output_tpr_path
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_bin}     
    fit_selection: Protein
    center_selection: Protein
    output_selection: Protein
    center: False
    fit: rot+trans
"""

def create_config_file(config_path: str, **config_args) -> bool:
    """
    Create a YAML configuration file for the workflow if needed.
    Check if the config_path is None or does not exist before creating the file.
    If the file already exists, it will not be overwritten.
    Return a boolean indicating whether the file was created or not.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file to be created.
    config_args : dict
        Arguments to be used in the configuration file.
    
    Returns
    -------
    
    bool
        True if the file was created, False if it already exists.
    """
    
    # Check if the config_path is None
    if config_path is not None:
        # Check if the file already exists
        if os.path.exists(config_path):
            print(f"Configuration file already exists at {config_path}.")
            return False
        else:
            print(f"Warning: Configuration file path is set to {config_path}, but it does not exist. Creating a new config.yml file.")
        
    config_path = 'config.yml'
        
    # Write the contents to the file
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))
        
    print(f"Configuration file created at {config_path}.")
    
    # Return True indicating the file was created
    return True
    
    
# Main workflow
def main_wf(input_pdb_path: Optional[str] = None, 
            ligands_top_folder: Optional[str] = None, 
            input_gro_path: Optional[str] = None, 
            input_top_path: Optional[str] = None, 
            input_tpr_path: Optional[str] = None, 
            input_cpt_path: Optional[str] = None,
            input_ndx_path: Optional[str] = None,
            input_plumed_path: Optional[str] = None,
            input_plumed_folder: Optional[str] = None,
            configuration_path: Optional[str] = None, 
            gmx_bin: Optional[str] = 'gmx',
            mpi_bin: Optional[str] = None,
            mpi_np: Optional[int] = None,
            num_threads_mpi: Optional[int] = 0,
            num_threads_omp: Optional[int] = 0,
            use_gpu: Optional[bool] = False,
            restart: Optional[bool] = False,
            forcefield: Optional[str] = 'amber99sb-ildn', 
            ions_concentration:  Optional[float] = 0.15,
            temperature:  Optional[float] = 300.0,
            random_seed: Optional[int] = None,
            setup_only: Optional[bool] = False, 
            dt: Optional[float] = 2.0,
            equil_time: Optional[float] = 1,
            equil_frames: Optional[int] = 500,
            equil_only: Optional[bool] = False, 
            prod_time: Optional[float] = 100,
            prod_frames: Optional[int] = 2000,
            debug: Optional[bool] = False,
            output_path: Optional[str] = None
    ):
    '''
    Main MD setup and run workflow with GROMACS. Can be used to prepare and launch an MD simulation.

    Inputs
    ------

        input_pdb_path: 
            path to input PDB file
        ligands_top_folder:
            path to folder with ligands topology files (.itp) and coordinates (.gro)
        input_gro_path: 
            path to already-prepared input structure file (.gro)
        input_top_path: 
            path to already-prepared input topology file (.zip)
        input_tpr_path:
            path to already-prepared binary input run file (.tpr)
        input_cpt_path:
            path to input checkpoint file (.cpt)
        input_ndx_path:
            path to input index file (.ndx)
        input_plumed_path:
            path to the main PLUMED input file. If provided, PLUMED will be used during the simulation. (.dat)
        input_plumed_folder:
            path to folder with PLUMED input files if needed
        configuration_path: 
            path to YAML configuration file
        gmx_bin:
            path to GROMACS binary, either gmx (single-node) or gmx_mpi (multi-node)
        mpi_bin:
            path to MPI binary (e.g. mpirun, srun) if needed. Use only when running multi-node with gmx_mpi
        mpi_np:
            number of MPI ranks to start by mpi binary (e.g. mpirun -np N). Only used when mpi_bin is provided. Use only when running multi-node with gmx_mpi.
            When using srun, this parameter is not needed as srun automatically allocates the number of ranks specified by the job scheduler.
        num_threads_mpi:
            number of thread-MPI ranks to start by gmx mdrun (-ntmpi). Use only with gmx binary. (0 is guess)
        num_threads_omp:
            number of OpenMP threads per MPI rank to start by gmx mdrun (-ntomp) (0 is guess)
        use_gpu:
            whether to use GPU support or not
        restart:
            whether to restart the workflow from the last completed step or start from the beginning.
        forcefield: 
            forcefield to be used by pdb2gmx to generate the topology. Default: amber99sb-ildn. 
            See values supported by pdb2gmx (gromos45a3, charmm27, gromos53a6, amber96, amber99, 
            gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, amber99sb, amber03, amber99sb-ildn, 
            oplsaa, amber94, amber99sb-star-ildn-mut). 
        ions_concentration: 
            salt concentration to be used in the simulation in mols/L. Default: 0.15.
        temperature:
            temperature to be used in the simulation. Default: 300.0.
        random_seed:
            random seed to be used to generate velocities and control the temperature. Default: -1.
        setup_only: 
            whether to only setup the system or also run the simulations
        dt:
            time step to be used in the simulation. Default: 2 fs.
        equil_time:
            Time of each equilibration simulation in ns. Default: 1 ns.
        equil_frames:
            Number of frames to save during the equilibration steps. Default: 500 frames.
        equil_only: 
            whether to only run the equilibration or also run the production simulations
        prod_time: 
            Time of production simulation in ns. Default: 100 ns.
        prod_frames:
            Number of frames to save during the production step. Default: 2000 frames.
        debug:
            whether to run the workflow in debug mode (keep temporary files)
        output_path: 
            path to output folder

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

    # Create a default configuration file if needed
    config_args = {
        'gmx_bin': gmx_bin,
        'mpi_bin': mpi_bin,
        'mpi_np': mpi_np,
        'num_threads_mpi': num_threads_mpi,
        'num_threads_omp': num_threads_omp,
        'use_gpu': use_gpu,
        'restart': restart,
        'forcefield': forcefield,
        'ions_concentration': ions_concentration,
        'temp': temperature,
        'seed': random_seed,
        'dt': dt,
        'equil_time': equil_time,
        'equil_frames': equil_frames,
        'prod_time': prod_time,
        'prod_frames': prod_frames,
        'debug': debug
    }
    default_config = create_config_file(configuration_path, **config_args)
    if default_config:
        configuration_path = 'config.yml'  # Use the default config file if it was created
        
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Enforce output_path if provided
    if output_path is not None:
        output_path = fu.get_working_dir_path(output_path, restart = conf.properties.get('restart', 'False'))
        conf.working_dir_path = output_path
    else:
        output_path = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()
    
    # Check input files
    input_mode = check_inputs(input_pdb_path, 
                 input_gro_path, 
                 input_top_path, 
                 ligands_top_folder, 
                 input_tpr_path,
                 input_cpt_path,
                 input_ndx_path,
                 input_plumed_path,
                 input_plumed_folder,
                 setup_only,
                 equil_only,
                 global_log)
    
    # Get ligands information
    ligands_dict = get_ligands(ligands_top_folder, global_log)
    chains_dict = {}

    # If prepared structure is not provided
    setup_needed = input_mode == 'input_pdb'
    if setup_needed:
        
        ###########################################
        # Prepare topology and coordinates for MD #
        ###########################################
        
        setup_prefix = "step1_setup"
        setup_prop = conf.get_prop_dic(prefix=setup_prefix)
        setup_paths = conf.get_paths_dic(prefix=setup_prefix)

        setup_paths["step1_pdb2gmx"]["input_pdb_path"] = input_pdb_path
        setup_prop["step1_pdb2gmx"]["force_field"]=forcefield
        # NOTE: I don't see the HIS protonation states
        # Determine protonation state of titratable residues from resname in the input PDB file
        global_log.info(f"step1_pdb2gmx: Reading protonation states for titratable residues (0: deprotonated, 1: protonated)")
        for residue in gmx_titra_resnames.keys():
            setup_prop["step1_pdb2gmx"][residue.lower()] = read_protonation_states(input_pdb_path, residue, global_log)
            global_log.info(f"step1_pdb2gmx: {residue} protonation state: {setup_prop['step1_pdb2gmx'][residue.lower()]}")
        
        # STEP 1: add H atoms, generate coordinate (.gro) and topology (.top) file for the PDB structure
        global_log.info("step1_pdb2gmx: Generate the topology")
        pdb2gmx(**setup_paths["step1_pdb2gmx"], properties=setup_prop["step1_pdb2gmx"])
        
        master_index_file = ""
        chains_dict = get_chains_dict(input_pdb_path) 

        # STEP 2: Add chain groups to a master index file to be used for position restraints
        for chain_id in chains_dict:

            prefix = f"{setup_prefix}/step2_adding_chain_{chain_id}_group"
            chain_prop = conf.get_prop_dic(prefix=prefix)
            chain_paths = conf.get_paths_dic(prefix=prefix)
            
            structure_path = setup_paths["step1_pdb2gmx"]["output_gro_path"]
            
            # If not the first chain, we need to append the group to the master index file
            if master_index_file:
                chain_paths["step1_make_ref_group"]["input_ndx_path"] = master_index_file

            # Create index file for chain
            global_log.info(f"{chain_id} > Create index group for the chain")
            chain_paths["step1_make_ref_group"]["input_structure_path"] = structure_path
            chain_prop["step1_make_ref_group"]["selection"] = f"ri {chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
            make_ndx(**chain_paths["step1_make_ref_group"], properties=chain_prop["step1_make_ref_group"])
            
            global_log.info(f"{chain_id} > Create index group for the chain's C-alpha atoms")
            chain_paths["step2_make_rest_group"]["input_structure_path"] = structure_path
            chain_prop["step2_make_rest_group"]["selection"] = f"a CA & ri {chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
            make_ndx(**chain_paths["step2_make_rest_group"], properties=chain_prop["step2_make_rest_group"])  
            
            # Save group names for each chain
            chains_dict[chain_id]["reference_group"] = f"r_{chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
            chains_dict[chain_id]["restrain_group"] = f"CA_&_r_{chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
            
            # Save POSRES names for each chain
            chains_dict[chain_id]["posres_name"] = f"CHAIN_{chain_id}_POSRES"
            
            # Update master index file
            master_index_file = chain_paths["step2_make_rest_group"]["output_ndx_path"]

        # Reference, restraint and chain triplet list to add restraints to the topology with ndx2resttop
        ref_rest_chain_triplet_list = ", ".join([f"({chains_dict[chain]['reference_group']}, {chains_dict[chain]['restrain_group']}, {chain})" for chain in chains_dict])
        
        # POSRES names for each chain
        posres_names = " ".join([chains_dict[chain]["posres_name"] for chain in chains_dict])
        
        # STEP 3: Append position restraints to the topology file
        global_log.info(f"step3_append_posres: Append restraints to the topology")
        setup_prop["step3_append_posres"]["ref_rest_chain_triplet_list"] = ref_rest_chain_triplet_list
        setup_prop["step3_append_posres"]["posres_names"] = posres_names
        setup_paths["step3_append_posres"]["input_ndx_path"] = master_index_file
        ndx2resttop(**setup_paths["step3_append_posres"], properties=setup_prop["step3_append_posres"])
        
        if ligands_dict:
            
            # STEP 4: Convert gro of main structure to pdb - to concatenate with ligands
            global_log.info("step4_structure_pdb: Convert GRO to PDB")
            gmx_trjconv_str(**setup_paths["step4_structure_pdb"], properties=setup_prop["step4_structure_pdb"])
            
            complex_topology_path = setup_paths["step3_append_posres"]["output_top_zip_path"]
            complex_pdb_path = setup_paths["step4_structure_pdb"]["output_str_path"]
            
            # STEP 5: Add each ligand to the PDB and topology files
            for ligand_id in ligands_dict:

                prefix = f"{setup_prefix}/step5_add_ligand_{ligand_id}"
                ligand_prop = conf.get_prop_dic(prefix=prefix)
                ligand_paths = conf.get_paths_dic(prefix=prefix)
                
                ligand_gro_path = ligands_dict[ligand_id]["coordinates"]
                ligand_itp_path = ligands_dict[ligand_id]["topology"]
                
                ligands_dict[ligand_id]["posres_name"] = f"LIGAND_{ligand_id}_POSRES"
            
                # STEP 1: Convert ligand coordinates from gro to pdb
                global_log.info(f"{ligand_id} > Convert ligand GRO to PDB")
                ligand_paths["step1_create_ligand_pdb"]["input_structure_path"] = ligand_gro_path
                ligand_paths["step1_create_ligand_pdb"]["input_top_path"] = ligand_gro_path
                gmx_trjconv_str(**ligand_paths["step1_create_ligand_pdb"], properties=ligand_prop["step1_create_ligand_pdb"])
                
                # STEP 2: Create complex pdb file concatenating the current complex and the new ligand
                global_log.info(f"{ligand_id} > Create complex PDB file")
                ligand_paths["step2_create_complex_pdb"]["input_structure1"] = complex_pdb_path
                ligand_paths["step2_create_complex_pdb"]["output_structure_path"] = os.path.join(str(Path(ligand_paths["step2_create_complex_pdb"]["output_structure_path"]).parent), f"{ligand_id}_complex.pdb")
                cat_pdb(**ligand_paths["step2_create_complex_pdb"], properties=ligand_prop["step2_create_complex_pdb"])
                
                # Update complex pdb path for the next ligands
                complex_pdb_path = ligand_paths["step2_create_complex_pdb"]["output_structure_path"]
                
                # STEP 3: Make ndx file for the ligand's heavy atoms
                global_log.info(f"{ligand_id} > Create index file for the ligand's heavy atoms")
                ligand_paths["step3_make_ligand_ndx"]["input_structure_path"] = ligand_gro_path
                make_ndx(**ligand_paths["step3_make_ligand_ndx"], properties=ligand_prop["step3_make_ligand_ndx"])
                
                ligand_restraints_path = os.path.join(str(Path(ligand_paths["step4_create_ligand_restraints"]["output_itp_path"]).parent), f"{ligand_id}_posre.itp")
                
                # STEP 4: Generate restraints for the ligand's heavy atoms
                global_log.info(f"{ligand_id} > Generate restraints for ligand")
                ligand_paths["step4_create_ligand_restraints"]["input_structure_path"] = ligand_gro_path
                ligand_paths["step4_create_ligand_restraints"]["output_itp_path"] = ligand_restraints_path
                genrestr(**ligand_paths["step4_create_ligand_restraints"], properties=ligand_prop["step4_create_ligand_restraints"])
                
                # STEP 5: Append parameterized ligand to the current complex topology zip file
                ligand_paths["step5_append_ligand_topology"]["input_top_zip_path"] = complex_topology_path
                ligand_paths["step5_append_ligand_topology"]["input_itp_path"] = ligand_itp_path
                ligand_paths["step5_append_ligand_topology"]["input_posres_itp_path"] = ligand_restraints_path
                ligand_prop["step5_append_ligand_topolofgy"]["posres_name"] = ligands_dict[ligand_id]["posres_name"]
                global_log.info(f"{ligand_id} > Append ligand to the topology")
                append_ligand(**ligand_paths["step5_append_ligand_topology"], properties=ligand_prop["step5_append_ligand_topology"])
                
                # Update complex topology path for the next ligands
                complex_topology_path = ligand_paths["step5_append_ligand_topology"]["output_top_zip_path"]
                
            # Modify paths for the next steps
            setup_paths["step6_editconf"]["input_gro_path"] = complex_pdb_path
            setup_paths["step7_solvate"]["input_top_zip_path"] = complex_topology_path
            
        # STEP 6: Create simulation box
        global_log.info("step6_editconf: Create the solvent box")
        editconf(**setup_paths["step6_editconf"], properties=setup_prop["step6_editconf"])

        # STEP 7: Add solvent molecules
        global_log.info("step7_solvate: Fill the solvent box with water molecules")
        solvate(**setup_paths["step7_solvate"], properties=setup_prop["step7_solvate"])

        # STEP 8: ion generation pre-processing
        global_log.info("step8_grompp_genion: Preprocess ion generation")
        grompp(**setup_paths["step8_grompp_genion"], properties=setup_prop["step8_grompp_genion"])

        # STEP 9: ion generation
        global_log.info("step9_genion: Ion generation")
        setup_prop["step9_genion"]["concentration"] = ions_concentration
        genion(**setup_paths["step9_genion"], properties=setup_prop["step9_genion"])

        if setup_only:
            global_log.info("Set up only: setup_only flag is set to True! Exiting...")
            return
        
    equil_needed = (input_mode == 'input_pdb') or (input_mode == 'prepared_system') 
    if equil_needed:
        
        #######################################
        # Minimize and equilibrate the system #
        #######################################
        
        equilibrate_prefix = "step2_equil"
        equil_prop = conf.get_prop_dic(prefix=equilibrate_prefix)
        equil_paths = conf.get_paths_dic(prefix=equilibrate_prefix)

        if input_mode == 'input_pdb':
            input_gro_path = setup_paths["step9_genion"]["output_gro_path"]
            input_top_path = setup_paths["step9_genion"]["output_top_zip_path"]
            
        # Connect equil steps with previous ones
        equil_paths['step1_grompp_min']['input_gro_path'] = input_gro_path
        equil_paths['step1_grompp_min']['input_top_zip_path'] = input_top_path
        equil_paths['step5_grompp_nvt']['input_top_zip_path'] = input_top_path
        equil_paths['step8_grompp_npt']['input_top_zip_path'] = input_top_path
    
        # STEP 1: minimization pre-processing
        global_log.info("step1_grompp_min: Preprocess energy minimization")
        grompp(**equil_paths["step1_grompp_min"], properties=equil_prop["step1_grompp_min"])
        
        # STEP 2: minimization
        global_log.info("step2_mdrun_min: Execute energy minimization")
        mdrun(**equil_paths["step2_mdrun_min"], properties=equil_prop["step2_mdrun_min"])

        # STEP 3: create index file
        global_log.info("step3_make_ndx: Create index file")
        make_ndx(**equil_paths["step3_make_ndx"], properties=equil_prop["step3_make_ndx"])

        # STEP 4: dump potential energy evolution
        global_log.info("step4_energy_min: Compute potential energy during minimization")
        gmx_energy(**equil_paths["step4_energy_min"], properties=equil_prop["step4_energy_min"])

        # Prepare position restraints definition for equilibration steps
        if chains_dict:
            chain_posres = " ".join([f"-D{chains_dict[chain]['posres_name']}" for chain in chains_dict])
        else:
            chain_posres = ""
        if ligands_dict:
            ligand_posres = " ".join([f"-D{ligands_dict[ligand]['posres_name']}" for ligand in ligands_dict])
        else:
            ligand_posres = ""
        eq_posres = f"{chain_posres} {ligand_posres}"
        
        # STEP 5: NVT equilibration pre-processing
        # Add position restraints if any
        if eq_posres != " ":
            equil_prop["step5_grompp_nvt"]["mdp"]["define"] = eq_posres
        global_log.info("step5_grompp_nvt: Preprocess NVT equilibration")
        grompp(**equil_paths["step5_grompp_nvt"], properties=equil_prop["step5_grompp_nvt"])

        # STEP 6: NVT equilibration
        global_log.info("step6_mdrun_nvt: Execute NVT equilibration")
        mdrun(**equil_paths["step6_mdrun_nvt"], properties=equil_prop["step6_mdrun_nvt"])

        # STEP 7: dump temperature evolution
        global_log.info("step7_temp_nvt: Compute temperature during NVT equilibration")
        gmx_energy(**equil_paths["step7_temp_nvt"], properties=equil_prop["step7_temp_nvt"])

        # STEP 8: NPT equilibration pre-processing
        # Add position restraints if any
        if eq_posres != " ":
            equil_prop["step8_grompp_npt"]["mdp"]["define"] = eq_posres
        global_log.info("step8_grompp_npt: Preprocess NPT equilibration")
        grompp(**equil_paths["step8_grompp_npt"], properties=equil_prop["step8_grompp_npt"])

        # STEP 9: NPT equilibration
        global_log.info("step9_mdrun_npt: Execute NPT equilibration")
        mdrun(**equil_paths["step9_mdrun_npt"], properties=equil_prop["step9_mdrun_npt"])

        # STEP 10: dump density and pressure evolution
        global_log.info("step10_density_npt: Compute Density & Pressure during NPT equilibration")
        gmx_energy(**equil_paths["step10_density_npt"], properties=equil_prop["step10_density_npt"])
        
        # NOTE: add free equilibration removing those restraints that are not needed in the production run - none by default

        if equil_only:
            global_log.info("Equilibration only: equil_only flag is set to True! Exiting...")
            return
    
    ##########################
    # Production simulations #
    ##########################

    production_prefix = "step3_prod"
    prod_prop = conf.get_prop_dic(prefix=production_prefix)
    prod_paths = conf.get_paths_dic(prefix=production_prefix)
    
    # Connect production steps with previous ones
    # STEP 1: Prepare the production run
    if input_mode == 'restart_simulation':           # NOTE: remember to include -noapend when restarting from a cpt file!
                                                     # NOTE: also remember to adapt the paths to the traj files from the plugin - names will change with replicas
        # Extend the simulation time
        prod_paths['step1B_convert_tpr']['input_tpr_path'] = input_tpr_path
        global_log.info("step1B_convert_tpr: Extend the simulation time of the input TPR file")
        convert_tpr(**prod_paths['step1B_convert_tpr'], properties=prod_prop['step1B_convert_tpr'])

        prod_paths['step2_mdrun_prod']['input_tpr_path'] = prod_paths['step1B_convert_tpr']['output_tpr_path']
        prod_paths['step2_mdrun_prod']['input_cpt_path'] = input_cpt_path
        prod_prop['step2_mdrun_prod']['noappend'] = True
    else:
        prod_paths['step1_grompp_md']['input_gro_path'] = equil_paths["step9_mdrun_npt"]['output_gro_path']
        prod_paths['step1_grompp_md']['input_cpt_path'] = equil_paths["step9_mdrun_npt"]['output_cpt_path']
        prod_paths['step1_grompp_md']['input_top_zip_path'] = input_top_path
        prod_paths['step1_grompp_md']['input_ndx_path'] = equil_paths["step3_make_ndx"]['output_ndx_path']
        input_ndx_path = equil_paths["step3_make_ndx"]['output_ndx_path']
        
        # Modify default temperature coupling groups
        prod_prop["step1_grompp_md"]["mdp"]["define"] = "" # NOTE: here restraint what is asked by the user
        global_log.info("step1_grompp_md: Preprocess production simulation")
        grompp(**prod_paths['step1_grompp_md'], properties=prod_prop["step1_grompp_md"])
        
        input_tpr_path = prod_paths['step1_grompp_md']['output_tpr_path']

    # STEP 2: free NPT production run
    prod_paths['step2_mdrun_prod']['input_plumed_path'] = input_plumed_path
    prod_paths['step2_mdrun_prod']['input_plumed_folder'] = input_plumed_folder
    prod_paths['step2_mdrun_prod']['output_plumed_folder'] = os.path.join(prod_prop['step2_mdrun_prod']['path'], 'plumed_outputs')
    global_log.info("step2_mdrun_prod: Execute production simulation")
    mdrun(**prod_paths['step2_mdrun_prod'], properties=prod_prop['step2_mdrun_prod'])
    
    ############################
    # Post-processing analysis #
    ############################
    
    # NOTE: Make these steps fault-tolerant
    # NOTE: Change the selections to make them more general
    # NOTE: Improve the centering and fitting if needed using biopython to find specific atoms
    
    analysis_prefix = "step4_analysis"
    analysis_prop = conf.get_prop_dic(prefix=analysis_prefix)
    analysis_paths = conf.get_paths_dic(prefix=analysis_prefix)

    analysis_paths['step1_gro2pdb']['input_top_path'] = prod_paths["step2_mdrun_prod"]['output_gro_path']
    analysis_paths['step1_gro2pdb']['input_structure_path'] = prod_paths["step2_mdrun_prod"]['output_gro_path']
    analysis_paths['step2_rmsd_equilibrated']['input_traj_path'] = prod_paths["step2_mdrun_prod"]['output_xtc_path']
    analysis_paths['step3_rmsd_experimental']['input_traj_path'] = prod_paths["step2_mdrun_prod"]['output_xtc_path']
    analysis_paths['step4_rgyr']['input_traj_path'] = prod_paths["step2_mdrun_prod"]['output_xtc_path']
    analysis_paths['step5_rmsf']['input_traj_path'] = prod_paths["step2_mdrun_prod"]['output_xtc_path']
    analysis_paths['step6_dry_str']['input_top_path'] = input_tpr_path
    analysis_paths['step7_dry_traj']['input_top_path'] = input_tpr_path
    analysis_paths['step7_dry_traj']['input_traj_path'] = prod_paths["step2_mdrun_prod"]['output_xtc_path']
    analysis_paths['step8_center']['input_top_path'] = input_tpr_path
    analysis_paths['step9_image_traj']['input_top_path'] = input_tpr_path
    analysis_paths['step10_fit_traj']['input_top_path'] = input_tpr_path

    if input_ndx_path:
        analysis_paths['step6_dry_str']['input_index_path'] = input_ndx_path
        analysis_paths['step7_dry_traj']['input_index_path'] = input_ndx_path
        analysis_paths['step8_center']['input_index_path'] = input_ndx_path
        analysis_paths['step9_image_traj']['input_index_path'] = input_ndx_path
        analysis_paths['step10_fit_traj']['input_index_path'] = input_ndx_path
    
    # STEP 1: conversion of topology from gro to pdb
    global_log.info("step1_gro2pdb: Convert topology from GRO to PDB")
    gmx_trjconv_str(**analysis_paths["step1_gro2pdb"], properties=analysis_prop["step1_gro2pdb"])

    # STEP 2: compute the RMSD with respect to equilibrated structure
    global_log.info("step2_rmsd_equilibrated: Compute Root Mean Square deviation against equilibrated structure")
    gmx_rms(**analysis_paths['step2_rmsd_equilibrated'], properties=analysis_prop['step2_rmsd_equilibrated'])
    
    # STEP 3: compute the RMSD with respect to minimized structure
    global_log.info("step3_rmsd_experimental: Compute Root Mean Square deviation against minimized structure (exp)")
    gmx_rms(**analysis_paths['step3_rmsd_experimental'], properties=analysis_prop['step3_rmsd_experimental'])

    # STEP 4: compute the Radius of gyration
    global_log.info("step4_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
    gmx_rgyr(**analysis_paths['step4_rgyr'], properties=analysis_prop['step4_rgyr'])

    # STEP 5: compute the RMSF
    global_log.info("step5_rmsf: Compute Root Mean Square Fluctuation to measure the protein flexibility during the free MD simulation")
    cpptraj_rmsf(**analysis_paths['step5_rmsf'], properties=analysis_prop['step5_rmsf'])

    # STEP 6: obtain dry structure
    try:
        global_log.info("step6_dry_str: Obtain dry structure")
        gmx_trjconv_str(**analysis_paths["step6_dry_str"], properties=analysis_prop["step6_dry_str"])
        global_log.info("step6_dry_str: Completed successfully")
    
    except SystemExit as e:
        global_log.error(
            f"step6_dry_str failed due to GMX exit "
            f"(SystemExit, code={e.code})"
        )
    except Exception:
        global_log.exception("step6_dry_str failed with unexpected exception")

    # STEPS 7-10: process trajectory: dry, center, image, fit
    try:
        # STEP 7: obtain dry trajectory
        global_log.info("step7_dry_traj: Obtain dry trajectory")
        gmx_trjconv_trj(**analysis_paths["step7_dry_traj"], properties=analysis_prop["step7_dry_traj"])

        # STEP 8: center the trajectory
        global_log.info("step8_center: Center the trajectory")
        gmx_image(**analysis_paths['step8_center'], properties=analysis_prop['step8_center'])

        # Remove intermediate trajectory
        os.remove(analysis_paths["step7_dry_traj"]["output_traj_path"])

        # STEP 9: image the trajectory
        global_log.info("step9_image_traj: Imaging the trajectory")
        gmx_image(**analysis_paths['step9_image_traj'], properties=analysis_prop['step9_image_traj'])

        # Remove intermediate trajectory
        os.remove(analysis_paths["step8_center"]["output_traj_path"])
        
        # STEP 10: fit the trajectory
        global_log.info("step10_fit_traj: Fit the trajectory")
        gmx_image(**analysis_paths['step10_fit_traj'], properties=analysis_prop['step10_fit_traj'])

        # Remove intermediate trajectory
        os.remove(analysis_paths["step9_image_traj"]["output_traj_path"])

    except SystemExit as e:
        global_log.error(
            f"steps 7 to 10 failed due to GMX exit "
            f"(SystemExit, code={e.code})"
        )
    except Exception:
        global_log.exception("steps 7 to 10 failed with unexpected exception")

    if default_config:
        # Move the default configuration file to the output path
        shutil.move(configuration_path, os.path.join(output_path, 'config.yml'))
        configuration_path = os.path.join(output_path, 'config.yml')
        
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

    parser = argparse.ArgumentParser("MD Simulation with GROMACS")

    ###############
    # Input files #
    ###############
    parser.add_argument('--input_pdb', dest='input_pdb_path', type=str,
                        help="""Input PDB file. The workflow assumes the protonation state specified by the residue names
                        is the correct one. Default: None""",
                        required=False)

    parser.add_argument('--ligands_folder', dest='ligands_top_folder', type=str,
                        help="""Path to folder with .itp and .gro files for the ligands that 
                        should be included in the simulation. Make sure the coordinates of the 
                        ligands correspond to the PDB used. Only compatible with '--input_pdb'. Default: None""",
                        required=False)

    parser.add_argument('--input_gro', dest='input_gro_path', type=str,
                        help="""Input structure file ready to minimize (.gro). To provide an externally prepared system, 
                        use together with '--input_top'. Default: None""",
                        required=False)

    parser.add_argument('--input_top', dest='input_top_path', type=str,
                        help="""Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, 
                        use together with '--input_gro'. Default: None""",
                        required=False)
    
    parser.add_argument('--input_tpr', dest='input_tpr_path', type=str,
                        help="""Input portable binary run input file (.tpr) to restart a simulation. Use together with '--input_cpt'.
                        Default: None""",
                        required=False)
    
    parser.add_argument('--input_cpt', dest='input_cpt_path', type=str,
                        help="""Input checkpoint file (.cpt) to restart a simulation. Use together with '--input_tpr'.
                        Default: None""",
                        required=False)

    parser.add_argument('--input_ndx', dest='input_ndx_path', type=str,
                        help="""Input index file (.ndx) to use in the simulation. If not provided no index file will be used
                        for the analysis of the trajectories and only standard groups will be available.
                        Default: None""",
                        required=False)

    parser.add_argument('--input_plumed_path', dest='input_plumed_path', type=str,
                        help="""Path to the main PLUMED input file (plumed.dat). If provided, PLUMED will be used during the production run.
                        Default: None""",
                        required=False)
    
    parser.add_argument('--input_plumed_folder', dest='input_plumed_folder', type=str,
                        help="""Path to the folder with all files needed by the main PLUMED input file, see input_plumed_path.
                        Default: None""",
                        required=False)

    # NOTE: Add an option to remove raw data and just keep prepared traj
    # NOTE: Add an option to leave waters and ions in the prepared traj and top
    # NOTE: Add option for H mass repartitioning
    # NOTE: Add flag to determine what should remain restrained during the production run - currently everything is free always
    # NOTE: Add progressive release of position restraints during equilibration (additional steps if needed)

    #########################
    # Configuration options #
    #########################

    # General configuration file
    parser.add_argument('--config', dest='config_path', type=str,
                        help="Configuration file (YAML)",
                        required=False)

    parser.add_argument('--gmx_bin', dest='gmx_bin', type=str,
                        help="Path to GROMACS binary (gmx for single node and gmx_mpi for multi-node). Default: gmx",
                        required=False, default='gmx')
    
    parser.add_argument('--mpi_bin', dest='mpi_bin', type=str,
                        help="Path to MPI binary. Default: null",
                        required=False, default='null')
    
    parser.add_argument('--mpi_np', dest='mpi_np', type=int,
                        help="Number of MPI processes given to the mpi_bin. Default: None",
                        required=False)
    
    parser.add_argument('--num_threads_mpi', dest='num_threads_mpi', type=int,
                        help="Number of MPI threads. Default: 0 (Let GROMACS guess)",
                        required=False, default=0)
    
    parser.add_argument('--num_threads_omp', dest='num_threads_omp', type=int,
                        help="Number of OpenMP threads. Default: 0 (Let GROMACS guess)",
                        required=False, default=0)

    parser.add_argument('--use_gpu', action='store_true',
                        help="""Calculate non-bonding interactions and particle-mesh ewald in GPU by adding '-nb gpu -pme gpu' 
                        to mdrun call. If not used, gmx will still use a GPU for these calculations if available. Default: False""",
                        required=False, default=False)
    
    parser.add_argument('--restart', action='store_true',
                        help="Restart the workflow from the last completed step. Default: False",
                        required=False, default=False)

    parser.add_argument('--forcefield', dest='forcefield', type=str,
                        help="Forcefield to use. Default: amber99sb-ildn",
                        required=False, default='amber99sb-ildn')

    parser.add_argument('--ions_concentration', dest='ions_concentration', type=float,
                        help="Concentration of ions in the system in mol/L (M). Default: 0.15 M",
                        required=False, default=0.15)
    
    parser.add_argument('--temp', dest='temperature', type=float,
                        help="Temperature of the system in K. Default: 300",
                        required=False, default=300)
    
    parser.add_argument('--seed', dest='random_seed', type=int,
                        help="Random seed for the simulations. If given, new velocities will be generated with this seed. Default: -1",
                        required=False, default=-1)
    
    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--dt', dest='dt', type=float,
                        help="Time step in fs. Default: 2 fs",
                        required=False, default=2)
    
    parser.add_argument('--equil_time', dest='equil_time', type=float,
                        help="Time of each equilibration step in ns. Default: 1.0 ns",
                        required=False, default=1.0)

    parser.add_argument('--equil_frames', dest='equil_frames', type=int,
                        help="Number of frames to save during the equilibration steps. Default: 500 frames",
                        required=False, default=500)
    
    parser.add_argument('--equil_only', action='store_true',
                        help="Only run the equilibration steps. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--prod_time', dest='prod_time', type=float,
                        help="Total time of the production simulation in ns. Default: 100.0 ns",
                        required=False, default=100.0)

    parser.add_argument('--prod_frames', dest='prod_frames', type=int,
                        help="Number of frames to save during the production steps. Default: 2000 frames",
                        required=False, default=2000)
    
    parser.add_argument('--debug', action='store_true',
                        help="Activate debug mode with more verbose logging. Default: False",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path', type=str,
                        help="Output path. Default: 'output' in the current working directory",
                        required=False, default='output')

    args = parser.parse_args()

    # Convert to corresponding types
    if args.ions_concentration:
        args.ions_concentration = float(args.ions_concentration)
    if args.temperature:
        args.temperature = float(args.temperature)
    if args.prod_time:
        args.prod_time = float(args.prod_time)
    if args.equil_time:
        args.equil_time = float(args.equil_time)
    if args.dt:
        args.dt = float(args.dt)
    if args.equil_frames:
        args.equil_frames = int(args.equil_frames)
    if args.prod_frames:
        args.prod_frames = int(args.prod_frames)
    if args.random_seed:
        args.random_seed = int(args.random_seed)
        
    # Run the main workflow
    main_wf(input_pdb_path=args.input_pdb_path, 
            ligands_top_folder=args.ligands_top_folder, 
            input_gro_path=args.input_gro_path, 
            input_top_path=args.input_top_path,
            input_tpr_path=args.input_tpr_path,
            input_cpt_path=args.input_cpt_path,
            input_ndx_path=args.input_ndx_path,
            input_plumed_path=args.input_plumed_path,
            input_plumed_folder=args.input_plumed_folder,
            configuration_path=args.config_path, 
            gmx_bin=args.gmx_bin,
            mpi_bin=args.mpi_bin,
            mpi_np=args.mpi_np,
            num_threads_mpi=args.num_threads_mpi,
            num_threads_omp=args.num_threads_omp,
            use_gpu=args.use_gpu,
            restart=args.restart,
            forcefield=args.forcefield, 
            ions_concentration=args.ions_concentration,
            temperature=args.temperature,
            random_seed=args.random_seed,
            setup_only=args.setup_only, 
            dt=args.dt, 
            equil_time=args.equil_time,
            equil_frames=args.equil_frames,
            equil_only=args.equil_only, 
            prod_time=args.prod_time, 
            prod_frames=args.prod_frames,
            debug=args.debug,
            output_path=args.output_path)