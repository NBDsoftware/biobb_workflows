#!/usr/bin/env python3

# Importing all the needed libraries
from typing import List, Dict, Union, Optional
from pathlib import Path
import argparse
import random
import shutil
import time
import os

from Bio.PDB import PDBParser

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
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
    global_log
    ) -> None:
    """
    Check the inputs for the workflow. The function checks if either the input PDB or the gro and top files
    exist. If the ligands_top_folder is provided, it checks if the folder exists and if it's not empty. 
    If any of the checks fail, an error is raised.
    
    Inputs
    ------
    
        input_pdb_path (str): Path to the input PDB file.
        input_gro_path (str): Path to the input GRO file.
        input_top_path (str): Path to the input topology file.
        ligands_top_folder (str): Path to the folder with the ligands .itp files.
        global_log: Logger object for logging messages.
    """
    
    # Check if the input PDB file exists
    if input_pdb_path is not None:
        if not os.path.exists(input_pdb_path):
            global_log.error(f"Input PDB file {input_pdb_path} not found")
            raise FileNotFoundError(f"Input PDB file {input_pdb_path} not found")
    
    # Check if the input GRO file exists
    if input_gro_path is not None:
        if not os.path.exists(input_gro_path):
            global_log.error(f"Input GRO file {input_gro_path} not found")
            raise FileNotFoundError(f"Input GRO file {input_gro_path} not found")
    
    # Check if the input topology file exists
    if input_top_path is not None:
        if not os.path.exists(input_top_path):
            global_log.error(f"Input topology file {input_top_path} not found")
            raise FileNotFoundError(f"Input topology file {input_top_path} not found")
    
    # Check if the user is not providing both input PDB and GRO/Top files
    if input_pdb_path is not None and (input_gro_path is not None or input_top_path is not None):
        global_log.error("You cannot provide both input PDB and GRO/Top files. Please provide only one of them.")
        raise ValueError("You cannot provide both input PDB and GRO/Top files. Please provide only one of them.")

    # Check if the user is not providing any input PDB or GRO/Top files
    if input_pdb_path is None and (input_gro_path is None or input_top_path is None):
        global_log.error("You must provide either an input PDB file or both input GRO and Top files.")
        raise ValueError("You must provide either an input PDB file or both input GRO and Top files.")

    # Check if the ligands topology folder exists
    if ligands_top_folder is not None:
        if not os.path.exists(ligands_top_folder):
            global_log.error(f"Ligands topology folder {ligands_top_folder} not found")
            raise FileNotFoundError(f"Ligands topology folder {ligands_top_folder} not found")
        
        # Check if the ligands topology folder is empty
        if not os.listdir(ligands_top_folder):
            global_log.error(f"Ligands topology folder {ligands_top_folder} is empty")
            raise ValueError(f"Ligands topology folder {ligands_top_folder} is empty")

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
    
    HETEROATOM chains are not considered. 
    
    Parameters
    ----------
    
    pdb_file : str
        Path to the PDB file.
    
    Returns
    -------
    
    Dict:
        A dictionary with the chain IDs as keys and a list with the residue numbers as values.
        
        Example:
        
        {
            'A': 
                {
                    'residues': [1, 100]
                },
            'B': 
                {
                    'residues': [1, 200]
                }
        }
    """
    
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    # Initialize the dictionary with the chain IDs and residue numbers
    chains_dict = {}
    
    # Iterate over all residues in the structure
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            residue_numbers = [residue.get_id()[1] for residue in chain]
            chains_dict[chain_id] = {'residues':[min(residue_numbers), max(residue_numbers)]}
    
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
    
# Process analysis files
def merge_xvgtimeseries_files(file_paths: list, output_xvg_path: str) -> None:
    """
    Merges multiple XVG files containing time series data into a single file.
    
    The function reads each file, handles time data by adjusting the time
    of subsequent files, ensures no duplicated points at the boundary between files,
    and writes the merged data into a new XVG file. It also preserves the header 
    from the first input file.
    
    Inputs
    ------
    filepaths :
        A list of file paths for the input XVG files to be merged. 
        Each file should follow the format of time series data with comments starting 
        with '#' followed by a header section starting with '@' with plot details 
    
    output_filepath :
        The path for the output XVG file that will contain the merged data.

    Returns:
    -------
    None
        The merged file is written to the specified output filepath.

    Notes:
    -----
    - The function assumes that the time in each file starts at zero. 
    - The last data point of each file is compared to the first data point of the next file 
      to avoid duplication of the overlapping entry.
    - The function adjusts the time values of the second and subsequent files by adding 
      the final time value of the previous file to maintain a continuous time series.
    
    Example:
    -------
    To merge two files `part1.xvg` and `part2.xvg` into `merged_output.xvg`, you would do:

    >>> filepaths = ["part1.xvg", "part2.xvg"]
    >>> output_filepath = "merged_output.xvg"
    >>> merge_xvg_files(filepaths, output_filepath)
    
    After execution, a merged XVG file with a continuous time series will be saved at the
    specified output path.
    """

    # Initialize the time offset
    time_offset = 0
    header_written = False
    merged_data = []

    for i, file_path in enumerate(file_paths):
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Temporary variables for the current file
        file_data = []
        file_last_time = 0
        file_header = []

        # Read current file
        for line in lines:
            line = line.strip()

            if line.startswith('#'):
                # Ignore comments
                continue
            elif line.startswith('@'):
                # Store headers only from the first file
                if not header_written:
                    file_header.append(line)
            else:
                # Parse data (assuming it's space-separated values)
                data_parts = line.split()
                time = float(data_parts[0])
                values = [float(val) for val in data_parts[1:]]

                # Adjust time based on the current offset
                adjusted_time = time + time_offset

                # Add the data point to the current file data
                if len(merged_data) > 0 and len(file_data) == 0:
                    # If we are adding the first point of the next file, check for duplicates of the last point in the merged data
                    last_time_in_merged = merged_data[-1][0]
                    if adjusted_time <= last_time_in_merged:
                        # Skip this point as it is a duplicate
                        continue
                    
                # Add the data point to the current file data
                file_data.append([adjusted_time] + values)

        # If the file contains data, adjust the time offset and prepare for merging
        if len(file_data) > 0:
            file_last_time = file_data[-1][0]
            merged_data.extend(file_data)
            time_offset = file_last_time

        # Write header only once
        if not header_written and len(file_header) > 0:
            with open(output_xvg_path, 'w') as output_file:
                output_file.write('\n'.join(file_header) + '\n')
            header_written = True

    # Write the merged data to the output file
    with open(output_xvg_path, 'a') as output_file:
        for data_row in merged_data:
            output_file.write(" ".join(f"{val:.6f}" for val in data_row) + '\n')
    
def concatenate_gmx_analysis(conf, simulation_folders: List[str], output_path: str) -> None:
    """
    Concatenates the analysis files for each step of the GROMACS analysis including
    RMSD and Rgyr. The function reads the analysis files for each step from the
    simulation folders, merges them into a single file, and writes the merged data
    
    Inputs:
    -------
    
    conf :
        Configuration file reader object
    simulation_folders : 
        List of folder names containing the simulation data for each part or replica
    output_path :
        Path to the output folder where the concatenated files will be saved
    """
    
    # Construct a dictionary with the paths to the analysis files
    gmx_analysis_steps = ['step6A_rmsd_equilibrated', 'step6B_rmsd_experimental', 'step6C_rgyr']
    analysis_files = {step: [] for step in gmx_analysis_steps}        
    for simulation in simulation_folders:
        traj_paths = conf.get_paths_dic(prefix=simulation)
        for step in gmx_analysis_steps:
            analysis_files[step].append(traj_paths[step]['output_xvg_path'])
        
    # Concatenate the analysis files for each step
    for step, file_paths in analysis_files.items():
        
        # Output path for the concatenated file
        output_xvg_path = os.path.join(output_path, f"{step}.xvg")
        merge_xvgtimeseries_files(file_paths, output_xvg_path)

# YML construction
def config_contents(
    gmx_binary_path: str = 'gmx',
    mpi_bin: Optional[str] = None,
    num_threads_mpi: Optional[int] = None,
    num_threads_omp: Optional[int] = None,
    use_gpu: bool = False,
    restart: bool = False,
    forcefield: str = 'amber99sb-ildn',
    salt_conc: float = 0.15,
    temp: float = 300.0,
    dt: float = 0.002,
    equil_nsteps: int = 500000,
    equil_traj_freq: Optional[int] = None,
    nsteps: int = 50000000,
    traj_freq: Optional[int] = None
    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.
    
    The YAML file contains the configuration for the protein preparation workflow.
    
    Parameters
    ----------
    
    gmx_binary_path : str
        Path to the GROMACS binary. Default is 'gmx'.
    mpi_bin : str
        Path to the MPI binary. Default is 'null'.
    num_threads_mpi : int
        Number of threads for MPI. Default is None.
    num_threads_omp : int
        Number of threads for OpenMP. Default is None.
    use_gpu : bool
        Whether to use GPU for GROMACS. Default is False.
    restart : bool
        Whether to skip steps already performed. Default is True.
    forcefield : str
        Force field to use. Default is 'amber99sb-ildn'.
    salt_conc : float
        Concentration of salt in mols/L. Default is 0.15.
    temp : float
        Temperature in Kelvin. Default is 300.0.
    dt : float
        Time step in picoseconds. Default is 0.002.
    equil_nsteps : int
        Number of steps for equilibration. Default is 500000.
    equil_traj_freq : int
        Saving frequency of the trajectory during equilibration in number of steps. 
        Default is equil_nsteps // 500
    nsteps : int
        Number of steps for production. Default is 50000000
    traj_freq : int
        Saving frequency of the trajectory during production in number of steps.
        Default is nsteps // 1000
    
    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    if equil_traj_freq is None:
        equil_traj_freq = equil_nsteps // 500
    
    if traj_freq is None:
        traj_freq = nsteps // 1000
    
    if mpi_bin is None:
        mpi_bin = 'null'
        
    if num_threads_mpi is None:
        num_threads_mpi_config = ''
        mpi_np_config = ''
    else:
        num_threads_mpi_config = f"num_threads_mpi: {num_threads_mpi}"
        mpi_np_config = f"mpi_np: {num_threads_mpi}"

    if num_threads_omp is None:
        num_threads_omp_config = ''
    else:
        num_threads_omp_config = f"num_threads_omp: {num_threads_omp}"
    
    return f""" 
# Global properties (common for all steps)
global_properties:
  working_dir_path: output                                          # Workflow default output directory
  can_write_console_log: False                                      # Verbose writing of log information
  restart: {restart}                                                # Skip steps already performed
  remove_tmp: True                                                  # Remove temporal files

##################################################################
# Section 3 (Steps A-I): Prepare topology and coordinates for MD #
##################################################################

# Generate the topology of the structure with pdb2gmx
step3B_structure_topology:
  tool: pdb2gmx
  paths:
    input_pdb_path: path/to/input.pdb                     # Will be set by the workflow
    output_gro_path: structure.gro
    output_top_zip_path: structure_top.zip
  properties:
    binary_path: {gmx_binary_path}    # GROMACS binary path
    force_field: {forcefield}     # Will be set by the workflow
    water_type: tip3p             # spc, spce, tip3p, tip4p, tip5p, tips3p
    ignh: False                   # Ignore hydrogens in the input structure
    merge: False                  # Merge all chains into one molecule

# Add the reference group to the index file
step3C_make_ref_group:
  tool: make_ndx
  paths:
    input_structure_path: dependency/step3B_structure_topology/output_gro_path
    output_ndx_path: chain.ndx
  properties:
    binary_path: {gmx_binary_path}   # GROMACS binary path
    selection:  "System"         # Will be set by the workflow

# Add the restrained group to the index file
step3C_make_rest_group:
  tool: make_ndx
  paths:
    input_structure_path: dependency/step3B_structure_topology/output_gro_path
    input_ndx_path: dependency/step3C_make_ref_group/output_ndx_path
    output_ndx_path: calpha.ndx
  properties:
    binary_path: {gmx_binary_path}   # GROMACS binary path
    selection: "a CA"            # Will be set by the workflow

# Append position restraints to the topology file using the reference and restrained groups of the index file
step3D_append_posres:
  tool: ndx2resttop
  paths:
    input_ndx_path: dependency/step3C_make_rest_group/output_ndx_path
    input_top_zip_path: dependency/step3B_structure_topology/output_top_zip_path
    output_top_zip_path: structure_top.zip
  properties:
    force_constants: "500 500 500"

step3E_structure_pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: dependency/step3B_structure_topology/output_gro_path
    input_structure_path: dependency/step3B_structure_topology/output_gro_path
    output_str_path: structure.pdb
  properties:
    binary_path: {gmx_binary_path}   # GROMACS binary path

step3F_ligand_pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: path/to/ligand.gro        # Will be set by the workflow
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    output_str_path: ligand.pdb
  properties:   
    binary_path: {gmx_binary_path}   # GROMACS binary path

step3G_complex_pdb:
  tool: cat_pdb
  paths: 
    input_structure1: dependency/step3E_structure_pdb/output_str_path  # Will be set by the workflow
    input_structure2: dependency/step3F_ligand_pdb/output_str_path
    output_structure_path: complex.pdb

step3H_make_ligand_ndx:
  tool: make_ndx
  paths:
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    output_ndx_path: ligand_heavy_atoms.ndx
  properties:
    binary_path: {gmx_binary_path}  # GROMACS binary path
    selection: "0 & ! a H*"

step3I_ligand_restraints:
  tool: genrestr
  paths:
    input_structure_path: path/to/ligand.gro  # Will be set by the workflow
    input_ndx_path: dependency/step3H_make_ligand_ndx/output_ndx_path
    output_itp_path: ligand_restraints.itp    # Will be set by the workflow
  properties:
    binary_path: {gmx_binary_path} 
    restrained_group: "System_&_!H*"
    force_constants: "500 500 500"

step3J_append_ligand:
  tool: append_ligand
  paths:
    input_top_zip_path: dependency/step3D_append_posres/output_top_zip_path  # Will be set by the workflow
    input_itp_path: path/to/ligand.itp                                            # Will be set by the workflow
    input_posres_itp_path: dependency/step3I_ligand_restraints/output_itp_path    # Path to ligand position restraint topology file
    output_top_zip_path: complex_top.zip

step3K_editconf:
  tool: editconf
  paths:
    input_gro_path: dependency/step3B_structure_topology/output_gro_path
    output_gro_path: editconf.gro
  properties:
    binary_path: {gmx_binary_path}   # GROMACS binary path
    box_type: octahedron         # cubic, triclinic, octahedron, dodecahedron
    distance_to_molecule: 1.0    # Distance of the box from the outermost atom in nm

step3L_solvate:
  tool: solvate
  paths:
    input_solute_gro_path: dependency/step3K_editconf/output_gro_path
    input_top_zip_path: dependency/step3D_append_posres/output_top_zip_path
    output_top_zip_path: solvate_top.zip
    output_gro_path: solvate.gro
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    
step3M_grompp_genion:
  tool: grompp
  paths:
    input_gro_path: dependency/step3L_solvate/output_gro_path
    input_top_zip_path: dependency/step3L_solvate/output_top_zip_path
    output_tpr_path: gppion.tpr
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    simulation_type: minimization
    maxwarn: 10                    # NOTE: this will be ligand dependent!! :o - the warning is for each atom with redefined parameters

step3N_genion:
  tool: genion
  paths:
    input_tpr_path: dependency/step3M_grompp_genion/output_tpr_path
    input_top_zip_path: dependency/step3L_solvate/output_top_zip_path
    output_top_zip_path: genion_top.zip
    output_gro_path: genion.gro
  properties:
    binary_path: {gmx_binary_path}    # GROMACS binary path
    neutral: True                 # Neutralize charge of the system
    concentration: {salt_conc}    # Concentration of ions in mols/L

step3O_gro2pdb:
  tool: gmx_trjconv_str
  paths: 
    input_top_path: dependency/step3N_genion/output_gro_path
    input_structure_path: dependency/step3N_genion/output_gro_path
    output_str_path: topology.pdb
  properties:
    binary_path: {gmx_binary_path}   # GROMACS binary path

#############################################################################
# Section 4 (Steps A-I): Minimize and equilibrate the initial configuration #
#############################################################################

# NOTE: Leverage hierarchy in mdp when creating new eq steps without restraints - you can overwrite the default mdp config from the type of simulation

step4A_grompp_min:
  tool: grompp
  paths:
    input_gro_path: dependency/step3N_genion/output_gro_path
    input_top_zip_path: dependency/step3N_genion/output_top_zip_path
    output_tpr_path: gppmin.tpr
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    simulation_type: minimization
    mdp:
      integrator: steep
      nsteps: 10000
      emtol: 500
      emstep: 0.01
    
step4B_mdrun_min:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step4A_grompp_min/output_tpr_path
    output_trr_path: min.trr
    output_gro_path: min.gro
    output_edr_path: min.edr
    output_log_path: min.log
  properties:
    binary_path: {gmx_binary_path} # GROMACS binary path
    mpi_bin: {mpi_bin}             # MPI binary path, e.g. mpirun, srun... (should be null if gromacs already includes MPI support)                                               
    {mpi_np_config}                # Number of processors for MPI selected in the MPI call
    {num_threads_mpi_config}
    {num_threads_omp_config}
    
step4C_make_ndx:
  tool: make_ndx 
  paths:
    input_structure_path: dependency/step4B_mdrun_min/output_gro_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: '"System"'

step4D_energy_min:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step4B_mdrun_min/output_edr_path
    output_xvg_path: min_ene.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    terms: ["Potential"]
    xvg: xmgr 

step4E_grompp_nvt:
  tool: grompp
  paths:
    input_gro_path: dependency/step4B_mdrun_min/output_gro_path
    input_ndx_path: dependency/step4C_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step3N_genion/output_top_zip_path
    output_tpr_path: gppnvt.tpr
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    simulation_type: nvt
    mdp:
      ref-t: {temp} {temp}
      tc-grps: "Protein Water_and_ions"
      nsteps: {equil_nsteps} 
      dt: {dt}
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {equil_traj_freq}
      
step4F_mdrun_nvt:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step4E_grompp_nvt/output_tpr_path
    output_xtc_path: nvt.xtc
    output_trr_path: nvt.trr
    output_gro_path: nvt.gro
    output_edr_path: nvt.edr
    output_log_path: nvt.log
    output_cpt_path: nvt.cpt
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    mpi_bin: {mpi_bin}             # MPI binary path, e.g. mpirun, srun... (should be null if gromacs already includes MPI support)      
    {mpi_np_config}                # Number of processors for MPI selected in the MPI call
    use_gpu: {use_gpu}             # Wether to use GPU support or not
    {num_threads_mpi_config}
    {num_threads_omp_config}
    
step4G_temp_nvt:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step4F_mdrun_nvt/output_edr_path
    output_xvg_path: nvt_temp.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    terms: ["Temperature"]
    xvg: xmgr 

step4H_grompp_npt:
  tool: grompp
  paths:
    input_gro_path: dependency/step4F_mdrun_nvt/output_gro_path
    input_ndx_path: dependency/step4C_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step3N_genion/output_top_zip_path
    input_cpt_path: dependency/step4F_mdrun_nvt/output_cpt_path
    output_tpr_path: gppnpt.tpr
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    simulation_type: npt
    mdp:
      pcoupltype: isotropic
      nsteps: {equil_nsteps} 
      dt: {dt}
      ref-t: {temp} {temp}
      tc-grps: "Protein Water_and_ions"
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {equil_traj_freq}

step4I_mdrun_npt:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step4H_grompp_npt/output_tpr_path
    output_xtc_path: npt.xtc
    output_trr_path: npt.trr
    output_gro_path: npt.gro
    output_edr_path: npt.edr
    output_log_path: npt.log
    output_cpt_path: npt.cpt
  properties: 
    binary_path: {gmx_binary_path}     # GROMACS binary path
    mpi_bin: {mpi_bin}             # MPI binary path, e.g. mpirun, srun... (should be null if gromacs already includes MPI support)      
    {mpi_np_config}                # Number of processors for MPI selected in the MPI call
    use_gpu: {use_gpu}             # Wether to use GPU support or not
    {num_threads_mpi_config}
    {num_threads_omp_config}

step4J_density_npt:
  tool: gmx_energy
  paths:
    input_energy_path: dependency/step4I_mdrun_npt/output_edr_path
    output_xvg_path: npt_press_den.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    terms: ["Pressure", "Density"]
    xvg: xmgr 

############################################
# Section 5 (Steps A-B): MD production run #
############################################

step5A_grompp_md:
  tool: grompp
  paths:
    input_gro_path: dependency/step4I_mdrun_npt/output_gro_path
    input_cpt_path: dependency/step4I_mdrun_npt/output_cpt_path 
    input_ndx_path: dependency/step4C_make_ndx/output_ndx_path
    input_top_zip_path: dependency/step3N_genion/output_top_zip_path
    output_tpr_path: gppmd.tpr
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    simulation_type: free
    mdp:
      nsteps: {nsteps}
      dt: {dt} 
      ref-t: {temp} {temp}
      tc-grps: "Protein Water_and_ions"
      nstxout: 0           
      nstvout: 0
      nstxout-compressed: {traj_freq}
      nstenergy: 500
      continuation: 'yes'
      gen-vel: 'no'          
      ld-seed: 1

step5B_mdrun_md:
  tool: mdrun
  paths:
    input_tpr_path: dependency/step5A_grompp_md/output_tpr_path
    output_xtc_path: md.xtc
    output_trr_path: md.trr
    output_gro_path: md.gro
    output_edr_path: md.edr
    output_log_path: md.log
    output_cpt_path: md.cpt
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    mpi_bin: {mpi_bin}             # MPI binary path, e.g. mpirun, srun... (should be null if gromacs already includes MPI support)      
    {mpi_np_config}                # Number of processors for MPI selected in the MPI call
    use_gpu: {use_gpu}             # Wether to use GPU support or not
    {num_threads_mpi_config}
    {num_threads_omp_config}
        
#########################################
# Section 6 (Steps A-D): Basic analysis #
#########################################

step6A_rmsd_equilibrated:
  tool: gmx_rms
  paths:
    input_structure_path: dependency/step4I_mdrun_npt/output_gro_path
    input_traj_path: dependency/step5B_mdrun_md/output_xtc_path
    output_xvg_path: md_rmsdfirst.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: Backbone
    xvg: xmgr 

step6B_rmsd_experimental:
  tool: gmx_rms
  paths:
    input_structure_path: dependency/step4A_grompp_min/input_gro_path
    input_traj_path: dependency/step5B_mdrun_md/output_xtc_path
    output_xvg_path: md_rmsdexp.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: Backbone
    xvg: xmgr 
  
step6C_rgyr:
  tool: gmx_rgyr
  paths:
    input_structure_path: dependency/step4I_mdrun_npt/output_gro_path
    input_traj_path: dependency/step5B_mdrun_md/output_xtc_path
    output_xvg_path: md_rgyr.xvg
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: Backbone
    xvg: xmgr 

step6D_rmsf:
  tool: cpptraj_rmsf
  paths:
    input_top_path: dependency/step3O_gro2pdb/output_str_path
    input_traj_path: dependency/step5B_mdrun_md/output_xtc_path
    output_cpptraj_path: md_rmsf.xmgr   # .dat, .agr, .xmgr, .gnu
  properties:
    start: 1
    end: -1
    steps: 1
    mask: "!@H=" # by default cpptraj already strips solvent atoms

#####################################################
# Section 7 (Steps A-D): Trajectory post-processing #
#####################################################

# Optional step: used only if trajectories are different parts
step7A_trjcat:
  tool: trjcat
  paths:
    input_trj_zip_path: all_trajectories.zip
    output_trj_path: all_trajectories.xtc
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    concatenate: True

step7B_dry_str:
  tool: gmx_trjconv_str
  paths: 
    input_structure_path: dependency/step4I_mdrun_npt/output_gro_path
    input_top_path: dependency/step5A_grompp_md/output_tpr_path
    input_index_path: dependency/step4C_make_ndx/output_ndx_path
    output_str_path: dry_structure.gro
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: Protein
    center: True
    pbc: mol
    ur: compact

step7C_dry_traj:
  tool: gmx_trjconv_trj
  paths: 
    input_traj_path: dependency/step5B_mdrun_md/output_xtc_path
    input_top_path: dependency/step5A_grompp_md/output_tpr_path
    input_index_path: dependency/step4C_make_ndx/output_ndx_path
    output_traj_path: dry_traj.xtc
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    selection: Protein

step7D_center:
  tool: gmx_image 
  paths:
    input_traj_path: dependency/step7C_dry_traj/output_traj_path
    input_top_path: dependency/step5A_grompp_md/output_tpr_path
    input_index_path: dependency/step4C_make_ndx/output_ndx_path
    output_traj_path: center_traj.xtc
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    center_selection: Protein
    output_selection: Protein
    center: True
    ur: compact
    pbc: none

step7E_image_traj:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7D_center/output_traj_path
    input_top_path: dependency/step5A_grompp_md/output_tpr_path
    input_index_path: dependency/step4C_make_ndx/output_ndx_path
    output_traj_path: imaged_traj.xtc
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
    output_selection: Protein
    cluster_selection: Protein
    center_selection: Protein  # NOTE: why is this used??
    center: False
    ur: compact
    pbc: mol

step7F_fit_traj:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7E_image_traj/output_traj_path
    input_top_path: dependency/step5A_grompp_md/output_tpr_path
    input_index_path: dependency/step4C_make_ndx/output_ndx_path
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_binary_path}     # GROMACS binary path
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
def main_wf(configuration_path: Optional[str] = None, 
            input_pdb_path: Optional[str] = None, 
            ligands_top_folder: Optional[str] = None, 
            gmx_binary_path: str = 'gmx',
            mpi_bin: Optional[str] = None,
            num_threads_mpi: Optional[int] = None,
            num_threads_omp: Optional[int] = None,
            use_gpu: bool = False,
            restart: bool = False,
            forcefield: str = 'amber99sb-ildn', 
            salt_concentration: float = 0.15,
            temperature: float = 300.0,
            setup_only: bool = False, 
            input_gro_path: str = None, 
            input_top_path: str = None, 
            dt: float = 0.002,
            equil_nsteps: int = 500000,
            equil_traj_freq: Optional[int] = None,
            equil_only: bool = False, 
            nsteps: float = 50000000,
            traj_freq: Optional[int] = None,
            num_parts: int = 1, 
            num_replicas: int = 1, 
            skip_traj_processing: bool = False, 
            output_path: str = None
    ):
    '''
    Main MD setup and run workflow with GROMACS. Can be used to prepare and launch an MD simulation.

    Inputs
    ------

        configuration_path: 
            path to YAML configuration file
        input_pdb_path: 
            path to input PDB file
        ligands_top_folder:
            path to folder with ligands topology files (.itp) and coordinates (.gro)
        gmx_binary_path:
            path to GROMACS binary
        mpi_bin:
            path to MPI binary (e.g. mpirun, srun...) (should be null if gromacs already includes MPI support)
        num_threads_mpi:
            number of threads to be used by MPI
        num_threads_omp:
            number of threads to be used by OpenMP
        use_gpu:
            whether to use GPU support or not
        restart:
            whether to restart the workflow from the last completed step or start from the beginning.
        forcefield: 
            forcefield to be used in the simulation. Default: amber99sb-ildn. 
            See values supported by pdb2gmx (gromos45a3, charmm27, gromos53a6, amber96, amber99, 
            gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, amber99sb, amber03, amber99sb-ildn, 
            oplsaa, amber94, amber99sb-star-ildn-mut). 
        salt_concentration: 
            salt concentration to be used in the simulation. Default: 0.15.
        temperature:
            temperature to be used in the simulation. Default: 300.0.
        setup_only: 
            whether to only setup the system or also run the simulations
        input_gro_path: 
            path to already-prepared input structure file (.gro)
        input_top_path: 
            path to already-prepared input topology file (.zip)
        dt:
            time step to be used in the simulation. Default: 0.002 (2 fs).
        equil_nsteps:
            number of steps of the equilibration simulations. For NVT and for NPT. Default: 500000 (1 ns).
        equil_traj_freq:
            Saving frequency of the trajectory during equilibration in number of steps.
            Default: equil_nsteps // 500.
        equil_only: 
            whether to only run the equilibration or also run the production simulations
        nsteps: 
            Total number of steps of the production simulation. Default: 50000000 (100 ns).
        traj_freq:
            Saving frequency of the trajectory during production in number of steps.
            Default: nsteps // 1000.
        num_parts: 
            number of parts of the trajectory 
        num_replicas: 
            number of replicas of the trajectory
        skip_traj_processing:
            Skip the trajectory post-processing. Otherwise the trajectory will be dried, centered,
            imaged and fitted using 'gmx trjconv'. Default: False
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
        'gmx_binary_path': gmx_binary_path,
        'mpi_bin': mpi_bin,
        'num_threads_mpi': num_threads_mpi,
        'num_threads_omp': num_threads_omp,
        'use_gpu': use_gpu,
        'restart': restart,
        'forcefield': forcefield,
        'salt_conc': salt_concentration,
        'temp': temperature,
        'dt': dt,
        'equil_nsteps': equil_nsteps,
        'equil_traj_freq': equil_traj_freq,
        'nsteps': nsteps,
        'traj_freq': traj_freq
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
    
    # Check input files
    check_inputs(input_pdb_path, 
                 input_gro_path, 
                 input_top_path, 
                 ligands_top_folder, 
                 global_log)
    
    # If prepared structure is not provided
    if input_gro_path is None:
            
        # Get the chain residues
        # NOTE: we are assuming this info is not changed by pdb2gmx
        chains_dict = get_chains_dict(input_pdb_path) 
        
        ###########################################
        # Prepare topology and coordinates for MD #
        ###########################################
               
        # STEP 3B: add H atoms, generate coordinate (.gro) and topology (.top) file for the system
        # Histidine protonation state are determined from resname
        # Other protonation states are chosen interactively
        global_log.info("step3B_structure_topology: Generate the topology")
        global_paths["step3B_structure_topology"]["input_pdb_path"] = input_pdb_path
        global_prop["step3B_structure_topology"]["force_field"]=forcefield
        global_log.info(f"step3B_structure_topology: Reading protonation states for titratable residues (0: deprotonated, 1: protonated)")
        for residue in gmx_titra_resnames.keys():
            global_prop["step3B_structure_topology"][residue.lower()] = read_protonation_states(input_pdb_path, residue, global_log)
            global_log.info(f"step3B_structure_topology: {residue} protonation state: {global_prop['step3B_structure_topology'][residue.lower()]}")
        pdb2gmx(**global_paths["step3B_structure_topology"], properties=global_prop["step3B_structure_topology"])
        
        master_index_file = ""
        if len(chains_dict) > 1:
            
            # Add chain and chain CA groups to a master index file 
            for chain_id in chains_dict:
                
                prefix = f"step3B_Chain_{chain_id}"
                chain_prop = conf.get_prop_dic(prefix=prefix)
                chain_paths = conf.get_paths_dic(prefix=prefix)
                
                structure_path = global_paths["step3B_structure_topology"]["output_gro_path"]
                
                # If not the first chain, we need to append the group to the master index file
                if master_index_file:
                    chain_paths["step3C_make_ref_group"]["input_ndx_path"] = master_index_file
                
                # STEP 3B: Create index file for chain
                global_log.info(f"{chain_id} > Create index group for the chain")
                chain_paths["step3C_make_ref_group"]["input_structure_path"] = structure_path
                chain_prop["step3C_make_ref_group"]["selection"] = f"ri {chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
                make_ndx(**chain_paths["step3C_make_ref_group"], properties=chain_prop["step3C_make_ref_group"])
                
                global_log.info(f"{chain_id} > Create index group for the chain's C-alpha atoms")
                chain_paths["step3C_make_rest_group"]["input_structure_path"] = structure_path
                chain_prop["step3C_make_rest_group"]["selection"] = f"a CA & ri {chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
                make_ndx(**chain_paths["step3C_make_rest_group"], properties=chain_prop["step3C_make_rest_group"])  
                
                # Save group names for each chain
                chains_dict[chain_id]["reference_group"] = f"r_{chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
                chains_dict[chain_id]["restrain_group"] = f"CA_&_r_{chains_dict[chain_id]['residues'][0]}-{chains_dict[chain_id]['residues'][1]}"
                
                # Save POSRES names for each chain
                chains_dict[chain_id]["posres_name"] = f"CHAIN_{chain_id}_POSRES"
                
                # Update master index file
                master_index_file = chain_paths["step3C_make_rest_group"]["output_ndx_path"]

            # Reference, restraint and chain triplet list for ndx2resttop
            ref_rest_chain_triplet_list = ", ".join([f"({chains_dict[chain]['reference_group']}, {chains_dict[chain]['restrain_group']}, {chain})" for chain in chains_dict])
            
            # POSRES names for each chain
            posres_names = " ".join([chains_dict[chain]["posres_name"] for chain in chains_dict])
        else:
            
            chain_id = list(chains_dict.keys())[0]
            
            # Add chain and chain CA groups to a master index file
            global_log.info("step3C_make_ref_group: Create index file")
            make_ndx(**global_paths["step3C_make_ref_group"], properties=global_prop["step3C_make_ref_group"])
            
            # Save POSRES names for each chain
            chains_dict[chain_id]["posres_name"] = f"CHAIN_{chain_id}_POSRES"
             
            # Update master index file
            master_index_file = global_paths["step3C_make_ref_group"]["output_ndx_path"]
            
            # Reference, restraint and chain triplet list for ndx2resttop
            ref_rest_chain_triplet_list = f"(Protein, C-alpha, {chain_id})"
            
            # POSRES names for each chain
            posres_names = chains_dict[chain_id]["posres_name"]
        
        global_log.info(f"step3D_append_posres: Append restraints to the topology")
        global_prop["step3D_append_posres"]["ref_rest_chain_triplet_list"] = ref_rest_chain_triplet_list
        global_prop["step3D_append_posres"]["posres_names"] = posres_names
        global_paths["step3D_append_posres"]["input_ndx_path"] = master_index_file
        ndx2resttop(**global_paths["step3D_append_posres"], properties=global_prop["step3D_append_posres"])
        
        ligands_dict = get_ligands(ligands_top_folder, global_log)
        
        if ligands_dict:
            
            # STEP 3B: Convert gro of main structure to pdb - to concatenate with ligands
            global_log.info("step3E_structure_pdb: Convert GRO to PDB")
            gmx_trjconv_str(**global_paths["step3E_structure_pdb"], properties=global_prop["step3E_structure_pdb"])
            
            complex_pdb_path = global_paths["step3E_structure_pdb"]["output_str_path"]
            complex_topology_path = global_paths["step3D_append_posres"]["output_top_zip_path"]
            
            # For each ligand in the ligands folder
            for ligand_id in ligands_dict:
                
                prefix = f"step3FJ_{ligand_id}"
                ligand_prop = conf.get_prop_dic(prefix=prefix)
                ligand_paths = conf.get_paths_dic(prefix=prefix)
                
                ligand_gro_path = ligands_dict[ligand_id]["coordinates"]
                ligand_itp_path = ligands_dict[ligand_id]["topology"]
                ligand_restraints_path = os.path.join(str(Path(ligand_paths["step3I_ligand_restraints"]["output_itp_path"]).parent), f"{ligand_id}_posre.itp")
                
                ligands_dict[ligand_id]["posres_name"] = f"LIGAND_{ligand_id}_POSRES"
            
                # STEP 3C: Convert ligand coordinates from gro to pdb
                global_log.info(f"{ligand_id} > Convert ligand GRO to PDB")
                ligand_paths["step3F_ligand_pdb"]["input_structure_path"] = ligand_gro_path
                ligand_paths["step3F_ligand_pdb"]["input_top_path"] = ligand_gro_path
                gmx_trjconv_str(**ligand_paths["step3F_ligand_pdb"], properties=ligand_prop["step3F_ligand_pdb"])
                
                # STEP 3D: Create complex pdb file concatenating the current complex and the new ligand
                global_log.info(f"{ligand_id} > Create complex PDB file")
                ligand_paths["step3G_complex_pdb"]["input_structure1"] = complex_pdb_path
                ligand_paths["step3G_complex_pdb"]["output_structure_path"] = os.path.join(str(Path(ligand_paths["step3G_complex_pdb"]["output_structure_path"]).parent), f"{ligand_id}_complex.pdb")
                cat_pdb(**ligand_paths["step3G_complex_pdb"], properties=ligand_prop["step3G_complex_pdb"])
                
                # Update complex pdb path for the next ligands
                complex_pdb_path = ligand_paths["step3G_complex_pdb"]["output_structure_path"]
                
                # STEP 3E: Make ndx file for the ligand's heavy atoms
                global_log.info(f"{ligand_id} > Create index file for the ligand's heavy atoms")
                ligand_paths["step3H_make_ligand_ndx"]["input_structure_path"] = ligand_gro_path
                make_ndx(**ligand_paths["step3H_make_ligand_ndx"], properties=ligand_prop["step3H_make_ligand_ndx"])
                
                # STEP 3F: Generate restraints for the ligand's heavy atoms
                global_log.info(f"{ligand_id} > Generate restraints for ligand")
                ligand_paths["step3I_ligand_restraints"]["input_structure_path"] = ligand_gro_path
                ligand_paths["step3I_ligand_restraints"]["output_itp_path"] = ligand_restraints_path
                genrestr(**ligand_paths["step3I_ligand_restraints"], properties=ligand_prop["step3I_ligand_restraints"])
                
                # STEP 3G: Append parameterized ligand to the current complex topology zip file
                ligand_paths["step3J_append_ligand"]["input_top_zip_path"] = complex_topology_path
                ligand_paths["step3J_append_ligand"]["input_itp_path"] = ligand_itp_path
                ligand_paths["step3J_append_ligand"]["input_posres_itp_path"] = ligand_restraints_path
                ligand_prop["step3J_append_ligand"]["posres_name"] = ligands_dict[ligand_id]["posres_name"]
                global_log.info(f"{ligand_id} > Append ligand to the topology")
                append_ligand(**ligand_paths["step3J_append_ligand"], properties=ligand_prop["step3J_append_ligand"])
                
                # Update complex topology path for the next ligands
                complex_topology_path = ligand_paths["step3J_append_ligand"]["output_top_zip_path"]
                
            # Modify paths for the next steps
            global_paths["step3K_editconf"]["input_gro_path"] = complex_pdb_path
            global_paths["step3L_solvate"]["input_top_zip_path"] = complex_topology_path
            
        # STEP 3H: Create simulation box
        global_log.info("step3K_editconf: Create the solvent box")
        editconf(**global_paths["step3K_editconf"], properties=global_prop["step3K_editconf"])
        
        # STEP 3I: Add solvent molecules
        global_log.info("step3L_solvate: Fill the solvent box with water molecules")
        solvate(**global_paths["step3L_solvate"], properties=global_prop["step3L_solvate"])

        # STEP 3J: ion generation pre-processing
        global_log.info("step3M_grompp_genion: Preprocess ion generation")
        grompp(**global_paths["step3M_grompp_genion"], properties=global_prop["step3M_grompp_genion"])

        # STEP 3K: ion generation
        global_log.info("step3N_genion: Ion generation")
        global_prop["step3N_genion"]["concentration"] = salt_concentration
        genion(**global_paths["step3N_genion"], properties=global_prop["step3N_genion"])
        
        # Step 3L: conversion of topology from gro to pdb
        global_log.info("step3O_gro2pdb: Convert topology from GRO to PDB")
        gmx_trjconv_str(**global_paths["step3O_gro2pdb"], properties=global_prop["step3O_gro2pdb"])

        if setup_only:
            global_log.info("Set up only: setup_only flag is set to True! Exiting...")
            return

    else:
        
        global_log.info("Using prepared structure for MD")
        global_log.info(f"Input GRO file: {input_gro_path}")
        global_log.info(f"Input TOP file: {input_top_path}")
        
        # If prepared structure is provided, update the global paths
        global_paths['step4A_grompp_min']['input_gro_path'] = input_gro_path
        global_paths['step4A_grompp_min']['input_top_zip_path'] = input_top_path
        global_paths['step4E_grompp_nvt']['input_top_zip_path'] = input_top_path
        global_paths['step4H_grompp_npt']['input_top_zip_path'] = input_top_path
        global_paths['step5A_grompp_md']['input_top_zip_path'] = input_top_path

    # STEP 4A: minimization pre-processing
    global_log.info("step4A_grompp_min: Preprocess energy minimization")
    grompp(**global_paths["step4A_grompp_min"], properties=global_prop["step4A_grompp_min"])

    # STEP 4B: minimization
    global_log.info("step4B_mdrun_min: Execute energy minimization")
    mdrun(**global_paths["step4B_mdrun_min"], properties=global_prop["step4B_mdrun_min"])

    # STEP 4C: create index file
    if ligands_dict:
        global_prop["step4C_make_ndx"]["selection"] = f'"Protein" | r {" | r ".join(list(ligands_dict.keys()))}'
    global_log.info("step4C_make_ndx: Create index file")
    make_ndx(**global_paths["step4C_make_ndx"], properties=global_prop["step4C_make_ndx"])

    # STEP 4D: dump potential energy evolution
    global_log.info("step4D_energy_min: Compute potential energy during minimization")
    gmx_energy(**global_paths["step4D_energy_min"], properties=global_prop["step4D_energy_min"])

    chain_posres = " ".join([f"-D{chains_dict[chain]['posres_name']}" for chain in chains_dict])
    if ligands_dict:
        ligand_posres = " ".join([f"-D{ligands_dict[ligand]['posres_name']}" for ligand in ligands_dict])
    else:
        ligand_posres = ""
    eq_posres = f"{chain_posres} {ligand_posres}"
    
    # STEP 4E: NVT equilibration pre-processing
    if ligands_dict:
        global_prop["step4E_grompp_nvt"]["mdp"]["tc-grps"] = f"Protein_{'_'.join(list(ligands_dict.keys()))} Water_and_ions"
    global_prop["step4E_grompp_nvt"]["mdp"]["define"] = eq_posres
    global_prop["step4E_grompp_nvt"]["mdp"]["ref-t"] = f"{temperature} {temperature}"
    global_log.info("step4E_grompp_nvt: Preprocess system temperature equilibration")
    grompp(**global_paths["step4E_grompp_nvt"], properties=global_prop["step4E_grompp_nvt"])

    # STEP 4F: NVT equilibration
    global_log.info("step4F_mdrun_nvt: Execute system temperature equilibration")
    mdrun(**global_paths["step4F_mdrun_nvt"], properties=global_prop["step4F_mdrun_nvt"])

    # STEP 4G: dump temperature evolution
    global_log.info("step4G_temp_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step4G_temp_nvt"], properties=global_prop["step4G_temp_nvt"])

    # STEP 4H: NPT equilibration pre-processing
    if ligands_dict:
        global_prop["step4H_grompp_npt"]["mdp"]["tc-grps"] = f"Protein_{'_'.join(list(ligands_dict.keys()))} Water_and_ions"
    global_prop["step4H_grompp_npt"]["mdp"]["define"] = eq_posres
    global_prop["step4H_grompp_npt"]["mdp"]["ref-t"] = f"{temperature} {temperature}"
    global_log.info("step4H_grompp_npt: Preprocess system pressure equilibration")
    grompp(**global_paths["step4H_grompp_npt"], properties=global_prop["step4H_grompp_npt"])

    # STEP 4I: NPT equilibration
    global_log.info("step4I_mdrun_npt: Execute system pressure equilibration")
    mdrun(**global_paths["step4I_mdrun_npt"], properties=global_prop["step4I_mdrun_npt"])

    # STEP 4J: dump density and pressure evolution
    global_log.info("step4J_density_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step4J_density_npt"], properties=global_prop["step4J_density_npt"])
    
    # NOTE: add free equilibration removing those restraints that are not needed in the production run - none by default

    if equil_only:
        global_log.info("Equilibration only: equil_only flag is set to True! Exiting...")
        return
    
    ##########################
    # Production simulations #
    ##########################
    
    if num_replicas:
        # Folder names for replicas
        simulation_folders = [f"replica_{i}" for i in range(num_replicas)]
        global_log.info(f"Number of replicas: {num_replicas}")
    elif num_parts:
        # Folder names for parts
        simulation_folders = [f"parts_{i}" for i in range(num_parts)]
        global_log.info(f"Number of parts: {num_parts}")

    # Run each simulation (replica or part)
    traj_list = []
    for simulation in simulation_folders:
        
        traj_prop = conf.get_prop_dic(prefix=simulation)
        traj_paths = conf.get_paths_dic(prefix=simulation)

        # Update previous global paths needed by simulation-specific steps
        traj_paths['step5A_grompp_md']['input_gro_path'] = global_paths["step5A_grompp_md"]['input_gro_path']
        traj_paths['step5A_grompp_md']['input_cpt_path'] = global_paths["step5A_grompp_md"]['input_cpt_path']
        traj_paths['step5A_grompp_md']['input_top_zip_path'] = global_paths["step5A_grompp_md"]['input_top_zip_path']
        traj_paths['step5A_grompp_md']['input_ndx_path'] = global_paths["step5A_grompp_md"]['input_ndx_path']
        traj_paths['step6A_rmsd_equilibrated']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
        traj_paths['step6B_rmsd_experimental']['input_structure_path'] = global_paths["step4A_grompp_min"]['input_gro_path']
        traj_paths['step6C_rgyr']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
        traj_paths['step6D_rmsf']['input_top_path'] = global_paths["step3O_gro2pdb"]['output_str_path']
        
        # NOTE: num_replicas and num_parts will always have a value - review this
        
        # Simulations are replicas
        if num_replicas:
            # Change seed and velocities for each replica
            traj_prop['step5A_grompp_md']['mdp']['ld-seed'] = random.randint(1, 1000000)
            traj_prop['step5A_grompp_md']['mdp']['continuation'] = 'no'
            traj_prop['step5A_grompp_md']['mdp']['gen-vel'] = 'yes'

        # Simulations are parts of a single trajectory 
        if num_parts:
            # Divide the number of steps by the number of parts
            traj_prop['step5A_grompp_md']['mdp']['nsteps']=int(traj_prop['step5A_grompp_md']['mdp']['nsteps']/num_parts)
            # For all parts except the first one, use the previous gro and cpt files
            if simulation != simulation_folders[0]:
                traj_paths['step5A_grompp_md']['input_gro_path'] = previous_gro_path
                traj_paths['step5A_grompp_md']['input_cpt_path'] = previous_cpt_path

        # STEP 17: free NPT production run pre-processing
        if ligands_dict:
            traj_prop["step5A_grompp_md"]["mdp"]["tc-grps"] = f"Protein_{'_'.join(list(ligands_dict.keys()))} Water_and_ions"
        traj_prop["step5A_grompp_md"]["mdp"]["define"] = "" # NOTE: here restraint what is asked by the user
        traj_prop["step5A_grompp_md"]["mdp"]["ref-t"] = f"{temperature} {temperature}"
        global_log.info(f"{simulation} >  step5A_grompp_md: Preprocess free dynamics")
        grompp(**traj_paths['step5A_grompp_md'], properties=traj_prop["step5A_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{simulation} >  step5B_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step5B_mdrun_md'], properties=traj_prop['step5B_mdrun_md'])
        
        # Append the trajectory to the list
        traj_list.append(traj_paths['step5B_mdrun_md']['output_xtc_path'])
        
        # Update the previous gro and cpt files
        previous_gro_path = traj_paths['step5B_mdrun_md']['output_gro_path']
        previous_cpt_path = traj_paths['step5B_mdrun_md']['output_cpt_path']
        
        # STEP 19: compute the RMSD with respect to equilibrated structure
        global_log.info("step6A_rmsd_equilibrated: Compute Root Mean Square deviation against equilibrated structure")
        gmx_rms(**traj_paths['step6A_rmsd_equilibrated'], properties=traj_prop['step6A_rmsd_equilibrated'])
        
        # STEP 20: compute the RMSD with respect to minimized structure
        global_log.info("step6B_rmsd_experimental: Compute Root Mean Square deviation against minimized structure (exp)")
        gmx_rms(**traj_paths['step6B_rmsd_experimental'], properties=traj_prop['step6B_rmsd_experimental'])
        
        # STEP 21: compute the Radius of gyration
        global_log.info("step6C_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
        gmx_rgyr(**traj_paths['step6C_rgyr'], properties=traj_prop['step6C_rgyr'])

        # STEP 22: compute the RMSF
        global_log.info("step6D_rmsf: Compute Root Mean Square Fluctuation to measure the protein flexibility during the free MD simulation")
        cpptraj_rmsf(**traj_paths['step6D_rmsf'], properties=traj_prop['step6D_rmsf'])
        
    # Do the final analysis with all the previous parts or replicas
    if not skip_traj_processing:
        
        # If simulations are different parts of a single trajectory
        if num_parts:
            
            # Concatenate the analysis files that can be concatenated
            concatenate_gmx_analysis(conf, simulation_folders, output_path)
        
            # STEP 23: concatenate trajectories
            global_log.info("step7A_trjcat: Concatenate trajectories")
            fu.zip_list(zip_file=global_paths["step7A_trjcat"]['input_trj_zip_path'], file_list=traj_list)
            trjcat(**global_paths["step7A_trjcat"], properties=global_prop["step7A_trjcat"])
            
            # STEP 24: obtain dry structure
            global_log.info("step7B_dry_str: Obtain dry structure")
            gmx_trjconv_str(**global_paths["step7B_dry_str"], properties=global_prop["step7B_dry_str"])
    
            # STEP 25: obtain dry trajectory
            global_log.info("step7C_dry_traj: Obtain dry trajectory")
            traj_paths['step7C_dry_traj']['input_traj_path'] = traj_paths['step7A_trjcat']['output_trj_path']
            gmx_trjconv_trj(**global_paths["step7C_dry_traj"], properties=global_prop["step7C_dry_traj"])
        
            # Remove unused trajectory
            os.remove(global_paths["step7A_trjcat"]["output_trj_path"])
    
            # STEP 26: center the trajectory
            global_log.info(f"{simulation} > step7D_center: Center the trajectory")
            gmx_image(**traj_paths['step7D_center'], properties=traj_prop['step7D_center'])

            # Remove unused trajectory
            os.remove(global_paths["step7C_dry_traj"]["output_traj_path"])
            
            # STEP 27: image the trajectory
            global_log.info("step7E_image_traj: Imaging the trajectory")
            gmx_image(**global_paths['step7E_image_traj'], properties=global_prop['step7E_image_traj'])

            # STEP 28: fit the trajectory
            global_log.info("step7F_fit_traj: Fit the trajectory")
            gmx_image(**global_paths['step7F_fit_traj'], properties=global_prop['step7F_fit_traj'])
            
            # Remove unused trajectory
            os.remove(global_paths["step7E_image_traj"]["output_traj_path"])
            
        # If simulations are replicas
        if num_replicas:
            
            # For each replica, do the final analysis
            for simulation in simulation_folders:
                
                # Get the properties and paths for the replica
                traj_prop = conf.get_prop_dic(prefix=simulation)
                traj_paths = conf.get_paths_dic(prefix=simulation)

                # Update previous global paths needed by simulation-specific steps
                traj_paths['step7B_dry_str']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7B_dry_str']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7C_dry_traj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7D_center']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']  
                traj_paths['step7E_image_traj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7F_fit_traj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                
                # STEP 24: obtain dry structure
                global_log.info(f"{simulation} > step7B_dry_str: Obtain dry structure")
                gmx_trjconv_str(**traj_paths["step7B_dry_str"], properties=traj_prop["step7B_dry_str"])
                
                # STEP 25: obtain dry trajectory
                global_log.info(f"{simulation} > step7C_dry_traj: Obtain dry trajectory")
                gmx_trjconv_trj(**traj_paths["step7C_dry_traj"], properties=traj_prop["step7C_dry_traj"])
                
                # STEP 26: center the trajectory
                global_log.info(f"{simulation} > step7D_center: Center the trajectory")
                gmx_image(**traj_paths['step7D_center'], properties=traj_prop['step7D_center'])
                
                # Remove unused trajectory
                os.remove(traj_paths["step7C_dry_traj"]["output_traj_path"])
                
                # STEP 27: image the trajectory
                global_log.info(f"{simulation} > step7E_image_traj: Imaging the trajectory")
                gmx_image(**traj_paths['step7E_image_traj'], properties=traj_prop['step7E_image_traj'])
                
                # STEP 28: fit the trajectory
                global_log.info(f"{simulation} > step7F_fit_traj: Fit the trajectory")
                gmx_image(**traj_paths['step7F_fit_traj'], properties=traj_prop['step7F_fit_traj'])
                
                # Remove unused trajectory
                os.remove(traj_paths["step7E_image_traj"]["output_traj_path"])

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

    parser.add_argument('--config', dest='config_path', type=str,
                        help="Configuration file (YAML)",
                        required=False)

    parser.add_argument('--input_pdb', dest='input_pdb_path', type=str,
                        help="""Input PDB file. The workflow assumes the protonation state specified by the residue names
                        is the correct one. Default: None""",
                        required=False)

    parser.add_argument('--ligands_folder', dest='ligands_top_folder', type=str,
                        help="""Path to folder with .itp and .gro files for the ligands that 
                        should be included in the simulation. Make sure the coordinates of the 
                        ligands correspond to the PDB used. Default: None""",
                        required=False)

    parser.add_argument('--gmx_bin', dest='gmx_binary_path', type=str,
                        help="Path to GROMACS binary. Default: gmx",
                        required=False, default='gmx')
    
    parser.add_argument('--mpi_bin', dest='mpi_bin', type=str,
                        help="Path to MPI binary. Default: null",
                        required=False, default='null')
    
    parser.add_argument('--num_threads_mpi', dest='num_threads_mpi', type=int,
                        help="Number of MPI threads. Default: Let GROMACS guess",
                        required=False, default=1)
    
    parser.add_argument('--num_threads_omp', dest='num_threads_omp', type=int,
                        help="Number of OpenMP threads. Default: Let GROMACS guess",
                        required=False) 
    
    parser.add_argument('--use_gpu', action='store_true',
                        help="Use GPU for GROMACS. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--restart', action='store_true',
                        help="Restart the workflow from the last completed step. Default: False",
                        required=False, default=False)

    parser.add_argument('--forcefield', dest='forcefield', type=str,
                        help="Forcefield to use. Default: amber99sb-ildn",
                        required=False, default='amber99sb-ildn')

    parser.add_argument('--salt_conc', dest='salt_concentration', type=float,
                        help="Concentration of ions in the system. Default: 0.15",
                        required=False, default=0.15)
    
    parser.add_argument('--temp', dest='temperature', type=float,
                        help="Temperature of the system in K. Default: 300",
                        required=False, default=300)
    
    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system. Default: False",
                        required=False, default=False)

    parser.add_argument('--input_gro', dest='input_gro_path', type=str,
                        help="""Input structure file ready to minimize (.gro). To provide an externally prepared system, 
                        use together with --input_top (default: None)""",
                        required=False)

    parser.add_argument('--input_top', dest='input_top_path', type=str,
                        help="""Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, 
                        use together with --input_gro (default: None)""",
                        required=False)
    # NOTE: using input gro and top we don't have access to the pdb and thus we don't know which POSRES to apply - 
    # chains_dict and ligands_dict are not created
    
    parser.add_argument('--dt', dest='dt', type=float,
                        help="Time step in ps. Default: 0.002",
                        required=False, default=0.002)
    
    parser.add_argument('--equil_nsteps', dest='equil_nsteps', type=int,
                        help="Number of steps of the equilibration simulations. For NVT and for NPT. Default: 500000",
                        required=False, default=500000)
    
    parser.add_argument('--equil_traj_freq', dest='equil_traj_freq', type=int,
                        help="Saving frequency of the trajectory during equilibration in number of steps. Default: equil_nsteps // 500.",
                        required=False)
    
    parser.add_argument('--equil_only', action='store_true',
                        help="Only run the equilibration steps. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--nsteps', dest='nsteps', type=int,
                        help="Number of steps of the production simulation. Default: 50000000",
                        required=False, default=50000000)
    
    parser.add_argument('--traj_freq', dest='traj_freq', type=int,
                        help="Saving frequency of the trajectory during production in number of steps. Default: nsteps // 1000.",
                        required=False)
    
    parser.add_argument('--num_parts', dest='num_parts', type=int,
                        help="Number of parts to divide the simulation into. Default: 1",
                        required=False)
    
    parser.add_argument('--num_replicas', dest='num_replicas', type=int,
                        help="Number of replicas with different seeds to run the simulation. Default: 1",
                        required=False)

    parser.add_argument('--skip_traj_processing', action='store_true', dest='skip_traj_processing', 
                        help="""Skip the trajectory post-processing. Otherwise the trajectory will be dried, centered,
                        imaged and fitted using 'gmx trjconv'. Default: False""",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path', type=str,
                        help="Output path. Default: 'output' in the current working directory",
                        required=False, default='output')
    
    # NOTE: add flag to determine what should remain restrained during the production run - currently everything is free always

    args = parser.parse_args()

    # Check .pdb structure and .gro/.zip topology are not given together
    both_pdb_and_gro = args.input_pdb_path is not None and args.input_gro_path is not None
    if both_pdb_and_gro:
        raise Exception("Both --input_pdb and --input_gro are provided. Please provide only one of them")
    
    # Check .gro structure and .zip topology are given together
    only_topology = args.input_gro_path is None and args.input_top_path is not None
    only_coordinates = args.input_gro_path is not None and args.input_top_path is None
    if only_topology and only_coordinates:
        raise Exception("Both --input_gro and --input_top must be provided together")

    # Convert to corresponding types
    if args.salt_concentration:
        args.salt_concentration = float(args.salt_concentration)
    if args.temperature:
        args.temperature = float(args.temperature)
    if args.nsteps:
        args.nsteps = int(args.nsteps)
    if args.num_parts:
        args.num_parts = int(args.num_parts)
    if args.num_replicas:
        args.num_replicas = int(args.num_replicas)
    
    main_wf(configuration_path=args.config_path, 
            input_pdb_path=args.input_pdb_path, 
            ligands_top_folder=args.ligands_top_folder, 
            gmx_binary_path=args.gmx_binary_path,
            mpi_bin=args.mpi_bin,
            num_threads_mpi=args.num_threads_mpi,
            num_threads_omp=args.num_threads_omp,
            use_gpu=args.use_gpu,
            restart=args.restart,
            forcefield=args.forcefield, 
            salt_concentration=args.salt_concentration,
            temperature=args.temperature,
            setup_only=args.setup_only, 
            input_gro_path=args.input_gro_path, 
            input_top_path=args.input_top_path,
            dt=args.dt, 
            equil_nsteps=args.equil_nsteps,
            equil_traj_freq=args.equil_traj_freq,
            equil_only=args.equil_only, 
            nsteps=args.nsteps, 
            traj_freq=args.traj_freq, 
            num_parts=args.num_parts, 
            num_replicas=args.num_replicas, 
            skip_traj_processing=args.skip_traj_processing, 
            output_path=args.output_path)