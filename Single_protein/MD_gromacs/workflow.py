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
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf
from biobb_structure_utils.utils.cat_pdb import cat_pdb
from biobb_structure_utils.utils.extract_molecule import extract_molecule
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

def get_ligands(ligands_folder: Union[str, None], global_log) -> List[Dict[str, str]]:
    """
    Get a list of available ligands in the ligands folder. The function searches for all the .itp and .gro files in the ligands folder
    If the ligands folder is provided but doesn't exist or any of the ligands is missing the topology or coordinate file, an error is raised.
    
    Inputs
    ------
    
        ligands_folder (str): Path to the folder with the ligand .itp and .gro files.
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
    if ligands_folder is None:
        return ligands
    
    # Check if the ligands folder exists
    if not os.path.exists(ligands_folder):
        global_log.error(f"Folder {ligands_folder} not found")
        return ligands
    
    # Search for ligands in the ligands folder
    for file in os.listdir(ligands_folder):
        
        # Check if the file is a .itp or .gro file
        if file.endswith(".itp") or file.endswith(".gro"):
            
            # Get the file name without extension
            ligand_id = Path(file).stem   
            
            # Get the file extension
            file_extension = Path(file).suffix
            
            # Check if the ligand name is already in the dictionary
            if ligand_id in ligands:
                if file_extension == ".itp":
                    ligands[ligand_id]['topology'] = os.path.join(ligands_folder, file)
                elif file_extension == ".gro":
                    ligands[ligand_id]['coordinates'] = os.path.join(ligands_folder, file)
            else:
                if file_extension == ".itp":
                    ligands[ligand_id] = {'topology': os.path.join(ligands_folder, file)}
                elif file_extension == ".gro":
                    ligands[ligand_id] = {'coordinates': os.path.join(ligands_folder, file)}
    
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

def find_amber_his(pdb_file: str, global_log) -> List[str]:
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

def get_pdb2gmx_his(his_residues: List[str]) -> str:
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

# Set additional general properties not considered by the configuration reader
def set_gromacs_path(global_properties: dict, binary_path: str) -> None:
    """
    Set the path to the GROMACS binary for all steps using GROMACS.

    Inputs
    ------

        global_properties (dict): Dictionary containing the global_properties.
        binary_path (str): Path to the GROMACS binary.
    """

    list_of_steps = ['step3B_structure_topology', 'step3K_editconf', 'step3L_solvate', 'step3M_grompp_genion', 'step3N_genion',
                        'step4A_grompp_min', 'step4B_mdrun_min', 'step4E_grompp_nvt', 'step4F_mdrun_nvt',
                        'step4H_grompp_npt', 'step4I_mdrun_npt', 'step5A_grompp_md', 'step5B_mdrun_md']

    for step in list_of_steps:
        global_properties[step]['binary_path'] = binary_path

def set_mpi_path(global_properties: dict, mpi_bin: str, mpi_np: int) -> None:
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

def set_gpu_use(global_properties: dict, gpu_use: bool) -> None:
    """
    Set the use of GPU for all steps using GROMACS that support it.

    Inputs
    ------

        global_properties (dict): Dictionary containing the global_properties.
        gpu_use (bool): Whether to use GPU or not.
    """

    list_of_steps = ['step4F_mdrun_nvt', 'step4I_mdrun_npt', 'step5B_mdrun_md']

    for step in list_of_steps:
        global_properties[step]['use_gpu'] = gpu_use

def set_gmx_properties(global_properties: dict, gmx_properties: dict, global_log) -> None:
    """
    Set all the gmx global properties of this workflow, i.e. those global properties included at the beginning of the YAML configuration file that
    are general to some gmx steps.
    
    Inputs
    ------
    
        global_properties (dict): Dictionary containing the global_properties.
        gmx_properties (dict): Dictionary containing the gmx properties.
        global_log (Logger): Logger object for logging messages.
    """
    
    # Enforce gromacs binary path for all steps using gromacs
    if gmx_properties.get('binary_path'):
        global_log.info(f"Using GROMACS binary path: {gmx_properties['binary_path']}")
        set_gromacs_path(global_properties, gmx_properties['binary_path'])

    # Enforce mpi binary path for all steps using mpi
    if gmx_properties.get('mpi_bin'):
        global_log.info(f"Using MPI binary path: {gmx_properties['mpi_bin']}")
        set_mpi_path(global_properties, gmx_properties['mpi_bin'], gmx_properties.get('mpi_np'))

    # Enforce gpu use for all steps using gromacs that support it
    if gmx_properties.get('use_gpu'):
        global_log.info(f"Using GPU for GROMACS steps")
        set_gpu_use(global_properties, gmx_properties['use_gpu'])


# Process topology - temporal solution 
def process_ligand_top(input_path, output_path) -> None:
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
    filepaths : list of str
        A list of file paths for the input XVG files to be merged. 
        Each file should follow the format of time series data with comments starting 
        with '#' followed by a header section starting with '@' with plot details 
    
    output_filepath : str
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
    
def concatenate_gmx_analysis(conf, simulation_folders, output_path) -> None:
    """
    Concatenates the analysis files for each step of the GROMACS analysis including
    RMSD and Rgyr. The function reads the analysis files for each step from the
    simulation folders, merges them into a single file, and writes the merged data
    
    Inputs:
    -------
    
    conf : class settings.ConfReader
        Configuration file reader object
    simulation_folders : list of str
        List of folder names containing the simulation data for each part or replica
    output_path : str
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
    

def main_wf(configuration_path, input_pdb_path = None, pdb_code = None, pdb_chains = None, mutation_list = None, ligands_folder = None, skip_fix_backbone = None, skip_fix_side_chain = None, 
            fix_ss = None, fix_amide_clashes = None, his_protonation_tool = "pdb4amber", his = None, forcefield = 'amber99sb-ildn', setup_only = False, input_gro_path = None, input_top_path = None, 
            equil_only = False, nsteps = None, num_parts = 1, num_replicas = 1, final_analysis = None, output_path = None):
    '''
    Main setup, mutation and MD run workflow with GROMACS. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system, minimize it, equilibrate it and finally do N production runs (either replicas or parts).

    Inputs
    ------

        configuration_path   (str): path to YAML configuration file
        input_pdb_path       (str): (Optional) path to input PDB file
        pdb_code             (str): (Optional) PDB code to be used to get the canonical FASTA sequence
        pdb_chains           (str): (Optional) list of chains to be extracted from the PDB file and fixed
        mutation_list        (str): (Optional) list of mutations to be introduced in the structure
        ligands_folder       (str): (Optional) path to the folder containing the ligand .itp and .gro files
        skip_fix_backbone   (bool): (Optional) whether to skip the fix of the backbone atoms
        skip_fix_side_chain (bool): (Optional) whether to skip the fix of the side chain atoms
        fix_ss              (bool): (Optional) wether to add disulfide bonds
        fix_amide_clashes   (bool): (Optional) wether to flip clashing amides to relieve the clashes
        his_protonation_tool (str): (Optional) histidine protonation tool to be used (pdb4amber or pdb2gmx). Default: pdb4amber.
        his                  (str): (Optional) histidine protonation states to be used in the simulation. Default: None. See values supported by pdb2gmx
        forcefield           (str): (Optional) forcefield to be used in the simulation. Default: amber99sb-ildn. See values supported by pdb2gmx 
                                    (gromos45a3, charmm27, gromos53a6, amber96, amber99, gromos43a2, gromos54a7, gromos43a1, amberGS, gromos53a5, 
                                    amber99sb, amber03, amber99sb-ildn, oplsaa, amber94, amber99sb-star-ildn-mut). 
        setup_only          (bool): (Optional) whether to only setup the system or also run the simulations
        input_gro_path       (str): (Optional) path to already-prepared input structure file (.gro)
        input_top_path       (str): (Optional) path to already-prepared input topology file (.zip)
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
        
    # Remove gmx-specific properties from global properties - otherwise they will be applied to all steps
    gmx_properties = conf.global_properties.pop('gmx', None)

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
    set_gmx_properties(global_prop, gmx_properties, global_log)

    ##############################################
    # Extract atoms and prepare structure for MD #
    ##############################################
    
    # If prepared structure is not provided
    if input_gro_path is None:

        # If input PDB is given as argument
        if input_pdb_path is not None:
            global_paths["step1A_extractAtoms"]["input_structure_path"] = input_pdb_path

        # If chains are given as argument
        if pdb_chains is not None:
            global_prop["step1A_extractAtoms"]["molecule_type"] = "chains"
            global_prop["step1A_extractAtoms"]["chains"] = pdb_chains
        
        # STEP 1A: extract main structure of interest while removing water and ligands (heteroatoms)
        global_log.info("step1A_extractAtoms: extract chain of interest")
        extract_molecule(**global_paths["step1A_extractAtoms"], properties=global_prop["step1A_extractAtoms"])
        
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

        # STEP 2I: Renumber structure atoms and residues
        global_log.info("step2I_renumberstructure: renumber structure")
        renumber_structure(**global_paths["step2I_renumberstructure"], properties=global_prop["step2I_renumberstructure"])
        
        # Get the chain residues
        chains_dict = get_chains_dict(global_paths["step2I_renumberstructure"]["output_structure_path"])
        
        ###########################################
        # Prepare topology and coordinates for MD #
        ###########################################
        
        # NOTE: if we have a gap that we are not modeling (e.g. a missing loop), pdb2gmx will find terminal atoms OXT in non-terminal residues and will return an error
        # STEP 3A: add H atoms, generate coordinate (.gro) and topology (.top) file for the system
        if (his_protonation_tool == "pdb4amber") and (his is None):
            
            global_log.info("step3A_amber_reduce: Determine the protonation states of histidine residues")
            
            # Find the histidine protonation states with pdb4amber
            pdb4amber_run(**global_paths["step3A_amber_reduce"], properties=global_prop["step3A_amber_reduce"])
            
            # Read the histidine protonation states from the pdb4amber output residue names
            his_residues = find_amber_his(global_paths["step3A_amber_reduce"]["output_pdb_path"], global_log)
            
            # Convert from histidine names to pdb2gmx numbering convention
            global_prop["step3B_structure_topology"]["his"]=get_pdb2gmx_his(his_residues)
            
        elif his is not None:
            # Set user-defined histidine protonation states
            global_prop["step3B_structure_topology"]["his"]=his
            
        global_log.info("step3B_structure_topology: Generate the topology")
        global_prop["step3B_structure_topology"]["force_field"]=forcefield
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
        
        ligands_dict = get_ligands(ligands_folder, global_log)
        
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
        genion(**global_paths["step3N_genion"], properties=global_prop["step3N_genion"])
        
        # Step 3L: conversion of topology from gro to pdb
        global_log.info("step3O_gro2pdb: Convert topology from GRO to PDB")
        gmx_trjconv_str(**global_paths["step3O_gro2pdb"], properties=global_prop["step3O_gro2pdb"])

        if setup_only:
            global_log.info("Set up only: setup_only flag is set to True! Exiting...")
            return

    else:
        
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
        set_gmx_properties(traj_prop, gmx_properties, global_log)

        # Update previous global paths needed by simulation-specific steps
        traj_paths['step5A_grompp_md']['input_gro_path'] = global_paths["step5A_grompp_md"]['input_gro_path']
        traj_paths['step5A_grompp_md']['input_cpt_path'] = global_paths["step5A_grompp_md"]['input_cpt_path']
        traj_paths['step5A_grompp_md']['input_top_zip_path'] = global_paths["step5A_grompp_md"]['input_top_zip_path']
        traj_paths['step5A_grompp_md']['input_ndx_path'] = global_paths["step5A_grompp_md"]['input_ndx_path']
        traj_paths['step6A_rmsd_equilibrated']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
        traj_paths['step6B_rmsd_experimental']['input_structure_path'] = global_paths["step4A_grompp_min"]['input_gro_path']
        traj_paths['step6C_rgyr']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
        traj_paths['step6D_rmsf']['input_top_path'] = global_paths["step3O_gro2pdb"]['output_str_path']
        
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
                traj_paths['step5A_grompp_md']['input_gro_path'] = previous_gro_path
                traj_paths['step5A_grompp_md']['input_cpt_path'] = previous_cpt_path

        # STEP 17: free NPT production run pre-processing
        if ligands_dict:
            traj_prop["step5A_grompp_md"]["mdp"]["tc-grps"] = f"Protein_{'_'.join(list(ligands_dict.keys()))} Water_and_ions"
        traj_prop["step5A_grompp_md"]["mdp"]["define"] = "" # NOTE: here restraint what is asked by the user
        global_log.info(f"{simulation} >  step5A_grompp_md: Preprocess free dynamics")
        grompp(**traj_paths['step5A_grompp_md'], properties=traj_prop["step5A_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{simulation} >  step5B_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step5B_mdrun_md'], properties=traj_prop['step5B_mdrun_md'])
        
        # Append the trajectory to the list
        traj_list.append(traj_paths['step5B_mdrun_md']['output_trr_path'])
        
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
                set_gmx_properties(traj_prop, gmx_properties, global_log)
                
                # NOTE: we are hard-coding the kind of traj that we are using with these paths: output_trr_path
                # Update previous global paths needed by simulation-specific steps
                traj_paths['step7B_dry_trj']['input_traj_path'] = traj_paths['step5B_mdrun_md']['output_trr_path']
                traj_paths['step7B_dry_trj']['input_top_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7B_dry_trj']['input_index_path'] = global_paths["step4C_make_ndx"]['output_ndx_path']
                traj_paths['step7C_dry_str']['input_structure_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
                traj_paths['step7C_dry_str']['input_top_path'] = global_paths["step4I_mdrun_npt"]['output_gro_path']
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

    parser = argparse.ArgumentParser("MD Simulation with GROMACS")

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

    parser.add_argument('--ligands_folder', dest='ligands_folder',
                        help="Folder with .itp and .gro files for the ligands that should be included in the simulation. Default: None",
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

    parser.add_argument('--his_protonation_tool', dest='his_protonation_tool',
                        help="Tool to use for histidine protonation (pdb4amber or pdb2gmx). Default: pdb4amber",
                        required=False, default='pdb4amber')
    
    parser.add_argument('--his', dest='his',
                        help="Histidine protonation states with pdb2gmx convention (HID: 0, HIE: 1, HIP:2). Overrides his_protonation_tool. Default: None. Example: '0 1 1'",
                        required=False)
    
    parser.add_argument('--forcefield', dest='forcefield',
                        help="Forcefield to use. Default: amber99sb-ildn",
                        required=False, default='amber99sb-ildn')

    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system. Default: False",
                        required=False, default=False)

    parser.add_argument('--input_gro', dest='input_gro_path',
                        help="Input structure file ready to minimize (.gro). To provide an externally prepared system, use together with --input_top (default: None)",
                        required=False)

    parser.add_argument('--input_top', dest='input_top_path',
                        help="Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, use together with --input_gro (default: None)",
                        required=False)
    
    # NOTE: using input gro and top we don't have access to the pdb and thus we don't know which POSRES to apply - chains_dict and ligands_dict are not created
    
    parser.add_argument('--equil_only', action='store_true',
                        help="Only run the equilibration steps. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--nsteps', dest='nsteps',
                        help="Number of steps of the simulation",
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
    if (args.input_gro_path is None and args.input_top_path is not None) or (args.input_gro_path is not None and args.input_top_path is None):
        raise Exception("Both --input_gro and --input_top must be provided together")

    # Check .pdb structure and .gro/.zip topology are not given together
    if (args.input_pdb_path is not None and args.input_gro_path is not None):
        raise Exception("Both --input_pdb and --input_gro/--input_top are provided. Please provide only one of them")

    main_wf(configuration_path=args.config_path, input_pdb_path=args.input_pdb_path, pdb_code=args.pdb_code, pdb_chains=args.pdb_chains, mutation_list=args.mutation_list, 
            ligands_folder=args.ligands_folder, skip_fix_backbone=args.skip_fix_backbone, skip_fix_side_chain=args.skip_fix_side_chain, fix_ss=args.fix_ss, fix_amide_clashes=args.fix_amide_clashes, 
            his_protonation_tool=args.his_protonation_tool, his=args.his, forcefield=args.forcefield, setup_only=args.setup_only, input_gro_path=args.input_gro_path, 
            input_top_path=args.input_top_path, equil_only=args.equil_only, nsteps=args.nsteps,  num_parts=args.num_parts, num_replicas=args.num_replicas, final_analysis=args.final_analysis, 
            output_path=args.output_path)