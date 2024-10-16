#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
from Bio.SeqIO.PdbIO import PdbSeqresIterator
from Bio.PDB import PDBParser
from Bio import SeqIO
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
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_pdb_tools.pdb_tools.biobb_pdb_tofasta import biobb_pdb_tofasta

def highest_occupancy_altlocs(pdb_file, global_log) -> list:
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

def set_gromacs_path(properties: dict, binary_path: str) -> None:
    """
    Set the path to the GROMACS binary for all steps using GROMACS.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        binary_path (str): Path to the GROMACS binary.
    """

    list_of_steps = ['step3_pdb2gmx', 'step4_editconf', 'step5_solvate', 'step6_grompp_genion', 'step7_genion',
                        'step8_grompp_min', 'step9_mdrun_min', 'step11_grompp_nvt', 'step12_mdrun_nvt',
                        'step14_grompp_npt', 'step15_mdrun_npt', 'step17_grompp_md', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['binary_path'] = binary_path

def set_mpi_path(properties: dict, mpi_bin: str, mpi_np: int) -> None:
    """
    Set the path to the MPI binary for all steps using MPI.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        mpi_bin (str): Path to the MPI binary.
        mpi_np (int): Number of processors to be used.
    """

    list_of_steps = ['step9_mdrun_min', 'step12_mdrun_nvt', 'step15_mdrun_npt', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['mpi_bin'] = mpi_bin
        properties[step]['mpi_np'] = mpi_np

def set_gpu_use(properties: dict, gpu_use: bool) -> None:
    """
    Set the use of GPU for all steps using GROMACS that support it.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        gpu_use (bool): Whether to use GPU or not.
    """

    list_of_steps = ['step12_mdrun_nvt', 'step15_mdrun_npt', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['use_gpu'] = gpu_use

def set_general_properties(properties: dict, conf, global_log) -> None:
    """
    Set all the additional global properties of this workflow, i.e. those properties included at the beginning of the YAML configuration file that
    are general to all steps and are not included already when the global properties are parsed.
    
    Inputs
    ------
    
        properties (dict): Dictionary containing the properties.
        conf (class settings.ConfReader): Configuration file reader.
    """
    
    # Enforce gromacs binary path for all steps using gromacs
    if conf.properties.get('binary_path'):
        global_log.info(f"Using GROMACS binary path: {conf.properties['binary_path']}")
        set_gromacs_path(properties, conf.properties['binary_path'])

    # Enforce mpi binary path for all steps using mpi
    if conf.properties.get('mpi_bin'):
        global_log.info(f"Using MPI binary path: {conf.properties['mpi_bin']}")
        set_mpi_path(properties, conf.properties['mpi_bin'], conf.properties.get('mpi_np'))

    # Enforce gpu use for all steps using gromacs that support it
    if conf.properties.get('use_gpu'):
        global_log.info(f"Using GPU for GROMACS steps")
        set_gpu_use(properties, conf.properties['use_gpu'])

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
    gmx_analysis_steps = ['step19_rmsd_equilibrated', 'step20_rmsd_experimental', 'step21_rgyr']
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
    
def main_wf(configuration_path, setup_only, num_parts, num_replicas, output_path = None, input_pdb_path = None, pdb_chains = None,
            mutation_list = None, input_gro_path = None, input_top_path = None, fix_ss = None, fix_amide_clashes = None, 
            his = None, nsteps = None, final_analysis = None):
    '''
    Main setup, mutation and MD run workflow with GROMACS. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system, minimize it, equilibrate it and finally do N production runs (either replicas or parts).

    Inputs
    ------

        configuration_path (str): path to YAML configuration file
        setup_only        (bool): (Optional) whether to only setup the system or also run the simulations
        num_parts          (int): (Optional) number of parts of the trajectory 
        num_replicas       (int): (Optional) number of replicas of the trajectory
        output_path        (str): (Optional) path to output folder
        input_pdb_path     (str): (Optional) path to input PDB file
        pdb_chains         (str): (Optional) chains to be extracted from the input PDB file
        mutation_list      (str): (Optional) list of mutations to be introduced in the structure
        input_gro_path     (str): (Optional) path to input structure file (.gro)
        input_top_path     (str): (Optional) path to input topology file (.zip)
        fix_ss            (bool): (Optional) wether to add disulfide bonds
        fix_amide_clashes (bool): (Optional) wether to flip clashing amides to relieve the clashes
        his               (str): (Optional) histidine protonation states list
        nsteps            (int): (Optional) Total number of steps of the production simulation
        final_analysis    (bool): (Optional) whether to perform the final analysis or not
        
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
    
    # Check num parts and num of replicas
    if (num_replicas is None) and (num_parts is None):
        # Default behavior: 1 replica
        num_replicas=1
    elif (num_replicas is not None) and (num_parts is not None):
        # Cannot set both num parts and num replicas
        global_log.error("Number of trajectories and replicas cannot be set at the same time")

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Set general properties for all steps
    set_general_properties(global_prop, conf, global_log)

    # If prepared structure is not provided
    if input_gro_path is None:

        # If input PDB is given as argument
        if input_pdb_path is not None:
            global_paths["step1_extractMolecule"]["input_structure_path"] = input_pdb_path

        # If chains are given as argument
        if pdb_chains is not None:
            global_prop["step1_extractMolecule"]["chains"] = pdb_chains

        # STEP 1: extract molecules of interest
        global_log.info("step1_extractMolecule: extract molecule of interest (protein)")
        extract_molecule(**global_paths["step1_extractMolecule"], properties=global_prop["step1_extractMolecule"])

        # STEP 2 (A): Fix alternative locations
        global_log.info("step2A_fixaltlocs: Fix alternative locations")
        global_prop["step2A_fixaltlocs"]["altlocs"] = highest_occupancy_altlocs(global_paths["step1_extractMolecule"]["input_structure_path"], global_log)
        fix_altlocs(**global_paths["step2A_fixaltlocs"], properties=global_prop["step2A_fixaltlocs"])

        if mutation_list is not None:
            global_prop["step2B_mutations"]["mutation_list"] = ",".join(mutation_list)

        # STEP 2 (B): Add mutations if requested
        global_log.info("step2B_mutations: Preparing mutated structure")
        mutate(**global_paths["step2B_mutations"], properties=global_prop["step2B_mutations"])
        
        # STEP 2 (C): Get FASTA sequence to model the backbone
        try:
            # Try to get the canonical FASTA sequence with an http request from the PDB code
            global_log.info("step2C_canonical_fasta: Get canonical FASTA")
            canonical_fasta(**global_paths["step2C_canonical_fasta"], properties=global_prop["step2C_canonical_fasta"])
            fasta_available = True
        except:
            global_log.warning("step2C_canonical_fasta: Could not get canonical FASTA. Check the internet connection in the machine running the workflow. Trying to get the canonical FASTA from the PDB file...")
            fasta_available = False
            
        if not fasta_available:
            # Try to get the FASTA sequence from SEQRES records in the PDB file
            global_log.info("step2C_pdb_tofasta: Get FASTA from SEQRES of PDB file")
            fasta_available = fasta_from_pdb(global_paths["step1_extractMolecule"]["input_structure_path"], global_paths["step2C_pdb_tofasta"]["output_file_path"], global_log)

            # Update fix backbone input
            global_paths['step2D_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step2C_pdb_tofasta']['output_file_path']
            
        if not fasta_available:
            # Try to get the FASTA sequence from the PDB file
            global_log.info("step2C_pdb_tofasta: Get FASTA from PDB file")
            
            # Update the input file path
            global_paths['step2C_pdb_tofasta']['input_file_path'] = global_paths["step1_extractMolecule"]["input_structure_path"]
            
            # NOTE: this is not the canonical, only existing residues in the PDB file are included
            biobb_pdb_tofasta(**global_paths["step2C_pdb_tofasta"], properties=global_prop["step2C_pdb_tofasta"])
            
            # Update fix backbone input
            global_paths['step2D_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step2C_pdb_tofasta']['output_file_path']
            fasta_available = True
            
        if fasta_available:
            # STEP 2 (D): Model missing heavy atoms of backbone
            global_log.info("step2D_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
            fix_backbone(**global_paths["step2D_fixbackbone"], properties=global_prop["step2D_fixbackbone"])
        else:
            global_log.warning("step2D_fixbackbone: Could not get FASTA sequence. Skipping modeling of the missing heavy atoms in the backbone.")
            global_paths['step2E_fixsidechain']['input_pdb_path'] = global_paths['step2B_mutations']['output_pdb_path']

        # STEP 2 (E): model missing heavy atoms of side chains
        global_log.info("step2E_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        fix_side_chain(**global_paths["step2E_fixsidechain"], properties=global_prop["step2E_fixsidechain"])

        if fix_ss:
            # STEP 2 (F): model SS bonds (CYS -> CYX)
            global_log.info("step2F_fixssbonds: Fix SS bonds")
            fix_ssbonds(**global_paths["step2F_fixssbonds"], properties=global_prop["step2F_fixssbonds"])
        else:
            global_paths['step2G_fixamides']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

        if fix_amide_clashes:
            # STEP 2 (G): Fix amides
            global_log.info("step2G_fixamides: fix clashing amides")
            fix_amides(**global_paths["step2G_fixamides"], properties=global_prop["step2G_fixamides"])
        else:
            if fix_ss:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2F_fixssbonds']['output_pdb_path']
            else:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

        # STEP 2 (H): Fix chirality
        global_log.info("step2H_fixchirality: fix chirality of residues")
        fix_chirality(**global_paths["step2H_fixchirality"], properties=global_prop["step2H_fixchirality"])

        # STEP 2 (I): renumber structure atoms and residues
        global_log.info("step2I_renumberstructure: renumber structure")
        renumber_structure(**global_paths["step2I_renumberstructure"], properties=global_prop["step2I_renumberstructure"])
        
        # NOTE: Histidine protonation states come from external call to pdb4amber, should be done within the WF!
        # STEP 3: add H atoms, generate coordinate (.gro) and topology (.top) file
        global_log.info("step3_pdb2gmx: Generate the topology")
        if his:
            global_prop["step3_pdb2gmx"]["his"]=his
 
        pdb2gmx(**global_paths["step3_pdb2gmx"], properties=global_prop["step3_pdb2gmx"])

        # STEP 4: Create simulation box
        global_log.info("step4_editconf: Create the solvent box")
        editconf(**global_paths["step4_editconf"], properties=global_prop["step4_editconf"])

        # STEP 5: Add solvent molecules
        global_log.info("step5_solvate: Fill the solvent box with water molecules")
        solvate(**global_paths["step5_solvate"], properties=global_prop["step5_solvate"])

        # STEP 6: ion generation pre-processing
        global_log.info("step6_grompp_genion: Preprocess ion generation")
        grompp(**global_paths["step6_grompp_genion"], properties=global_prop["step6_grompp_genion"])

        # STEP 7: ion generation
        global_log.info("step7_genion: Ion generation")
        genion(**global_paths["step7_genion"], properties=global_prop["step7_genion"])
        
        # Step 7B: conversion of topology from gro to pdb
        global_log.info("step7B_gro2pdb: Convert topology from GRO to PDB")
        gmx_trjconv_str(**global_paths["step7B_gro2pdb"], properties=global_prop["step7B_gro2pdb"])

        if setup_only:
            global_log.info("setup_only: setup_only flag is set to True, exiting...")
            return

    else:

        global_paths['step8_grompp_min']['input_gro_path'] = input_gro_path
        global_paths['step8_grompp_min']['input_top_zip_path'] = input_top_path

    # STEP 8: minimization pre-processing
    global_log.info("step8_grompp_min: Preprocess energy minimization")
    grompp(**global_paths["step8_grompp_min"], properties=global_prop["step8_grompp_min"])

    # STEP 9: minimization
    global_log.info("step9_mdrun_min: Execute energy minimization")
    mdrun(**global_paths["step9_mdrun_min"], properties=global_prop["step9_mdrun_min"])

    # STEP 9B: create index file
    global_log.info("step9B_make_ndx: Create index file")
    make_ndx(**global_paths["step9B_make_ndx"], properties=global_prop["step9B_make_ndx"])

    # STEP 10: dump potential energy evolution
    global_log.info("step10_energy_min: Compute potential energy during minimization")
    gmx_energy(**global_paths["step10_energy_min"], properties=global_prop["step10_energy_min"])

    # STEP 11: NVT equilibration pre-processing
    global_log.info("step11_grompp_nvt: Preprocess system temperature equilibration")
    grompp(**global_paths["step11_grompp_nvt"], properties=global_prop["step11_grompp_nvt"])

    # STEP 12: NVT equilibration
    global_log.info("step12_mdrun_nvt: Execute system temperature equilibration")
    mdrun(**global_paths["step12_mdrun_nvt"], properties=global_prop["step12_mdrun_nvt"])

    # STEP 13: dump temperature evolution
    global_log.info("step13_temp_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step13_temp_nvt"], properties=global_prop["step13_temp_nvt"])

    # STEP 14: NPT equilibration pre-processing
    global_log.info("step14_grompp_npt: Preprocess system pressure equilibration")
    grompp(**global_paths["step14_grompp_npt"], properties=global_prop["step14_grompp_npt"])

    # STEP 15: NPT equilibration
    global_log.info("step15_mdrun_npt: Execute system pressure equilibration")
    mdrun(**global_paths["step15_mdrun_npt"], properties=global_prop["step15_mdrun_npt"])

    # STEP 16: dump density and pressure evolution
    global_log.info("step16_density_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step16_density_npt"], properties=global_prop["step16_density_npt"])

    if num_replicas:
        # Folder names for replicas
        simulation_folders = [f"replica_{i}" for i in range(int(num_replicas))]
    elif num_parts:
        # Folder names for parts
        simulation_folders = [f"parts_{i}" for i in range(int(num_parts))]

    global_log.info(f"Number of parts: {num_parts}")
    global_log.info(f"Number of replicas: {num_replicas}")

    # Run each simulation (replica or part)
    traj_list = []
    for simulation in simulation_folders:
        
        traj_prop = conf.get_prop_dic(prefix=simulation)
        traj_paths = conf.get_paths_dic(prefix=simulation)

        # Set general properties for all steps
        set_general_properties(traj_prop, conf, global_log)

        # Update previous global paths needed by simulation-specific steps
        traj_paths['step17_grompp_md']['input_gro_path'] = global_paths["step17_grompp_md"]['input_gro_path']
        traj_paths['step17_grompp_md']['input_cpt_path'] = global_paths["step17_grompp_md"]['input_cpt_path']
        traj_paths['step17_grompp_md']['input_top_zip_path'] = global_paths["step17_grompp_md"]['input_top_zip_path']
        traj_paths['step17_grompp_md']['input_ndx_path'] = global_paths["step17_grompp_md"]['input_ndx_path']
        traj_paths['step19_rmsd_equilibrated']['input_structure_path'] = global_paths["step15_mdrun_npt"]['output_gro_path']
        traj_paths['step20_rmsd_experimental']['input_structure_path'] = global_paths["step8_grompp_min"]['input_gro_path']
        traj_paths['step21_rgyr']['input_structure_path'] = global_paths["step15_mdrun_npt"]['output_gro_path']
        traj_paths['step22_rmsf']['input_top_path'] = global_paths["step7B_gro2pdb"]['output_str_path']
        
        # Enforce nsteps if provided
        if nsteps is not None:
            traj_prop['step17_grompp_md']['mdp']['nsteps']=int(nsteps)
            
        # Simulations are replicas
        if num_replicas:
            
            # Change seed and velocities for each replica
            traj_prop['step17_grompp_md']['mdp']['ld-seed'] = random.randint(1, 1000000)
            traj_prop['step17_grompp_md']['mdp']['continuation'] = 'no'
            traj_prop['step17_grompp_md']['mdp']['gen-vel'] = 'yes'

        # Simulations are parts of a single trajectory
        if num_parts:
            
            # Divide the number of steps by the number of parts
            traj_prop['step17_grompp_md']['mdp']['nsteps']=int(traj_prop['step17_grompp_md']['mdp']['nsteps']/int(num_parts))
            
            # For all parts except the first one, use the previous gro and cpt files
            if simulation != simulation_folders[0]:
                traj_paths['step17_grompp_md']['input_gro_path'] = previous_gro_path
                traj_paths['step17_grompp_md']['input_cpt_path'] = previous_cpt_path

        # STEP 17: free NPT production run pre-processing
        global_log.info(f"{simulation} >  step17_grompp_md: Preprocess free dynamics")
        grompp(**traj_paths['step17_grompp_md'], properties=traj_prop["step17_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{simulation} >  step18_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step18_mdrun_md'], properties=traj_prop['step18_mdrun_md'])
        
        # Append the trajectory to the list
        traj_list.append(traj_paths['step18_mdrun_md']['output_trr_path'])
        
        # Update the previous gro and cpt files
        previous_gro_path = traj_paths['step18_mdrun_md']['output_gro_path']
        previous_cpt_path = traj_paths['step18_mdrun_md']['output_cpt_path']
        
        # STEP 19: compute the RMSD with respect to equilibrated structure
        global_log.info("step19_rmsd_equilibrated: Compute Root Mean Square deviation against equilibrated structure")
        gmx_rms(**traj_paths['step19_rmsd_equilibrated'], properties=traj_prop['step19_rmsd_equilibrated'])
        
        # STEP 20: compute the RMSD with respect to minimized structure
        global_log.info("step20_rmsd_experimental: Compute Root Mean Square deviation against minimized structure (exp)")
        gmx_rms(**traj_paths['step20_rmsd_experimental'], properties=traj_prop['step20_rmsd_experimental'])
        
        # STEP 21: compute the Radius of gyration
        global_log.info("step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
        gmx_rgyr(**traj_paths['step21_rgyr'], properties=traj_prop['step21_rgyr'])

        # STEP 22: compute the RMSF
        global_log.info("step22_rmsf: Compute Root Mean Square Fluctuation to measure the protein flexibility during the free MD simulation")
        cpptraj_rmsf(**traj_paths['step22_rmsf'], properties=traj_prop['step22_rmsf'])
        
    # Do the final analysis with all the previous parts or replicas
    if final_analysis:
        
        # If simulations are different parts of a single trajectory
        if num_parts:
            
            # Concatenate the analysis files that can be concatenated
            concatenate_gmx_analysis(conf, simulation_folders, output_path)
        
            # STEP 23: concatenate trajectories
            global_log.info("step23_trjcat: Concatenate trajectories")
            fu.zip_list(zip_file=global_paths["step23_trjcat"]['input_trj_zip_path'], file_list=traj_list)
            trjcat(**global_paths["step23_trjcat"], properties=global_prop["step23_trjcat"])
            
            # STEP 24: obtain dry the merged trajectory
            global_log.info("step24_dry_trj: Obtain dry trajectory")
            gmx_trjconv_trj(**global_paths["step24_dry_trj"], properties=global_prop["step24_dry_trj"])
        
            #Remove unused trajectory
            os.remove(global_paths["step23_trjcat"]["output_trj_path"])
            
            # STEP 25: obtain dry structure
            global_log.info("step25_dry_str: Obtain dry structure")
            gmx_trjconv_str(**global_paths["step25_dry_str"], properties=global_prop["step25_dry_str"])

            # STEP 26: image the trajectory
            global_log.info("step26_image_traj: Imaging the trajectory")
            gmx_image(**global_paths['step26_image_traj'], properties=global_prop['step26_image_traj'])
            
            #Remove unused trajectory
            os.remove(global_paths["step24_dry_trj"]["output_traj_path"])

            # STEP 27: fit the trajectory
            global_log.info("step27_fit_traj: Fit the trajectory")
            gmx_image(**global_paths['step27_fit_traj'], properties=global_prop['step27_fit_traj'])
            
            #Remove unused trajectory
            os.remove(global_paths["step26_image_traj"]["output_traj_path"])
            
        # If simulations are replicas
        if num_replicas:
            
            # For each replica, do the final analysis
            for simulation in simulation_folders:
                
                # Get the properties and paths for the replica
                traj_prop = conf.get_prop_dic(prefix=simulation)
                traj_paths = conf.get_paths_dic(prefix=simulation)
                
                # Set general properties for all steps
                set_general_properties(traj_prop, conf, global_log)
                
                # NOTE: we are hard-coding the kind of traj that we are using with these paths: output_trr_path
                # Update previous global paths needed by simulation-specific steps
                traj_paths['step24_dry_trj']['input_traj_path'] = traj_paths['step18_mdrun_md']['output_trr_path']
                traj_paths['step24_dry_trj']['input_top_path'] = global_paths["step15_mdrun_npt"]['output_gro_path']
                traj_paths['step24_dry_trj']['input_index_path'] = global_paths["step9B_make_ndx"]['output_ndx_path']
                traj_paths['step25_dry_str']['input_structure_path'] = global_paths["step15_mdrun_npt"]['output_gro_path']
                traj_paths['step25_dry_str']['input_top_path'] = global_paths["step15_mdrun_npt"]['output_gro_path']
                traj_paths['step25_dry_str']['input_index_path'] = global_paths["step9B_make_ndx"]['output_ndx_path']
                traj_paths['step26_image_traj']['input_index_path'] = global_paths["step9B_make_ndx"]['output_ndx_path']
                traj_paths['step27_fit_traj']['input_index_path'] = global_paths["step9B_make_ndx"]['output_ndx_path']
                
                # STEP 24: obtain dry the trajectory
                global_log.info(f"{simulation} > step24_dry_trj: Obtain dry trajectory")
                gmx_trjconv_trj(**traj_paths["step24_dry_trj"], properties=traj_prop["step24_dry_trj"])
                
                # STEP 25: obtain dry structure
                global_log.info(f"{simulation} > step25_dry_str: Obtain dry structure")
                gmx_trjconv_str(**traj_paths["step25_dry_str"], properties=traj_prop["step25_dry_str"])
                
                # STEP 26: image the trajectory
                global_log.info(f"{simulation} > step26_image_traj: Imaging the trajectory")
                gmx_image(**traj_paths['step26_image_traj'], properties=traj_prop['step26_image_traj'])
                
                #Remove unused trajectory
                os.remove(traj_paths["step24_dry_trj"]["output_traj_path"])
                
                # STEP 27: fit the trajectory
                global_log.info(f"{simulation} > step27_fit_traj: Fit the trajectory")
                gmx_image(**traj_paths['step27_fit_traj'], properties=traj_prop['step27_fit_traj'])
                
                #Remove unused trajectory
                os.remove(traj_paths["step26_image_traj"]["output_traj_path"])

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

    parser = argparse.ArgumentParser("Simple MD Protein Setup")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system (default: False)",
                        required=False)

    parser.add_argument('--num_parts', dest='num_parts',
                        help="Number of parts of the trajectorie (default: 1)",
                        required=False)
    
    parser.add_argument('--num_replicas', dest='num_replicas',
                        help="Number of replicas (default: 1)",
                        required=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)

    parser.add_argument('--input_pdb', dest='input_pdb_path',
                        help="Input PDB file (default: input_structure_path in step 1 of configuration file)",
                        required=False)

    parser.add_argument('--his', dest='his',
                        help="Histidine protonation states list.",
                        required=False, type = str)

    parser.add_argument('--pdb_chains', nargs='+', dest='pdb_chains',
                        help="PDB chains to be extracted from PDB file (default: chains in properties of step 1)",
                        required=False)

    parser.add_argument('--mutation_list', nargs='+', dest='mutation_list',
                        help="List of mutations to be introduced in the protein (default: None, ex: A:Arg220Ala)",
                        required=False)

    parser.add_argument('--input_gro', dest='input_gro_path',
                        help="Input structure file ready to minimize (.gro). To provide an externally prepared system, use together with --input_top (default: None)",
                        required=False)

    parser.add_argument('--input_top', dest='input_top_path',
                        help="Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, use together with --input_gro (default: None)",
                        required=False)

    parser.add_argument('--fix_ss', action='store_true',
                        help="Add disulfide bonds to the protein. Use carefully! (default: False)",
                        required=False)

    parser.add_argument('--fix_amides', action='store_true', dest='fix_amide_clashes',
                        help="Flip clashing amides to relieve the clashes (default: False)",
                        required=False)

    parser.add_argument('--nsteps', dest='nsteps',
                        help="Number of steps of the simulation",
                        required=False)

    parser.add_argument('--final_analysis', action='store_true', dest='final_analysis',
                        help="Run the final analysis of the trajectory/ies. Concatenation of the analysis and trajectory, trajectory drying, imaging and fitting (default: False)",
                        required=False)

    args = parser.parse_args()

    # Check .gro structure and .zip topology are given together
    if (args.input_gro_path is None and args.input_top_path is not None) or (args.input_gro_path is not None and args.input_top_path is None):
        raise Exception("Both --input_gro and --input_top must be provided together")

    # Check .pdb structure and .gro/.zip topology are not given together
    if (args.input_pdb_path is not None and args.input_gro_path is not None):
        raise Exception("Both --input_pdb and --input_gro/--input_top are provided. Please provide only one of them")

    main_wf(configuration_path=args.config_path, setup_only=args.setup_only, num_parts=args.num_parts, num_replicas=args.num_replicas, output_path=args.output_path,
            input_pdb_path=args.input_pdb_path, pdb_chains=args.pdb_chains, mutation_list=args.mutation_list,
            input_gro_path=args.input_gro_path, input_top_path=args.input_top_path,
            fix_ss=args.fix_ss, fix_amide_clashes=args.fix_amide_clashes, his=args.his, nsteps=args.nsteps, final_analysis=args.final_analysis)