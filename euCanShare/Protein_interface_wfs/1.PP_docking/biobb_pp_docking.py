#!/usr/bin/env python3

# Importing all the needed libraries
import os
import re
import time
import shutil
import argparse
import numpy as np
from Bio import PDB
import pandas as pd
from pathlib import Path
from functools import partial
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.encore.confdistmatrix import conformational_distance_matrix
from MDAnalysis.analysis.encore.confdistmatrix import set_rmsd_matrix_elements

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_pydock.pydock.setup import setup
from biobb_pydock.pydock.ftdock import ftdock
from biobb_pydock.pydock.dockser import dockser
from biobb_pydock.pydock.makePDB import makePDB
from biobb_pydock.pydock.oda import oda

#############
# Functions #
#############

def launch_step_full(global_log, global_paths, output_path, previous_output_path, tool, step_name, description, global_prop, run_remaining_steps):
    """
    Function used to run one step of the workflow. It will check if the step is requested by the user and if it is, it will run it.

    Constant arguments (same for all steps):
    ----------------------------------------

        global_log            (Logger): logger object
        global_paths            (dict): dictionary with all the paths used in the workflow
        output_path              (str): path to output folder
        previous_output_path     (str): path to previous output folder
    
    Variable arguments (different for each step):
    --------------------------------------------

        tool                     (str): tool name
        step_name                (str): step name
        description              (str): step description
        global_prop             (dict): dictionary with all the properties used in the workflow
        run_remaining_steps     (bool): flag to run remaining steps
    
    Returns
    -------

        elapsed_step_time (float): elapsed time for the step
    """

    # Start time
    start_step_time = time.time()

    if run_remaining_steps:
        # Run step
        global_log.info(f"{step_name}: {description}")
        tool(**global_paths[step_name], properties=global_prop[step_name])
    else:
        # Link step folder in previous output path to output path
        link_previous_step(global_log, step_name, output_path, previous_output_path)

    # End time
    elapsed_step_time = time.time() - start_step_time

    return elapsed_step_time

def check_arguments(global_log, receptor_pdb_path, ligand_pdb_path, output_path, previous_output_path, skip_until):
    """
    Check the arguments provided by the user
    """

    # If the user provides receptor and ligands as well as previous output path, raise a warning
    if ((receptor_pdb_path is not None) or (ligand_pdb_path is not None)) and (previous_output_path is not None):
        global_log.warning("WARNING: You provided receptor and ligand PDB files as well as a previous output path. The receptor and ligand PDB files will not be used!")

    # If the user provides skip_until but not a previous output path, raise an error
    if (skip_until is not None) and (previous_output_path is None):
        global_log.error("ERROR: You provided --skip_until but not a --previous_output. Please provide a previous output path!")

    # If the user provides a previous output path but not skip_until, raise an error
    if (skip_until is None) and (previous_output_path is not None):
        global_log.error("ERROR: You provided --previous_output but not --skip_until. Please provide a step to skip until!")

    # If skip_until is provided, check if it is a valid option (makePDB, clustering, oda_filtering)
    if skip_until is not None:
        if skip_until not in ["makePDB", "clustering", "oda_filtering"]:
            global_log.error("ERROR: --skip_until must be one of the following options: makePDB, clustering, oda_filtering!")

    if previous_output_path is not None:

        # Check if previous output path exists
        if not os.path.exists(previous_output_path):
            global_log.error("ERROR: The previous output path provided does not exist!")

        # Check if previous output path is a folder
        if not os.path.isdir(previous_output_path):
            global_log.error("ERROR: The previous output path provided is not a folder!")

        # Check if output path and previous output path are the same
        if output_path == previous_output_path:
            global_log.error("ERROR: The output path and the previous output path are the same!")

def link_previous_step(global_log, step_name, output_path, previous_output_path):
    """
    Link step folders in previous output path to output path

    Inputs
    ------

        step_name              (str): step name to link
        output_path            (str): path to output folder
        previous_output_path   (str): path to previous output folder
    """

    # Check if previous output path is provided
    if previous_output_path is None:
        return

    # Check if step folder exists in previous output path
    if not os.path.exists(str(Path(previous_output_path).joinpath(step_name))):
        global_log.warning(f"WARNING: The step folder {step_name} does not exist in the previous output path provided!")
        return
    
    # Check if step folder exists in output path
    if os.path.exists(str(Path(output_path).joinpath(step_name))):
        global_log.warning(f"WARNING: The step folder {step_name} already exists in the output path provided!")
        return
    
    # Link step folder in previous output path to output path
    os.symlink(str(Path(previous_output_path).joinpath(step_name)), str(Path(output_path).joinpath(step_name)))
    
    return

def prepare_step(input_zip_path: str, prop: dict):
    """
    1. Create step folder
    2. Create docking poses subfolder and unzip docking poses into it

    Inputs
    ------

        input_zip_path          (str): path to zip file with docking poses
        prop                   (dict): dictionary with oda_filtering step properties
    
    Returns
    -------

        poses_paths (list): list with paths to extracted docking poses
    """

    # Create oda filtering folder
    fu.create_dir(prop["path"])

    # Create subfolder for poses
    poses_folder_path = fu.create_dir(str(Path(prop["path"]).joinpath("poses")))

    # Extract the file with the docking poses into the subfolder
    poses_paths = fu.unzip_list(zip_file = input_zip_path, dest_dir = poses_folder_path)

    return poses_paths

def find_chain(structure, chain_id):
    """
    Find the Bio.PDB.Chain object with the given chain ID.

    Inputs
    ------

        structure (Bio.PDB.Structure): structure object
        chain_id               (str): chain ID
    
    Returns
    -------

        chain (Bio.PDB.Chain): chain object
    """

    # Iterate over all chains in the structure
    for chain in structure.get_chains():

        # If the chain ID matches the given chain ID, return the chain
        if chain.id == chain_id:
            return chain

    # If the chain ID does not match any chain ID, raise an error
    raise ValueError(f"Chain {chain_id} not found in structure!")


def rmsd_clustering(input_zip_path: str, output_zip_path: str, properties: dict):

    """ 
    Cluster the docking poses based on RMSD and save the best ranking pose from each cluster into a zip file.

    The input zip file contains docking poses in PDB format. The function will calculate the RMSD matrix
    between all poses using only the CA atoms of the ligand (as the receptor didn't move during docking). 
    Then it will calculate the linkage matrix and the clusters for a given RMSD threshold. Finally it will
    find the best ranking pose for each cluster and save it into a zip file in PDB format.

    Inputs
    ------

        input_zip_path          (str): path to zip file with docking poses
        output_zip_path         (str): path to output zip file with best ranking docking poses
        properties             (dict): dictionary with clustering step properties
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return
    
    # Prepare clustering step folder 
    poses_paths = prepare_step(input_zip_path, properties)

    # Load the ranking dataframe
    ranking_df = properties["ranking_df"]

    # Rename docking poses with rank and create a list with the corresponding rank for each file
    poses_ranked_paths, paths_ranking = rename_with_rank(poses_paths, ranking_df)

    # Calculate the RMSD matrix between all poses
    rmsd_matrix = compute_rmsd_matrix(poses_ranked_paths, properties["ligand_chain"])

    # Convert the uncondensed RMSD matrix into a condensed distance matrix
    rmsd_matrix = squareform(rmsd_matrix)

    # Calculate the linkage matrix
    Z = linkage(rmsd_matrix, method = "average")

    # Get the cluster labels from the linkage matrix
    cluster_labels = fcluster(Z, t = properties["rmsd_threshold"], criterion = "distance")

    plot_dendogram(Z, paths_ranking, properties["rmsd_threshold"], properties["path"])

    # Find the best ranking pose path for each cluster
    top_distinct_poses_paths = []
    for cluster in range(1, max(cluster_labels) + 1):

        # Get the ranking of the poses in the cluster
        cluster_ranking = [paths_ranking[i] for i in np.where(cluster_labels == cluster)[0]]

        # Get the paths of the poses in the cluster
        cluster_paths = [poses_ranked_paths[i] for i in np.where(cluster_labels == cluster)[0]]

        # Find the path of the pose with the best ranking in the cluster
        best_pose_path = cluster_paths[cluster_ranking.index(min(cluster_ranking))]
        top_distinct_poses_paths.append(best_pose_path)

    # Zip the best ranking pose paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = top_distinct_poses_paths)

    # Remove temporal folder with poses
    fu.rm(os.path.join(properties["path"], "poses"))

def compute_rmsd_matrix(poses_ranked_paths: list, ligand_chain: str):
    """
    Computes an uncondensed RMSD matrix between all poses using only the CA atoms of the ligand (as the receptor didn't move during docking).

    Inputs
    ------

        poses_ranked_paths (list): list with paths to docking poses
        ligand_chain        (str): chain of the ligand protein
    
    Returns
    -------

        rmsd_matrix (np.array): uncondensed RMSD matrix between all poses 
    """

    # Find the number of cores to use for parallelization (thread-based parallelization with shared memory, see conformational_distance_matrix source code)
    if os.getenv('SLURM_CPUS_PER_TASK') is None:
        # Running in local environment, use all available cores
        num_threads = os.cpu_count()
    else:
        # Running on HPC environment, use SLURM's allocated CPUs
        num_threads = int(os.getenv('SLURM_CPUS_PER_TASK'))

    # Max number of CA atoms to use for RMSD calculation 
    # (enough to describe the rigid ligand orientation)
    max_ca_atoms = 20

    # Create a universe with all poses, use the first pose as the topology
    poses_universe = mda.Universe(poses_ranked_paths[0])

    # Read coordinates of all poses into the universe
    poses_universe.load_new(poses_ranked_paths)

    # Select the CA atoms of the ligand
    ligand_ca_atoms = poses_universe.select_atoms(f"protein and chainID {ligand_chain} and name CA")

    # Find the total number of CA atoms in the ligand
    n_ligand_ca_atoms = len(ligand_ca_atoms)

    # If the number of CA atoms in the ligand is greater than the maximum number of CA atoms to use for RMSD calculation
    if n_ligand_ca_atoms > max_ca_atoms:

        # Create a list with max_ca_atoms indices uniformly spaced from 0 to n_ligand_ca_atoms
        indices = np.linspace(0, n_ligand_ca_atoms - 1, max_ca_atoms, dtype = int)

        # Select the CA atoms to use for RMSD calculation from the ligand_ca_atoms AtomGroup
        ligand_ca_atoms = ligand_ca_atoms[indices]

    # Transform the ligand_ca_atoms AtomGroup into a selection string
    rmsd_selection = "index " + " or index ".join([str(atom.index) for atom in ligand_ca_atoms])

    # Compute the RMSD matrix between all poses using only the CA atoms of the ligand
    rmsd_triangular_matrix = conformational_distance_matrix(poses_universe, 
                                                 conf_dist_function=set_rmsd_matrix_elements,
                                                 select = rmsd_selection,
                                                 pairwise_align = False,
                                                 metadata = False,
                                                 n_jobs = num_threads)

    return rmsd_triangular_matrix.as_array()

def rename_with_rank(poses_paths: list, ranking_df: pd.DataFrame):
    """ 
    Change pdb poses file names from conformation number to rank.

    Inputs
    ------

        poses_paths         (list): list with paths to docking poses
        ranking_df  (pd.DataFrame): dataframe with ranking of docking poses

    Returns
    -------

        poses_ranked_paths (list): list with paths to docking poses with ranked file names
        paths_ranking      (list): list with the corresponding rank for each docking pose
    
    """
    paths_ranking = []
    poses_ranked_paths = []

    # Find the parent path to the poses
    poses_parent_path = str(Path(poses_paths[0]).parent)

    # Change pdb file name for rank and create a list with the ranking
    for pdb_path in poses_paths:

        # Get the pdb name without the extension
        pdb_name = Path(pdb_path).stem

        # Look for a number (e.g. 13 or 123) in the pdb name    
        number_match = re.findall(r"\d+", pdb_name)

        # Get the conformation number from the match
        conformation_number = number_match[0]

        # Get the rank from the ranking Dataframe
        rank = ranking_df[ranking_df["Conf"] == int(conformation_number)]["RANK"].values[0]
        paths_ranking.append(int(rank))

        # Rename the file with the ranking
        new_pdb_path = str(Path(poses_parent_path).joinpath(f"rank_{rank}.pdb"))
        os.rename(pdb_path, new_pdb_path)

        # Add the new file path to the list
        poses_ranked_paths.append(new_pdb_path)
    
    return poses_ranked_paths, paths_ranking

def plot_dendogram(linkage_matrix, labels, rmsd_threshold, output_path):
    """ Plot the dendogram and save it as a png file """

    # Plot the dendrogram and color the leaves according to the clusters
    plt.figure(figsize = (25, 10))
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Rank of pose", fontsize=16)
    plt.ylabel("RMSD (Å)", fontsize=16)
    dendrogram(
        linkage_matrix,
        leaf_rotation = 90.,  # rotates the x axis labels
        leaf_font_size = 8.,  # font size for the x axis labels
        color_threshold = rmsd_threshold,
        labels = labels
        )
    plt.savefig(str(Path(output_path).joinpath("dendrogram.png")))
    plt.close() 


def oda_filtering(input_receptor_path: str, input_ligand_path: str, input_zip_path: str, output_zip_path: str, properties: dict):
    """
    Function that filters docking poses using ODA patches and saves the filtered poses into a zip file in PDB format.

    1. Find all residues pertaining to the surface of the receptor and ligand proteins (through their SASA values).
    2. Identify interface residues between the receptor and ligand proteins for each pose.
    3. Compute the percentage of interface residues that are covered by ODA patches for each pose. 
    4. Keep only the poses with a percentage of interface residues covered by ODA patches above a given threshold.

    Inputs
    ------

        input_receptor_path     (str): path to receptor pdb file with ODA values
        input_ligand_path       (str): path to ligand pdb file with ODA values
        input_zip_path          (str): path to zip file with docking poses
        output_zip_path         (str): path to output zip file with filtered docking poses
        properties             (dict): dictionary with oda_filtering step properties
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return
    
    # Prepare step folder and unzip docking poses into it
    poses_paths = prepare_step(input_zip_path, properties)

    # Find surface residues: dictionaries with residue ID as key and their CA atoms / ODA values as values
    receptor_surface_atoms, receptor_surface_oda = find_surface(input_receptor_path, properties["path"])
    ligand_surface_atoms, ligand_surface_oda = find_surface(input_ligand_path, properties["path"])

    # Create the Neighbor Search object for the receptor surface atoms (receptor remains fixed during docking) 
    Receptor_Neighbor_Search = PDB.NeighborSearch(list(receptor_surface_atoms.values())) # NOTE: optimize bucket size?

    # For each pose
    filtered_poses_paths = []
    for pose_path in poses_paths:

        # Identify interface residues for the pose
        receptor_interface_atoms, ligand_interface_atoms = find_interface_residues(pose_path, ligand_surface_oda.keys(), Receptor_Neighbor_Search, properties)

        # Filter the pose 
        accept = filter_pose_interface(receptor_interface_atoms, ligand_interface_atoms, 
                                  receptor_surface_oda, ligand_surface_oda, properties)

        if accept:
            filtered_poses_paths.append(pose_path)
    
    # Zip the filtered docking pose paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = filtered_poses_paths)

    # Remove temporal folders with poses and decorated poses
    fu.rm(os.path.join(properties["path"], "poses"))

def decorate_with_oda_values(pose_path, input_receptor_path, input_ligand_path, oda_poses_folder_path, receptor_chain, ligand_chain):
    """
    Takes a docking pose and decorates the receptor and ligand proteins with the corresponding ODA values.

    Inputs
    ------

        pose_path               (str): path to docking pose (contains receptor and ligand)
        input_receptor_path     (str): path to receptor pdb file with ODA values in the B-factor column
        input_ligand_path       (str): path to ligand pdb file with ODA values in the B-factor column
        oda_poses_folder_path   (str): path to folder where the decorated docking pose will be saved
        receptor_chain          (str): chain of the receptor protein
        ligand_chain            (str): chain of the ligand protein
    """

    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Load the docking pose
    pose_structure = parser.get_structure("pose", pose_path)

    # Select the receptor and ligand atoms from the docking pose
    receptor_atoms = []
    ligand_atoms = []

    for model in pose_structure:
        for chain in model:
            if chain.id == receptor_chain:
                receptor_atoms.extend(chain.get_atoms())
            elif chain.id == ligand_chain:
                ligand_atoms.extend(chain.get_atoms())

    # Load the receptor and ligand with the ODA values
    receptor_structure = parser.get_structure("receptor_oda", input_receptor_path)
    ligand_structure = parser.get_structure("ligand_oda", input_ligand_path)

    # Get the atoms with ODA values from the receptor and ligand structures
    receptor_atoms_oda = [atom for atom in receptor_structure.get_atoms() if atom.bfactor is not None]
    ligand_atoms_oda = [atom for atom in ligand_structure.get_atoms() if atom.bfactor is not None]

    # Iterate over all ligand atoms and decorate them with the corresponding ODA value
    for ligand_atom, ligand_oda_atom in zip(ligand_atoms, ligand_atoms_oda):
        ligand_atom.set_bfactor(ligand_oda_atom.bfactor)

    # Iterate over all receptor atoms and decorate them with the corresponding ODA value
    for receptor_atom, receptor_oda_atom in zip(receptor_atoms, receptor_atoms_oda):
        receptor_atom.set_bfactor(receptor_oda_atom.bfactor)

    # Find the name of the docking pose file
    pose_name = Path(pose_path).name

    # Create new path for the decorated docking pose
    decorated_pose_path = os.path.join(oda_poses_folder_path, pose_name)

    # Save the decorated docking pose
    io = PDB.PDBIO()
    io.set_structure(pose_structure)
    io.save(decorated_pose_path)

    return decorated_pose_path

def find_surface(pdb_path: str, output_path: str):
    """
    Find all surface residues of a pdb structure (using the Solvent Accessible Surface Area and a threshold).

    Return a list with the CA atoms of the surface residues (Bio.PDB.Atom.Atom objects) and a dictionary 
    containing the residue IDs of the surface residues (keys) and their ODA values (values). Here we are assuming the 
    beta factor column of the pdb file contains the ODA values.

    (Optionally) Save a structure showing the surface residues in the output path (for debugging purposes).

    Inputs
    ------

        pdb_path (str): path to pdb file

    Returns
    -------

        ca_atoms         (dict): dictionary with residue IDs of surface residues (keys) and their CA atoms (Bio.PDB.Atom.Atom objects) (values)
        surface_oda      (dict): dictionary with residue IDs of surface residues (keys) and their ODA values (values)
    """
    # SASA threshold to classify residues as surface or not
    sasa_threshold = 10 

    # Load the pdb structure
    structure = PDB.PDBParser(QUIET=True).get_structure("pdb", pdb_path)

    # Calculate the SASA for each residue
    sasa_calc = PDB.SASA.ShrakeRupley(probe_radius=1.4, n_points=100)
    sasa_calc.compute(structure, level="R")

    # Classify residues as surface or not
    ca_atoms = {}
    surface_oda = {}
    for residue in structure.get_residues():

        if residue.sasa > sasa_threshold:
            
            # Get the CA atom (Bio.PDB.Atom.Atom object) from residue (Bio.PDB.Residue.Residue object)
            ca_atom = residue["CA"]
            
            # Add the CA atom to the list of surface CA atoms
            ca_atoms[residue.get_id()] = ca_atom

            # Add the residue ID and ODA value to the dictionary
            surface_oda[residue.get_id()] = ca_atom.bfactor

    # Check if any surface residue was found
    if len(ca_atoms) == 0:
        raise ValueError("No surface residues found!")
    
    # Save a structure showing the surface residues (for debugging purposes) NOTE: make optional
    debug_find_surface(pdb_path, sasa_threshold, output_path)

    return ca_atoms, surface_oda

def find_interface_residues(pose_path: str, ligand_surface_ids: list, Receptor_Neighbor_Search: PDB.NeighborSearch, properties: dict):
    """
    Find the receptor and ligand interface residues for a given docking pose.

    Iterates over the ligand surface CA atoms and finds the receptor surface CA atoms that are within a given distance threshold.
    
    Inputs
    ------

        pose_path                               (str): path to docking pose with receptor and ligand
        ligand_surface_ids                     (list): list with residue IDs of ligand surface residues
        Receptor_Neighbor_Search (PDB.NeighborSearch): class with the receptor surface atoms to search for neighbors
        properties                             (dict): dictionary with step properties

    Returns
    -------

        receptor_interface_atoms (dict): dictionary with residue IDs of receptor interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
        ligand_interface_atoms   (dict): dictionary with residue IDs of ligand interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
    """
    # Global log
    global_log = properties["global_log"]

    # Load the docking pose using BioPython
    parser = PDB.PDBParser(QUIET=True)
    pose_structure = parser.get_structure("pose", pose_path)

    # Find the ligand chain
    ligand_chain = find_chain(pose_structure, properties["ligand_chain"])
    
    # Dictionaries with interface atoms
    receptor_interface_atoms = {}
    ligand_interface_atoms = {}

    # Iterate over all ligand surface atoms
    for ligand_residue in ligand_chain.get_residues():

        # Check if the ligand residue is a surface residue
        if ligand_residue.get_id() in ligand_surface_ids:

            # Get the ligand CA atom 
            residue_ca_atom = ligand_residue["CA"]

            # Find the receptor surface residues within a given distance threshold from the ligand surface residue
            receptor_neighbors = Receptor_Neighbor_Search.search(residue_ca_atom.get_coord(), properties["distance_threshold"], level="R")

            # Check if any receptor surface residue was found
            if len(receptor_neighbors) > 0:

                # Add receptor surface CA atoms to the interface atoms
                for receptor_residue in receptor_neighbors:
                    
                    # Check if we added the CA atom of the residue already
                    if receptor_residue.get_id() not in receptor_interface_atoms:
                        receptor_interface_atoms[receptor_residue.get_id()] = receptor_residue["CA"]

                # Add ligand surface CA atom to the interface atoms
                ligand_interface_atoms[ligand_residue.get_id()] = residue_ca_atom

    # Find the number of interface residues
    num_interface_residues = len(receptor_interface_atoms) + len(ligand_interface_atoms)

    # Check if any interface residue was found
    if num_interface_residues == 0:
        global_log.error("WARNING: No interface residues were found! Check the distance threshold and the ligand chain!")
    
    return receptor_interface_atoms, ligand_interface_atoms

def filter_pose_interface(receptor_interface_atoms, ligand_interface_atoms, receptor_surface_oda, ligand_surface_oda, properties: dict):
    """
    Filter the pose according to the ODA patches and the specific interface. Two main criteria:

        1. PIR COP above threshold (compulsory)
            
            Compute the "PIR COP" (Percentage of Interface Residues Covered by ODA Patches) for a given interface.
            This measures the overlap between the interface and the ODA patches. If the overlap is above a given threshold,
            the pose is accepted.
        
        2. Percentage of neighboring interface residues covered by oda patches between receptor and ligand.
        
            This measures overlap of ligand and receptor oda patches inside the interface. If the overlap is above a given threshold,
            the pose is accepted. This is optional and can be turned off by setting the threshold to 0 for faster filtering.

    Inputs
    ------

        receptor_interface_atoms (dict): dictionary with residue IDs of receptor interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
        ligand_interface_atoms   (dict): dictionary with residue IDs of ligand interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
        receptor_surface_oda     (dict): dictionary with residue IDs of receptor surface residues as keys and their ODA values as values
        ligand_surface_oda       (dict): dictionary with residue IDs of ligand surface residues as keys and their ODA values as values
        properties               (dict): dictionary with step properties
    
    Returns
    -------

        accept                   (bool): True if the pose is accepted, False otherwise
    """

    # Get the total number of interface residues
    n_interface_residues = len(receptor_interface_atoms) + len(ligand_interface_atoms)

    # Count the number of residues in the interface covered by ODA patches
    n_interface_residues_covered = 0

    # Count the number of residues in the interface covered by ODA patches and overlapping between receptor and ligand
    n_overlapping_residues = 0

    if properties["overlap_threshold"] > 0:

        # Create the Neighbor Search for the receptor interface surface atoms 
        Receptor_Neighbor_Search = PDB.NeighborSearch(list(receptor_interface_atoms.values())) # NOTE: optimize bucket size?

        # Create the Neighbor Search for the ligand interface surface atoms
        Ligand_Neighbor_Search = PDB.NeighborSearch(list(ligand_interface_atoms.values())) # NOTE: optimize bucket size?

    # Get the ODA threshold
    oda_threshold = properties["oda_threshold"]

    # Iterate over the receptor interface residues
    for residue_id, residue_ca_atom in receptor_interface_atoms.items():

        # Get the ODA value of the residue
        residue_oda = receptor_surface_oda[residue_id]

        # Check if the residue is covered by an ODA patch
        if residue_oda >= oda_threshold:
            n_interface_residues_covered += 1

    # Iterate over the ligand interface residues
    for residue_id in ligand_interface_atoms.keys():
        # Get the ODA value of the residue
        residue_oda = ligand_surface_oda[residue_id]

        # Check if the residue is covered by an ODA patch
        if residue_oda >= oda_threshold:
            n_interface_residues_covered += 1

    # Compute the "PIR COP" (Percentage of Interface Residues Covered by ODA Patches)
    pir_cop = n_interface_residues_covered / n_interface_residues

    # Check if the PIR COP is above the threshold
    if pir_cop < properties["pir_cop_threshold"]:
        return False
    
    # Get the number of neighboring interface residues
    return True

def debug_find_surface(pdb_path: str, sasa_threshold: float, output_path: str):
    """
    Function to debug residue surface finding. Saves a pdb file with the surface residues marked in the beta column.
    """

    # Load the pdb structure
    structure = PDB.PDBParser(QUIET=True).get_structure("pdb", pdb_path)

    # Find the pdb name
    pdb_name = Path(pdb_path).name

    # Calculate the SASA for each residue
    sasa_calc = PDB.SASA.ShrakeRupley(probe_radius=1.4, n_points=100)
    sasa_calc.compute(structure, level="R")

    for residue in structure.get_residues():

        if residue.sasa > sasa_threshold:
            # Mark surface residue
            for atom in residue:
                atom.set_bfactor(100)

        else:
            # Mark non-surface residue
            for atom in residue:
                atom.set_bfactor(0)

    # Save the structure with the surface residues marked in the beta column
    io = PDB.PDBIO()
    io.set_structure(structure)
    file_name = f"{pdb_name}_surface.pdb"
    file_path = os.path.join(output_path, file_name)
    io.save(file_path)


def distance_filtering(input_zip_path: str, output_zip_path: str, properties: dict):
    """
    Function that filters docking poses according to the distance between pairs of residues 
    (one from the receptor and one from the ligand) and saves the filtered poses into a zip file 
    in PDB format.

    1. Find the distance between pairs of residues (one from the receptor and one from the ligand)
    2. Keep only the poses with a distance between pairs of residues below a given threshold.

    Inputs
    ------
    
        input_zip_path          (str): path to zip file with docking poses
        output_zip_path         (str): path to output zip file with filtered docking poses
        properties             (dict): dictionary with distance_filtering step properties
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return

    # Prepare step folder and unzip docking poses into it
    poses_paths = prepare_step(input_zip_path, properties)
    
    # For each pose in chunk
    filtered_poses_paths = []
    for pose_path in poses_paths:

        # Find the distance between pairs of residues (one from the receptor and one from the ligand)
        accept_pose = calculate_distances(pose_path, properties)

        # Keep only the poses with a distance between pairs of residues below a given threshold.
        if accept_pose:
            filtered_poses_paths.append(pose_path)

    # Save the filtered poses paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = filtered_poses_paths)

    # Remove temporal folder with poses
    fu.rm(os.path.join(properties["path"], "poses"))

def calculate_distances(pose_path: str, properties: dict):
    """
    Compute distances between pairs of residues defined in properties and check
    they are below a given threshold. If all distances are below the threshold,
    return True. Otherwise return False.

    Inputs
    ------

        pose_path (str): path to docking pose
        properties (dict): dictionary with step properties
    
    Return
    ------

        accept_pose (bool): True if all distances are below the threshold, False otherwise
    """

    # Read the docking pose
    pose_universe = mda.Universe(pose_path)

    # Get the receptor and ligand chains 
    receptor_chain = properties["receptor_chain"]
    ligand_chain = properties["ligand_chain"]

    # For each distance in the properties
    for distance in properties["distances"]:

        # Get the receptor residue CA atom
        receptor_residue = pose_universe.select_atoms(f"protein and chainID {receptor_chain} and {distance['receptor_residue_selection']} and name CA")

        # Get the ligand residue CA atom
        ligand_residue = pose_universe.select_atoms(f"protein and chainID {ligand_chain} and {distance['ligand_residue_selection']} and name CA")

        # Get the distance between the receptor and ligand CA atoms
        distance_value = distance_array(receptor_residue.positions, ligand_residue.positions)[0][0]

        # Check if the distance is below the threshold
        if distance_value > distance["threshold"]:
            return False
    
    return True

########
# Main #
########

def main_wf(configuration_path, receptor_pdb_path, ligand_pdb_path, previous_output_path, skip_until, output_path):
    '''
    Main protein-protein docking workflow. Can be used to sample protein-protein docking poses,
    rank them according to pyDock score and cluster them based on RMSD. The best scoring pose from each cluster 
    is saved into a zip file in PDB format.

    Inputs
    ------

        configuration_path   (str): path to input.yml
        receptor_pdb_path    (str): path to receptor pdb file
        ligand_pdb_path      (str): path to ligand pdb file
        previous_output_path (str): path to previous output folder. Used if you want to re-use some steps of a previous run
        skip_until           (str): skip everything until this step. Options: dockrst, makePDB, clustering, oda_filtering 
        output_path          (str): path to output folder

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    # Start timer
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

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Check arguments
    check_arguments(global_log, receptor_pdb_path, ligand_pdb_path, output_path, previous_output_path, skip_until)
    
    # Enforce receptor_pdb_path and ligand_pdb_path if provided
    if receptor_pdb_path is not None:
        global_paths["step1_setup"]["input_rec_pdb_path"] = receptor_pdb_path
        global_paths["step2_oda_receptor"]["input_structure_path"] = receptor_pdb_path
    if ligand_pdb_path is not None:
        global_paths["step1_setup"]["input_lig_pdb_path"] = ligand_pdb_path
        global_paths["step3_oda_ligand"]["input_structure_path"] = ligand_pdb_path

    # Add receptor and ligand chains to oda_filtering and distance_filtering step properties
    global_prop["step8_oda_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step8_oda_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    global_prop["step9_distance_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step9_distance_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    # Add docking name and ligand chain to clustering step properties
    global_prop["step10_clustering"]["docking_name"] = global_prop["step7_makePDB"]["docking_name"]
    global_prop["step10_clustering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    
    # Initialize minimal launcher to avoid repeating constant arguments
    launch_step = partial(launch_step_full, global_log, global_paths, output_path, previous_output_path)

    # Skip steps only if requested
    run_remaining_steps = skip_until is None
    setup_time   = launch_step(setup,   "step1_setup", "setup receptor and ligand proteins for pyDock", global_prop, run_remaining_steps)
    oda_time     = launch_step(oda,     "step2_oda_receptor", "optimal docking area (ODA) analysis for the receptor",  global_prop, run_remaining_steps)
    oda_time    += launch_step(oda,     "step3_oda_ligand", "optimal docking area (ODA) analysis for the ligand",  global_prop, run_remaining_steps)
    ftdock_time  = launch_step(ftdock,  "step4_ftdock", "sample docking poses using ftdock (FFT-based algorithm)",  global_prop, run_remaining_steps)
    scoring_time =launch_step(dockser, "step5_dockser", "score docking poses using pyDock",  global_prop, run_remaining_steps)

    # Save the (pyDock score) ranking for later use 
    rank1 = global_prop["step7_makePDB"]["rank1"]
    rank2 = global_prop["step7_makePDB"]["rank2"]
    ranking_df = pd.read_csv(global_paths["step5_dockser"]["output_ene_path"], sep='\s+', skiprows=[1], header=0)
    top_ranking_df = ranking_df[(ranking_df["RANK"] >= rank1) & (ranking_df["RANK"] <= rank2)] 
    global_prop["step10_clustering"]["ranking_df"] = top_ranking_df

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "makePDB"
    makepdb_time = launch_step(makePDB, "step7_makePDB", "generate PDB files for top scoring docking poses",  global_prop, run_remaining_steps)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "oda_filtering"
    oda_filter_time = launch_step(oda_filtering,      "step8_oda_filtering", "filtering with ODA patches",  global_prop, run_remaining_steps)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "distance_filtering"
    dis_filter_time = launch_step(distance_filtering, "step9_distance_filtering", "filtering with distance between residues",  global_prop, run_remaining_steps)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "clustering"
    clu_filter_time = launch_step(rmsd_clustering,    "step10_clustering", "clustering with RMSD",  global_prop, run_remaining_steps)

    # NOTE: Decorate the receptor and the ligand proteins of each pose with the corresponding ODA values
    # oda_poses_folder_path = fu.create_dir(str(Path(properties["path"]).joinpath("oda_poses")))
    # oda_poses_paths = []
    # for pose_path in poses_paths:
         # Decorate the receptor and the ligand proteins of the pose with the corresponding ODA values
    #    decorated_pose_path = decorate_with_oda_values(pose_path, input_receptor_path, input_ligand_path, oda_poses_folder_path, properties["receptor_chain"], properties["ligand_chain"])
    #    oda_poses_paths.append(decorated_pose_path)

    # Print timing information to log file
    total_elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('')
    global_log.info('Total elapsed time: %.2f minutes' % (total_elapsed_time/60))
    global_log.info('')
    global_log.info('Time per step:')
    global_log.info('  setup: %.2f minutes' % (setup_time/60))
    global_log.info('  oda (both): %.2f minutes' % (oda_time/60))
    global_log.info('  ftdock: %.2f minutes' % (ftdock_time/60))
    global_log.info('  dockser: %.2f minutes' % (scoring_time/60))
    global_log.info('  makePDB: %.2f minutes' % (makepdb_time/60))
    global_log.info('  oda_filtering: %.2f minutes' % (oda_filter_time/60))
    global_log.info('  distance_filtering: %.2f minutes' % (dis_filter_time/60))
    global_log.info('  clustering: %.2f minutes' % (clu_filter_time/60))

    return global_paths, global_prop

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Protein-protein docking workflow")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    # Inputs
    parser.add_argument('--receptor_pdb', dest='receptor_pdb_path',
                        help="Input receptor PDB file path",
                        required=False)
    
    parser.add_argument('--ligand_pdb', dest='ligand_pdb_path',
                        help="Input ligand PDB file path (usually the smallest of the two proteins)",
                        required=False)
    
    parser.add_argument('--previous_output', dest='previous_output_path',
                        help="Path to previous output folder. Used if you want to re-use some steps of a previous run (default: None)",
                        required=False)
    
    # To skip all steps until a certain step.
    parser.add_argument('--skip_until', dest='skip_until',
                        help="Skip everything until this step (use together with previous_output to re-use previous results). Options: makePDB, oda_filtering, distance_filtering, clustering",
                        required=False)
    
    # Output                        
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path=args.config_path,
            receptor_pdb_path=args.receptor_pdb_path,
            ligand_pdb_path=args.ligand_pdb_path,
            previous_output_path=args.previous_output_path,
            skip_until=args.skip_until,
            output_path=args.output_path)
