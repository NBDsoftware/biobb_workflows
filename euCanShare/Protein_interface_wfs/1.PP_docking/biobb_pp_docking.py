#!/usr/bin/env python3

# Importing all the needed libraries
import os
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from Bio import PDB

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

def link_previous_steps(global_log, step_folders, output_path, previous_output_path):
    """
    Link step folders in previous output path to output path

    Inputs
    ------

        step_folders          (list): list of step folders to link
        output_path            (str): path to output folder
        previous_output_path   (str): path to previous output folder
    """

    # Check if previous output path is provided
    if previous_output_path is None:
        return

    # For each step folder, check if it exists in previous output path and link it to output path
    for step in step_folders:

        # Check if step folder exists in previous output path
        if not os.path.exists(str(Path(previous_output_path).joinpath(step))):
            global_log.warning(f"WARNING: The step folder {step} does not exist in the previous output path provided!")
            continue

        # Check if step folder exists in output path
        if os.path.exists(str(Path(output_path).joinpath(step))):
            global_log.warning(f"WARNING: The step folder {step} already exists in the output path provided!")
            continue

        # Link step folder in previous output path to output path
        os.symlink(str(Path(previous_output_path).joinpath(step)), str(Path(output_path).joinpath(step)))
    
    return


def rmsd_clustering(input_zip_path: str, output_zip_path: str, properties: dict, ranking_df: pd.DataFrame):

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
        ranking_df     (pd.Dataframe): dataframe with ranking of docking poses
    """

    # Max number of CA atoms to use for RMSD calculation 
    # (enough to describe the rigid ligand orientation)
    max_ca_atoms = 20

    # Prepare clustering step folder 
    poses_ranked_paths, paths_ranking = prepare_clustering_step(input_zip_path, ranking_df, properties)

    # Create a PDB parser 
    parser = PDB.PDBParser(QUIET=True)

    # Get the ligand chain from the properties
    ligand_chain = properties["ligand_chain"]

    # Get the number of poses
    n_poses = len(poses_ranked_paths)

    # Initialize the empty RMSD matrix (n_poses x n_poses)
    rmsd_matrix = np.zeros((n_poses, n_poses))

    # Debug
    global_log = properties["global_log"]
    # Start timer
    start_time = time.time()

    # Get the serials of some CA atoms of the ligand chain from the first pose
    ca_serials = sample_ca_atoms(poses_ranked_paths[0], ligand_chain, max_ca_atoms)

    # Iterate over pairs of poses without repetition
    for i in range(n_poses):

        # Get the ligand chain of the reference pose
        ref_path = poses_ranked_paths[i]
        ref_structure = parser.get_structure("reference", ref_path)
        ligand_ref = find_chain(ref_structure, ligand_chain)

        # Get the selected CA atom coordinates of the reference pose
        ref_atoms = [atom for atom in ligand_ref.get_atoms() if atom.get_serial_number() in ca_serials]
        ref_coords = np.array([atom.get_coord() for atom in ref_atoms])

        for j in range(i + 1, n_poses):
            
            # Get the target chain of the target pose
            target_path = poses_ranked_paths[j]
            target_structure = parser.get_structure("target", target_path)
            ligand_target = find_chain(target_structure, ligand_chain)
            
            # Get the selected CA atom coordinates of the target pose
            target_atoms = [atom for atom in ligand_target.get_atoms() if atom.get_serial_number() in ca_serials]
            target_coords = np.array([atom.get_coord() for atom in target_atoms])

            # Calculate RMSD between reference and target
            rmsd = np.sqrt(np.sum((ref_coords - target_coords) ** 2) / len(ref_coords))
        
            # Add RMSD to the matrix (symmetric)
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    
    # Debug
    elapsed_time = time.time() - start_time
    global_log.info(f"Elapsed time: {elapsed_time / 60} min")

    # Convert the RMSD matrix into a condensed distance matrix
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

def sample_ca_atoms(pdb_path, ligand_chain, max_ca_atoms):
    """
    Sample CA atoms of the ligand chain uniformly.

    Inputs
    ------

        pdb_path       (str): path to pdb file
        ligand_chain   (str): chain of the ligand protein
        max_ca_atoms   (int): maximum number of CA atoms to keep
    
    Returns
    -------

        ca_serial     (list): list with the serial numbers of the sampled CA atoms
    """
    
    # Create a PDB parser 
    parser = PDB.PDBParser(QUIET=True)

    # Get the structure of the first pose
    pose_structure = parser.get_structure("first_pose", pdb_path)

    # Find ligand chain
    ligand_chain = find_chain(pose_structure, ligand_chain)
    
    # Get the CA atoms of the ligand chain
    ligand_ca_atoms = [atom for atom in ligand_chain.get_atoms() if atom.get_id() == "CA"]

    # Get the number of CA atoms of the ligand chain
    total_ca_atoms = len(ligand_ca_atoms)

    # Sample the CA atoms of the ligand chain uniformly
    if total_ca_atoms > max_ca_atoms:

        # Get the indices of the CA atoms to keep
        indices = np.linspace(0, total_ca_atoms - 1, max_ca_atoms, dtype=int)

        # Get the CA atoms to keep
        ligand_ca_atoms = [ligand_ca_atoms[i] for i in indices]

    # Get the serial numbers of the CA atoms to keep
    ca_serial = [atom.get_serial_number() for atom in ligand_ca_atoms]

    return ca_serial

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

def prepare_clustering_step(input_zip_path: dict, ranking_df: pd.DataFrame, prop: dict):
    """ 
    1. Create clustering folder
    2. Create docking poses subfolder and unzip docking poses into it
    3. Rename docking poses with rank and create a list with the corresponding rank for each file

    Inputs
    ------

        input_zip_path          (str): path to zip file with docking poses
        ranking_df     (pd.Dataframe): dataframe with ranking of docking poses
        prop                   (dict): dictionary with clustering step properties

    Returns
    -------

        poses_ranked_paths     (list): list with paths to docking poses
        paths_ranking          (list): list with the corresponding rank for each docking pose
    """
    
    # Create clustering folder
    fu.create_dir(prop["path"])

    # Create subfolder for poses
    poses_folder_path = fu.create_dir(str(Path(prop["path"]).joinpath("poses")))

    # Extract the file with the docking poses into the subfolder
    poses_paths = fu.unzip_list(zip_file = input_zip_path, dest_dir = poses_folder_path)

    # Prepare docking poses - change file names from conformation number to rank and create a list with the corresponding rank for each file
    poses_ranked_paths, paths_ranking = rename_with_rank(poses_paths, ranking_df, prop['docking_name'])

    return poses_ranked_paths, paths_ranking

def rename_with_rank(poses_paths: list, ranking_df: pd.DataFrame, docking_name: str):
    """ 
    Change pdb poses file names from conformation number to rank.

    Inputs
    ------

        poses_paths         (list): list with paths to docking poses
        ranking_df  (pd.DataFrame): dataframe with ranking of docking poses
        docking_name         (str): name of the docking run

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

        # Get the conformation number from the file name
        conformation_number = Path(pdb_path).stem.strip(f"{docking_name}_")

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

    1. Decorate the receptor and ligand in all input docking poses with ODA values.
    2. Find all residues pertaining to the surface of the receptor and ligand proteins (through their SASA values).
    3. Identify interface residues between the receptor and ligand proteins for each pose.
    4. Compute the percentage of interface residues that are covered by ODA patches for each pose. 
    5. Keep only the poses with a percentage of interface residues covered by ODA patches above a given threshold.

    Inputs
    ------

        input_receptor_path     (str): path to receptor pdb file with ODA values
        input_ligand_path       (str): path to ligand pdb file with ODA values
        input_zip_path          (str): path to zip file with docking poses
        output_zip_path         (str): path to output zip file with filtered docking poses
        properties             (dict): dictionary with oda_filtering step properties
    """

    # Decorate the receptor and ligand in all input docking poses with ODA values.
    oda_poses_paths = prepare_oda_filtering_step(input_receptor_path, input_ligand_path, input_zip_path, properties)

    # Find all residues pertaining to the surface of the receptor and ligand proteins (through their SASA values).
    receptor_surface_residues = find_surface_residues(input_receptor_path, properties["path"])
    ligand_surface_residues = find_surface_residues(input_ligand_path, properties["path"])

    # Find the list of atoms defining the receptor surface (receptor remains fixed during docking)
    receptor_surface_atoms = []
    for residue in receptor_surface_residues:
        receptor_surface_atoms.extend(residue.get_atoms())

    # Create the Neighbour Search object for the receptor surface atoms
    Receptor_Neighbor_Search = PDB.NeighborSearch(receptor_surface_atoms)

    # For each pose
    filtered_poses_paths = []
    for pose_path in oda_poses_paths:

        # Identify interface residues for the pose
        interface_residues = find_interface_residues(pose_path, ligand_surface_residues, Receptor_Neighbor_Search, properties)

        # Compute the "PIR COP" (Percentage of Interface Residues Covered by ODA Patches) for the pose
        pir_cop = compute_pir_cop(interface_residues, properties["oda_threshold"])

        # Keep only the poses with a "PIR COP" above a given threshold
        if pir_cop >= properties["pir_cop_threshold"]:
             filtered_poses_paths.append(pose_path)
    
    # Zip the filtered docking pose paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = filtered_poses_paths)

def prepare_oda_filtering_step(input_receptor_path: str, input_ligand_path: str, input_zip_path: str, prop: dict):
    """
    1. Create oda_filtering folder
    2. Create docking poses subfolder and unzip docking poses into it
    3. Decorate the receptor and the ligand proteins of each pose with the corresponding ODA values

    Inputs
    ------

        input_receptor_path     (str): path to receptor pdb file with ODA values
        input_ligand_path       (str): path to ligand pdb file with ODA values
        input_zip_path          (str): path to zip file with docking poses
        prop                   (dict): dictionary with oda_filtering step properties
    
    Returns
    -------

        poses_paths            (list): list with paths to docking poses
    """

    # Create oda filtering folder
    fu.create_dir(prop["path"])

    # Create subfolder for poses
    poses_folder_path = fu.create_dir(str(Path(prop["path"]).joinpath("poses")))

    # Extract the file with the docking poses into the subfolder
    poses_paths = fu.unzip_list(zip_file = input_zip_path, dest_dir = poses_folder_path)

    # Create subfolder for poses with ODA values
    oda_poses_folder_path = fu.create_dir(str(Path(prop["path"]).joinpath("oda_poses")))

    # Iterate over the docking poses
    decorated_poses_paths = []
    for pose_path in poses_paths:

        # Decorate the receptor and the ligand proteins of the pose with the corresponding ODA values
        decorated_pose_path = decorate_with_oda_values(pose_path, input_receptor_path, input_ligand_path, oda_poses_folder_path, prop["receptor_chain"], prop["ligand_chain"])
        decorated_poses_paths.append(decorated_pose_path)

    return decorated_poses_paths

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
    pose_name = pose_path.split("/")[-1]

    # Create new file name for the decorated docking pose
    decorated_pose_name = f"oda_{pose_name}"

    # Create new path for the decorated docking pose
    decorated_pose_path = f"{oda_poses_folder_path}/{decorated_pose_name}"

    # Save the decorated docking pose
    io = PDB.PDBIO()
    io.set_structure(pose_structure)
    io.save(decorated_pose_path)

    return decorated_pose_path

def find_surface_residues(pdb_path: str, output_path: str):
    """
    Find all residues pertaining to the surface of the pdb structure (through the Solvent Accessible Surface Area).
    Save a structure showing the surface residues in the output path.
    Return a list with the surface residues (Bio.PDB.Residue objects)

    Inputs
    ------

        pdb_path (str): path to pdb file

    Returns
    -------

        surface_residues (list): list with surface residues (Bio.PDB.Residue objects)
    """
    # SASA threshold to classify residues as surface or not
    sasa_threshold = 10 

    # Load the pdb structure
    structure = PDB.PDBParser(QUIET=True).get_structure("pdb", pdb_path)

    # Find the pdb name
    pdb_name = Path(pdb_path).name

    # Calculate the SASA for each residue
    sasa_calc = PDB.SASA.ShrakeRupley(probe_radius=1.4, n_points=100)
    sasa_calc.compute(structure, level="R")

    # Classify residues as surface or not
    surface_residues = []
    for residue in structure.get_residues():

        if residue.sasa > sasa_threshold:
            surface_residues.append(residue)
            
            # Mark surface residue
            for atom in residue:
                atom.set_bfactor(100)
        else:
            # Mark non-surface residue
            for atom in residue:
                atom.set_bfactor(0)

    # Save the SASA values in the beta factor column of the pdb file
    io = PDB.PDBIO()
    io.set_structure(structure)
    file_name = f"{pdb_name}_surface.pdb"
    file_path = os.path.join(output_path, file_name)
    io.save(file_path)

    # Check if any surface residue was found
    if len(surface_residues) == 0:
        raise ValueError("No surface residues found!")

    return surface_residues

def find_interface_residues(pose_path: str, ligand_surface_residues: list, Receptor_Neighbor_Search: PDB.NeighborSearch, properties: dict):
    """
    Find the interface residues between the receptor and ligand proteins for a given docking pose.

    Iterates over the ligand surface residues and finds the receptor surface residues that are within a given distance threshold.
    
    Inputs
    ------

        pose_path                               (str): path to docking pose with receptor and ligand
        ligand_surface_residues                (list): list of Bio.PDB.Residue objects, with surface residues of the ligand protein
        Receptor_Neighbor_Search (PDB.NeighborSearch): class with the receptor surface atoms to search for neighbors
        properties                             (dict): dictionary with step properties

    Returns
    -------

        interface_residues (list): list with interface residues of the pose (Bio.PDB.Residue objects)
    """
    # Global log
    global_log = properties["global_log"]

    # Debug
    global_log.info("Finding interface residues...")
    global_log.info(f"Pose: {pose_path}")

    # Load the docking pose using bioPython
    parser = PDB.PDBParser(QUIET=True)
    pose_structure = parser.get_structure("pose", pose_path)

    # Find the ligand chain
    ligand_chain = find_chain(pose_structure, properties["ligand_chain"])

    # Find the ligand surface residues IDs
    ligand_surface_residues_ids = [residue.get_id() for residue in ligand_surface_residues]
    
    receptor_interface_residues = []
    receptor_interface_residues_ids = []
    ligand_interface_residues = []

    # Iterate over all ligand residues
    for ligand_residue in ligand_chain.get_residues():

        # ID of the ligand residue
        ligand_residue_id = ligand_residue.get_id()

        # Check if the ligand residue is a surface residue
        if ligand_residue_id in ligand_surface_residues_ids:

            # Debug 
            global_log.info(f"For ligand surface residue: {ligand_residue}")

            # Get the ligand residue atoms
            ligand_residue_atoms = [atom for atom in ligand_residue.get_atoms()]

            # Get the COM of the ligand residue
            ligand_residue_com = np.mean([atom.get_coord() for atom in ligand_residue_atoms], axis=0)

            # Find the receptor surface residues within a given distance threshold from the ligand surface residue
            receptor_neighbors = Receptor_Neighbor_Search.search(ligand_residue_com, properties["distance_threshold"], level="R")

            # Check if any receptor surface residue was found
            if len(receptor_neighbors) > 0:

                # Debug
                global_log.info(f"Found {len(receptor_neighbors)} receptor surface residues within {properties['distance_threshold']} Å from ligand surface residue {ligand_residue}")
                global_log.info(f"Receptor surface residues: {receptor_neighbors}")

                # Add non-repeated receptor surface residues to the list
                for receptor_residue in receptor_neighbors:
                    if receptor_residue.get_id() not in receptor_interface_residues_ids:
                        receptor_interface_residues.append(receptor_residue)
                        receptor_interface_residues_ids.append(receptor_residue.get_id())

                # Add the ligand surface residue to the interface residues
                ligand_interface_residues.append(ligand_residue)
    
    # Join the receptor and ligand interface residues
    interface_residues = receptor_interface_residues + ligand_interface_residues

    # Check if any interface residue was found
    if len(interface_residues) == 0:
        global_log.error("WARNING: No interface residues were found! Check the distance threshold and the ligand chain!")
    
    return interface_residues

def compute_pir_cop(interface_residues: list, oda_threshold: float):
    """
    Compute the "PIR COP" (Percentage of Interface Residues Covered by ODA Patches) for a given list of
    interface residues given as a list of Bio.PDB.Residue objects. In which the beta factor
    column is the ODA value.

    Inputs
    ------

        interface_residues (list): list with interface residues (Bio.PDB.Residue objects)
        oda_threshold     (float): ODA threshold to consider a residue covered by ODA patches
    
    Returns
    -------

        pir_cop           (float): Percentage of Interface Residues Covered by ODA Patches
    """

    # Get the number of interface residues
    n_interface_residues = len(interface_residues)

    # Iterate over the interface residues and count the number of residues covered by ODA patches
    n_interface_residues_covered = 0
    for residue in interface_residues:

        # Get the atoms of the residue
        atoms = [atom for atom in residue.get_atoms()]

        # Get the ODA value of the residue
        oda_value = atoms[0].bfactor

        # Check if the residue is covered by ODA patches
        if oda_value <= oda_threshold:
            n_interface_residues_covered += 1

    # Compute the "PIR COP" (Percentage of Interface Residues Covered by ODA Patches)
    pir_cop = n_interface_residues_covered / n_interface_residues

    return pir_cop

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

    # Skip steps only if requested
    run_remaining_steps = skip_until is None
    if run_remaining_steps:

        # STEP 1: Prepare receptor and ligand proteins for pyDock
        global_log.info("step1_setup: setup receptor and ligand proteins for pyDock")
        setup(**global_paths["step1_setup"], properties=global_prop["step1_setup"])
    
        # STEP 2: Optimal Docking Area (ODA) analysis for the receptor
        global_log.info("step2_oda_receptor: optimal docking area (ODA) analysis for the receptor")
        oda(**global_paths["step2_oda_receptor"], properties=global_prop["step2_oda_receptor"])
    
        # STEP 3: Optimal Docking Area (ODA) analysis for the ligand
        global_log.info("step3_oda_ligand: optimal docking area (ODA) analysis for the ligand")
        oda(**global_paths["step3_oda_ligand"], properties=global_prop["step3_oda_ligand"])

        # STEP 4: Sample docking poses using ftdock (FFT-based algorithm)
        global_log.info("step4_ftdock: sample docking poses using ftdock (FFT-based algorithm)")
        ftdock(**global_paths["step4_ftdock"], properties=global_prop["step4_ftdock"])

        # STEP 5: Score docking poses using pyDock
        global_log.info("step5_dockser: score docking poses using pyDock")
        dockser(**global_paths["step5_dockser"], properties=global_prop["step5_dockser"])

    else:

        # Link Steps 1-5 folders in previous output path to output path
        steps_to_link = ["step1_setup", "step2_oda_receptor", "step3_oda_ligand", "step4_ftdock", "step5_dockser"]
        link_previous_steps(global_log, steps_to_link, output_path, previous_output_path)

    # Read csv with ranking of docking poses (pydock score)
    ranking_df = pd.read_csv(global_paths["step5_dockser"]["output_ene_path"], sep='\s+', skiprows=[1], header=0)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "makePDB"
    if run_remaining_steps:

        # STEP 7: Generate PDB files for top scoring docking poses
        global_log.info("step7_makePDB: generate PDB files for top scoring docking poses")
        makePDB(**global_paths["step7_makePDB"], properties=global_prop["step7_makePDB"])
    
    else:

        # Link step folder in previous output path to output path
        link_previous_steps(global_log, ["step7_makePDB"], output_path, previous_output_path)

    # Keep only top scoring docking poses in the ranking data frame
    rank1 = global_prop["step7_makePDB"]["rank1"]
    rank2 = global_prop["step7_makePDB"]["rank2"]
    top_ranking_df = ranking_df[(ranking_df["RANK"] >= rank1) & (ranking_df["RANK"] <= rank2)] 

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "clustering"
    if run_remaining_steps:
       
        # Add docking name and ligand chain to clustering step properties
        global_prop["step8_clustering"]["docking_name"] = global_prop["step7_makePDB"]["docking_name"]
        global_prop["step8_clustering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]

        # STEP 8: Clustering with RMSD
        global_log.info("step8_clustering: clustering with RMSD")
        rmsd_clustering(**global_paths["step8_clustering"], properties=global_prop["step8_clustering"], ranking_df=top_ranking_df) 

    else:

        # Link step folder in previous output path to output path
        link_previous_steps(global_log, ["step8_clustering"], output_path, previous_output_path)

    # Add receptor and ligand chains to oda_filtering step properties
    global_prop["step9_oda_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step9_oda_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]

    # STEP 9: Filtering with ODA patches
    global_log.info("step9_oda_filtering: filtering with ODA patches")
    oda_filtering(**global_paths["step9_oda_filtering"], properties=global_prop["step9_oda_filtering"])

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
    
    # To skip all steps until a certain step. Options: makePDB, clustering, oda_filtering
    parser.add_argument('--skip_until', dest='skip_until',
                        help="Skip everything until this step (use together with previous_output to re-use previous results). Options: makePDB, clustering, oda_filtering",
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
