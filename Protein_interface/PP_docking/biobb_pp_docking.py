#!/shared/work/BiobbWorkflows/envs/biobb_pp_docking/bin/python

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
        num_filtered_poses = tool(**global_paths[step_name], properties=global_prop[step_name])

        if 'filtering' in step_name:
            global_log.info(f"  Number of filtered poses: {num_filtered_poses}")
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

    # If skip_until is provided, check if it is a valid option
    if skip_until is not None:
        if skip_until not in ["makePDB", "oda_filtering", "distance_filtering", "clustering"]:
            global_log.error("ERROR: --skip_until must be one of the following options: makePDB, oda_filtering, distance_filtering, clustering!")

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
    Prepares a given step by:
    
    1. Creating a step folder
    2. Creating a docking poses subfolder and unzipping input docking poses into it

    Inputs
    ------

        input_zip_path          (str): path to zip file with input docking poses
        prop                   (dict): dictionary with step properties

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
    
    Outputs
    -------

        num_filtered_poses        (int): number of filtered poses
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return 0
    
    # Prepare clustering step folder 
    poses_paths = prepare_step(input_zip_path, properties)

    # Load the poses table
    poses_table = properties["poses_table"]

    # Rename docking poses with rank and create a list with the corresponding rank for each file
    poses_ranking = get_poses_ranking(poses_paths, poses_table)

    # Calculate the RMSD matrix between all poses
    rmsd_matrix = compute_rmsd_matrix(poses_paths, properties["ligand_chain"])

    # Convert the uncondensed RMSD matrix into a condensed distance matrix
    rmsd_matrix = squareform(rmsd_matrix)

    # Calculate the linkage matrix
    Z = linkage(rmsd_matrix, method = "average")

    # Get the cluster labels from the linkage matrix
    cluster_labels = fcluster(Z, t = properties["rmsd_threshold"], criterion = "distance")

    plot_dendogram(Z, poses_ranking, properties["rmsd_threshold"], properties["path"])

    # Find the best ranking pose path for each cluster
    filtered_poses_paths = []
    for cluster in range(1, max(cluster_labels) + 1):

        # Get the ranking of the poses in the cluster
        cluster_ranking = [poses_ranking[i] for i in np.where(cluster_labels == cluster)[0]]

        # Get the paths of the poses in the cluster
        cluster_paths = [poses_paths[i] for i in np.where(cluster_labels == cluster)[0]]

        if properties["keep_all"]:
            # Add all poses in the cluster
            filtered_poses_paths.extend(cluster_paths)
            
            # For each pose, add the cluster label
            for pose_path in cluster_paths:

                # Get the conformation number from the pdb name
                conformation_number = get_conf_num(pose_path)

                # Add cluster label to this pose in the poses table
                poses_table.loc[poses_table["Conf"] == int(conformation_number), "cluster"] = cluster
        else:
            # Add the pose with the best ranking in the cluster
            best_pose_path = cluster_paths[cluster_ranking.index(min(cluster_ranking))]
            filtered_poses_paths.append(best_pose_path)

            # Get the conformation number from the pdb name
            conformation_number = get_conf_num(best_pose_path)

            # Add cluster label to this pose in the poses table
            poses_table.loc[poses_table["Conf"] == int(conformation_number), "cluster"] = cluster

    # Zip the best ranking pose paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = filtered_poses_paths)

    # Remove temporal folder with poses
    fu.rm(os.path.join(properties["path"], "poses"))

    return len(poses_paths)-len(filtered_poses_paths)

def compute_rmsd_matrix(poses_paths: list, ligand_chain: str):
    """
    Computes an uncondensed RMSD matrix between all poses using only the CA atoms of the ligand (as the receptor didn't move during docking).

    Inputs
    ------

        poses_paths        (list): list with paths to docking poses
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
    poses_universe = mda.Universe(poses_paths[0])

    # Read coordinates of all poses into the universe
    poses_universe.load_new(poses_paths)

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

def get_poses_ranking(poses_paths: list, poses_table: pd.DataFrame):
    """ 
    Get ranking of the poses in the input paths from the ranking dataframe

    Inputs
    ------

        poses_paths         (list): list with paths to docking poses
        poses_table (pd.DataFrame): dataframe with ranking of docking poses

    Returns
    -------

        poses_ranking      (list): list with the corresponding rank for each docking pose
    
    """
    # Create a list with the corresponding rank for each docking pose
    poses_ranking = []

    # Iterate over all poses
    for pdb_path in poses_paths:

        # Get the conformation number from the pdb name
        conformation_number = get_conf_num(pdb_path)

        # Get the rank from the ranking Dataframe using the conformation number
        rank = poses_table[poses_table["Conf"] == int(conformation_number)]["RANK"].values[0]
        poses_ranking.append(int(rank))
    
    return poses_ranking

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

def get_conf_num(pdb_path: str):
    """
    Get the conformation number from the pdb name.

    Inputs
    ------

        pdb_path (str): path to pdb file
    
    Returns
    -------

        conformation_number (str): conformation number
    """

    # Get the pdb name without the extension
    pdb_name = Path(pdb_path).stem

    # Look for a number (e.g. 13 or 123) in the pdb name    
    number_match = re.findall(r"\d+", pdb_name)

    # Get the conformation number from the match
    conformation_number = number_match[0]

    return conformation_number


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
    
    Outputs
    -------

        num_filtered_poses        (int): number of filtered poses
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return 0
    
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

        # Get the conformation number from the pdb name
        conformation_number = get_conf_num(pose_path)

        # Identify interface residues for the pose
        receptor_interface_atoms, ligand_interface_atoms = find_interface_residues(pose_path, ligand_surface_oda.keys(), Receptor_Neighbor_Search, properties)

        # Filter the pose 
        accept = filter_pose_interface(conformation_number, receptor_interface_atoms, ligand_interface_atoms, 
                                       receptor_surface_oda, ligand_surface_oda, properties)

        if accept:
            filtered_poses_paths.append(pose_path)
    
    # Zip the filtered docking pose paths into a zip file
    fu.zip_list(zip_file = output_zip_path, file_list = filtered_poses_paths)

    # Remove temporal folders with poses and decorated poses
    fu.rm(os.path.join(properties["path"], "poses"))

    return len(poses_paths)-len(filtered_poses_paths)

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

def filter_pose_interface(conformation_number, receptor_interface_atoms, ligand_interface_atoms, receptor_surface_oda, ligand_surface_oda, properties: dict):
    """
    Filter the pose according to the ODA patches and the specific interface. Two main criteria:

        1. ODA coverage above threshold (compulsory)
            
            Compute the ODA coverage (ratio of interface residues covered by ODA patches) for a given pose.
            This measures the overlap between the interface and the ODA patches. If the overlap is above a given threshold,
            the pose is accepted.
        
        2. Percentage of neighboring interface residues covered by oda patches between receptor and ligand.
        
            This measures overlap of ligand and receptor oda patches inside the interface. If the overlap is above a given threshold,
            the pose is accepted. This is optional and can be turned off by setting the threshold to 0 for faster filtering.

    It will also update the poses table in the properties dictionary with the ODA coverage and overlap values for this pose.

    Inputs
    ------

        conformation_number       (int): conformation number of the pose
        receptor_interface_atoms (dict): dictionary with residue IDs of receptor interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
        ligand_interface_atoms   (dict): dictionary with residue IDs of ligand interface residues as keys and their CA atoms (Bio.PDB.Atom.Atom objects) as values
        receptor_surface_oda     (dict): dictionary with residue IDs of receptor surface residues as keys and their ODA values as values
        ligand_surface_oda       (dict): dictionary with residue IDs of ligand surface residues as keys and their ODA values as values
        properties               (dict): dictionary with step properties
    
    Returns
    -------

        accept                   (bool): True if the pose is accepted, False otherwise
    """

    # Get the poses table 
    poses_table = properties["poses_table"]

    # Get the total number of interface residues
    num_interface_residues = len(receptor_interface_atoms) + len(ligand_interface_atoms)

    # Count the number of residues in the interface covered by ODA patches
    num_covered_residues = 0

    # Count the number of residues in the interface covered by ODA patches with at least one neighboring residue (from the other protein) covered by an ODA patch as well 
    num_overlap_residues = 0

    # Create Neighbor Search objects for the interface surface atoms 
    if properties["oda_overlap"] > 0:
        Receptor_Neighbor_Search = PDB.NeighborSearch(list(receptor_interface_atoms.values())) # NOTE: optimize bucket size?
        Ligand_Neighbor_Search = PDB.NeighborSearch(list(ligand_interface_atoms.values())) # NOTE: optimize bucket size?

    # Find covered and overlapping residues in the receptor interface
    for residue_id, residue_ca_atom in receptor_interface_atoms.items():

        # Get the ODA value of the residue
        residue_oda = receptor_surface_oda[residue_id]

        # Check if the residue is covered by an ODA patch
        if residue_oda >= properties["oda_threshold"]:
            num_covered_residues += 1

            if properties["oda_overlap"] > 0:

                # Find the neighboring ligand interface residues
                ligand_neighbors = Ligand_Neighbor_Search.search(residue_ca_atom.get_coord(), properties["distance_threshold"], level="R")

                # If any neighboring residues are found
                if len(ligand_neighbors) > 0:

                    # Check if any of them are covered by an ODA patch as well
                    for ligand_residue in ligand_neighbors:
                        if ligand_surface_oda[ligand_residue.get_id()] >= properties["oda_threshold"]:
                            num_overlap_residues += 1
                            break

    # Find covered and overlapping residues in the ligand interface
    for residue_id in ligand_interface_atoms.keys():

        # Get the ODA value of the residue
        residue_oda = ligand_surface_oda[residue_id]

        # Check if the residue is covered by an ODA patch
        if residue_oda >= properties["oda_threshold"]:
            num_covered_residues += 1

            if properties["oda_overlap"] > 0:

                # Find the neighboring receptor interface residues
                receptor_neighbors = Receptor_Neighbor_Search.search(residue_ca_atom.get_coord(), properties["distance_threshold"], level="R")

                # If any neighboring residues are found
                if len(receptor_neighbors) > 0:

                    # Check if any of them are covered by an ODA patch as well
                    for receptor_residue in receptor_neighbors:
                        if receptor_surface_oda[receptor_residue.get_id()] >= properties["oda_threshold"]:
                            num_overlap_residues += 1
                            break

    # Compute the oda coverage: ratio of interface residues covered by ODA patches
    oda_coverage = num_covered_residues / num_interface_residues

    # Update "oda_coverage" column of this conformation number in the poses table
    poses_table.loc[poses_table["Conf"] == int(conformation_number), "oda_coverage"] = round(oda_coverage,3)

    # Check if the oda coverage is above the threshold
    if oda_coverage < properties["oda_coverage"]:
        return False
    
    # If the overlap threshold is 0, return True
    if properties["oda_overlap"] == 0:
        return True
    
    # Compute the oda overlap: ratio of interface residues covered by ODA patches with at least one neighboring residue (from the other protein) covered by an ODA patch as well
    oda_overlap = num_overlap_residues / num_interface_residues

    # Update "oda_overlap" column of this conformation number in the poses table
    poses_table.loc[poses_table["Conf"] == int(conformation_number), "oda_overlap"] = round(oda_overlap,3)

    # Check if the oda overlap is above the threshold
    if oda_overlap < properties["oda_overlap"]:
        return False
    
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

    Outputs
    -------

        num_filtered_poses         (int): number of filtered poses
    """

    if properties["run_step"] is False:

        # Create step folder
        fu.create_dir(properties["path"])
        
        # Copy input zip file to output zip file
        shutil.copy(input_zip_path, output_zip_path)

        return 0

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

    return len(poses_paths)-len(filtered_poses_paths)

def calculate_distances(pose_path: str, properties: dict):
    """
    Compute distances between pairs of residues defined in properties and check
    they are below a given threshold. If all distances are below the threshold,
    return True. Otherwise return False.

    It also saves the value of the distance constraints in the poses table in the properties dictionary.

    Inputs
    ------

        pose_path (str): path to docking pose
        properties (dict): dictionary with step properties
    
    Return
    ------

        accept_pose (bool): True if all distances are below the threshold, False otherwise
    """

    # Get the conformation number from the pdb name
    conformation_number = get_conf_num(pose_path)

    # Get the poses table
    poses_table = properties["poses_table"]

    # Read the docking pose
    pose_universe = mda.Universe(pose_path)

    # Get the receptor and ligand chains 
    receptor_chain = properties["receptor_chain"]
    ligand_chain = properties["ligand_chain"]

    # For each distance in the properties
    for distance in properties["distances"]:

        # Get the receptor residue CA atom
        receptor_residue = pose_universe.select_atoms(f"protein and chainID {receptor_chain} and ({distance['receptor_residue_selection']}) and name CA")

        # Get the ligand residue CA atom
        ligand_residue = pose_universe.select_atoms(f"protein and chainID {ligand_chain} and ({distance['ligand_residue_selection']}) and name CA")

        # Get the distance between the receptor and ligand CA atoms
        distance_value = distance_array(receptor_residue.positions, ligand_residue.positions)[0][0]

        # Update the poses table with the distance value
        poses_table.loc[poses_table["Conf"] == int(conformation_number), distance["name"]] = round(distance_value, 3)

        # Check if the distance is below the threshold
        if distance_value > distance["threshold"]:
            return False
    
    return True


def decorate_with_oda(input_zip_path: str, input_receptor_path: str, input_ligand_path: str, output_zip_path: str, properties: dict):
    """
    Decorate all input docking poses with the ODA values of the receptor and ligand proteins and save them into a zip file. The decoration is
    done writing the ODA values in the B-factor column of each pdb file.

    Inputs
    ------

        input_zip_path          (str): path to input zip file with docking poses
        input_receptor_path     (str): path to receptor pdb file with ODA values
        input_ligand_path       (str): path to ligand pdb file with ODA values
        output_zip_path         (str): path to output zip file with decorated docking poses
        properties             (dict): dictionary with decorate_with_oda step properties
    """

    # Prepare step folder and unzip docking poses into it
    poses_paths = prepare_step(input_zip_path, properties)

    # Load the poses table
    poses_table = properties["poses_table"]

    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Load the receptor and ligand with the ODA values
    receptor_structure = parser.get_structure("receptor_oda", input_receptor_path)
    ligand_structure = parser.get_structure("ligand_oda", input_ligand_path)

    receptor_atoms_oda = []
    ligand_atoms_oda = []

    # Get the atoms with ODA values from the receptor and ligand structures
    for model in receptor_structure:
        for chain in model:
            if chain.id == properties["receptor_chain"]:
                for atom in chain.get_atoms():
                    if atom.bfactor is not None:
                        receptor_atoms_oda.append(atom)
    for model in ligand_structure:
        for chain in model:
            if chain.id == properties["ligand_chain"]:
                for atom in chain.get_atoms():
                    if atom.bfactor is not None:
                        ligand_atoms_oda.append(atom)

    # Iterate over all poses
    for pose_path in poses_paths:

        # Get the conformation number from the pdb name
        conformation_number = get_conf_num(pose_path)

        # Decorate the pose with the ODA values
        decorated_pose_path = decorate_pose_with_oda(pose_path, receptor_atoms_oda, ligand_atoms_oda, properties)

        # Update the poses table with the decorated pose path
        poses_table.loc[poses_table["Conf"] == int(conformation_number), "path"] = decorated_pose_path

    # Remove temporal folders with poses
    fu.rm(os.path.join(properties["path"], "poses"))

def decorate_pose_with_oda(pose_path, receptor_atoms_oda, ligand_atoms_oda, properties):
    """
    Takes a docking pose and decorates the receptor and ligand proteins with the corresponding ODA values.

    Inputs
    ------

        pose_path               (str): path to docking pose (contains receptor and ligand)
        receptor_atoms_oda     (list): list with Bio.PDB.Atom.Atom objects with the ODA values of the receptor protein
        ligand_atoms_oda       (list): list with Bio.PDB.Atom.Atom objects with the ODA values of the ligand protein
        properties             (dict): dictionary with decorate_with_oda step properties
    
    Outputs
    -------

        decorated_pose_path     (str): path to decorated docking pose
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
            if chain.id == properties["receptor_chain"]:
                receptor_atoms.extend(chain.get_atoms())
            elif chain.id == properties["ligand_chain"]:
                ligand_atoms.extend(chain.get_atoms())

    # Iterate over all ligand atoms and decorate them with the corresponding ODA value
    for ligand_atom, ligand_oda_atom in zip(ligand_atoms, ligand_atoms_oda):
        ligand_atom.set_bfactor(ligand_oda_atom.bfactor)

    # Iterate over all receptor atoms and decorate them with the corresponding ODA value
    for receptor_atom, receptor_oda_atom in zip(receptor_atoms, receptor_atoms_oda):
        receptor_atom.set_bfactor(receptor_oda_atom.bfactor)

    # Create new path for the decorated docking pose
    pose_name = Path(pose_path).name
    decorated_pose_path = os.path.join(properties["path"], pose_name)

    # Save the decorated docking pose
    io = PDB.PDBIO()
    io.set_structure(pose_structure)
    io.save(decorated_pose_path)

    return decorated_pose_path


def create_poses_table(global_prop: dict, global_paths: dict):
    """
    Create a dataframe with the ranking of the docking poses and their pydock energy score. Add the
    oda_coverage, oda_overlap, cluster and distance constraint columns with a default value of None. 
    Add the dataframe to the properties of the next step.

    Inputs
    ------

        global_prop (dict): dictionary with global properties
        global_paths (dict): dictionary with global paths
    """

    # Find the range of poses that have been extracted as pdb files
    rank1 = global_prop["step6_makePDB"]["rank1"]
    rank2 = global_prop["step6_makePDB"]["rank2"]

    # Read the whole .ene file with the pydock energy scores and ranking
    poses_table = pd.read_csv(global_paths["step5_dockser"]["output_ene_path"], sep='\s+', skiprows=[1], header=0)
    
    # Keep only the rows with the poses that have been extracted as pdb files
    top_poses_table = poses_table[(poses_table["RANK"] >= rank1) & (poses_table["RANK"] <= rank2)]
    
    # Add the oda_coverage, oda_overlap and cluster columns with None
    top_poses_table["oda_coverage"] = None
    top_poses_table["oda_overlap"] = None
    top_poses_table["cluster"] = None 

    # If there is any distance constraint defined
    if global_prop["step8_distance_filtering"].get("distances") is not None:

        # Add the distance constraint columns with None
        for distance in global_prop["step8_distance_filtering"]["distances"]:
            top_poses_table[distance["name"]] = None
    
    # Add the path column with the path to the pdb file (empty for now)
    top_poses_table["path"] = ""
    
    # Add the dataframe to the properties of the next step
    global_prop["step7_oda_filtering"]["poses_table"] = top_poses_table

def clean_poses_table(poses_table):
    """
    Clean the poses table:
        
        - Remove any column with only None values. Corresponds to filtering steps that have not been run.
        - Remove any row that contains a None value. Corresponds to poses that have been filtered and didn't reach that step.
        - Change some column names for clarity.

    Inputs
    ------

        poses_table (pd.DataFrame): poses table with ranking and pydock energy scores
    
    Returns
    -------

        poses_table (pd.DataFrame): cleaned poses table
    """

    # Columns to rename
    original_names = ['Conf',         'Ele',           'Desolv',      'VDW',           'Total', 'RANK', 'oda_coverage', 'oda_overlap', 'cluster']
    new_names =      ['Conformation', 'Electrostatic', 'Desolvation', 'Van Der Waals', 'Total', 'Rank', 'ODA coverage', 'ODA overlap', 'Cluster']

    # Remove any column that only contains None values
    poses_table = poses_table.dropna(axis=1, how='all')

    # Remove any row that contains a None value
    poses_table = poses_table.dropna(axis=0, how='any')

    # Rename some columns for clarity
    poses_table = poses_table.rename(columns=dict(zip(original_names, new_names)))

    return poses_table


# YML construction
def config_contents() -> str:
    """
    Returns the contents of the YAML configuration file as a string.
    
    The YAML file contains the configuration for the protein preparation workflow.
    
    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    return f""" 
working_dir_path: output      # Folder to write i/o files of the workflow steps
can_write_console_log: False  # Verbose writing of log information
restart: True                 # Skip steps already performed
remove_tmp: True

step1_setup:
  tool: setup
  paths:
    input_rec_pdb_path: /path/to/receptor.pdb
    input_lig_pdb_path: /path/to/ligand.pdb
    output_rec_path: prepared_receptor.pdb
    output_rec_H_path: prepared_receptor.pdb.H
    output_rec_amber_path: prepared_receptor.pdb.amber
    output_lig_path: prepared_ligand.pdb
    output_lig_H_path: prepared_ligand.pdb.H
    output_lig_amber_path: prepared_ligand.pdb.amber
    output_ref_path: prepared_reference.pdb
  properties:
    docking_name: docking_test
    receptor: 
      mol: A
      newmol: A
    ligand: 
      mol: A
      newmol: B
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec

step2_oda_receptor:
  tool: oda 
  paths:
    input_structure_path: dependency/step1_setup/input_rec_pdb_path
    output_oda_path: receptor.pdb.oda
    output_oda_H_path: receptor.pdb.oda.H
    output_oda_tab_path: receptor.pdb.oda.ODAtab
    output_oda_amber_path: receptor.oda.amber
  properties:
    subunit_name: receptor
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec

step3_oda_ligand:
  tool: oda 
  paths:
    input_structure_path: dependency/step1_setup/input_lig_pdb_path
    output_oda_path: ligand.pdb.oda
    output_oda_H_path: ligand.pdb.oda.H
    output_oda_tab_path: ligand.pdb.oda.ODAtab
    output_oda_amber_path: ligand.oda.amber
  properties:
    subunit_name: ligand
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec
    
step4_ftdock:
  tool: ftdock
  paths:
    input_rec_path: dependency/step1_setup/output_rec_path
    input_lig_path: dependency/step1_setup/output_lig_path
    output_ftdock_path: ftdock_output.ftdock
    output_rot_path: rotftdock_output.rot
  properties:
    docking_name: docking_test
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec

step5_dockser:
  tool: dockser
  paths:
    input_rec_path: dependency/step1_setup/output_rec_path
    input_rec_H_path: dependency/step1_setup/output_rec_H_path
    input_rec_amber_path: dependency/step1_setup/output_rec_amber_path
    input_lig_path: dependency/step1_setup/output_lig_path
    input_lig_H_path: dependency/step1_setup/output_lig_H_path
    input_lig_amber_path: dependency/step1_setup/output_lig_amber_path
    input_rot_path: dependency/step4_ftdock/output_rot_path
    output_ene_path: dockser_output.ene
  properties:
    docking_name: docking_test
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec

step6_makePDB:
  tool: makePDB
  paths:
    input_rec_path: dependency/step1_setup/output_rec_path
    input_rec_H_path: dependency/step1_setup/output_rec_H_path
    input_rec_amber_path: dependency/step1_setup/output_rec_amber_path
    input_lig_path: dependency/step1_setup/output_lig_path
    input_lig_H_path: dependency/step1_setup/output_lig_H_path
    input_lig_amber_path: dependency/step1_setup/output_lig_amber_path
    input_rot_path: dependency/step4_ftdock/output_rot_path
    input_ene_path: dependency/step5_dockser/output_ene_path
    output_zip_path: top_poses.zip
  properties:
    docking_name: docking_test
    rank1: 1
    rank2: 100
    container_path: singularity
    container_image: /shared/work/NBD_Utilities/pyDock3/3.6.1/nbd_pydock.sif
    container_volume_path: /data
    container_working_dir: /
    container_generic_command: exec

step7_oda_filtering:
  paths: 
    input_receptor_path: dependency/step2_oda_receptor/output_oda_path
    input_ligand_path: dependency/step3_oda_ligand/output_oda_path
    input_zip_path: dependency/step6_makePDB/output_zip_path
    output_zip_path: top_filtered_poses.zip
  properties:
    oda_threshold: -2                     # ODA score threshold for a residue to be considered as part of a hydrophobic patch (oda patch)
    distance_threshold: 8.0               # Distance threshold in A for two residues to be considered as neighbors
    oda_coverage: 0.1                     # Minimum ratio of interface residues covered by ODA Patches -> measures overlap between interface and oda patches
    oda_overlap: 0.05                     # Minimum ratio of neighboring interface residues covered by oda patches (between receptor and ligand) -> measures overlap of ligand and receptor oda patches inside the interface (0 to ignore) 
    run_step: True                        # Run this step

step8_distance_filtering:
  paths: 
    input_zip_path: dependency/step7_oda_filtering/output_zip_path
    output_zip_path: top_filtered_poses.zip
  properties:
    distances:
      - name: "dummy_constraint"
        receptor_residue_selection: "resnum 336 and resname GLN"
        ligand_residue_selection: "resnum 93 and resname THR"
        threshold: 12.0
    run_step: True

step9_rmsd_filtering:
  paths: 
    input_zip_path: dependency/step8_distance_filtering/output_zip_path
    output_zip_path: top_filtered_poses.zip
  properties:
    rmsd_threshold: 5.0          # RMSD threshold for the hierarchical clustering
    keep_all: True               # Keep all poses in the output zip file or just the best ranking one from each cluster
    run_step: True

step10_oda_decoration:
  paths:
    input_zip_path: dependency/step9_rmsd_filtering/output_zip_path
    input_receptor_path: dependency/step2_oda_receptor/output_oda_path
    input_ligand_path: dependency/step3_oda_ligand/output_oda_path
    output_zip_path: top_filtered_poses.zip
"""

def create_config_file(config_path: str) -> None:
    """
    Create a YAML configuration file for the workflow if needed.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file to be created.
    
    Returns
    -------
    None
    """
    
    # Check if the file already exists
    if os.path.exists(config_path):
        print(f"Configuration file already exists at {config_path}.")
        return
    
    # Write the contents to the file
    with open(config_path, 'w') as f:
        f.write(config_contents())
        
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

    if configuration_path is None:
        # Create a default configuration file
        configuration_path = "config.yml"
        create_config_file(configuration_path)
        
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
    global_prop["step7_oda_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step7_oda_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    global_prop["step8_distance_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step8_distance_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    # Add docking name and ligand chain to clustering step properties
    global_prop["step9_rmsd_filtering"]["docking_name"] = global_prop["step6_makePDB"]["docking_name"]
    global_prop["step9_rmsd_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]
    
    # Initialize minimal launcher to avoid repeating constant arguments
    launch_step = partial(launch_step_full, global_log, global_paths, output_path, previous_output_path)

    # Skip steps only if requested
    run_remaining_steps = skip_until is None
    setup_time   = launch_step(setup,   "step1_setup", "setup receptor and ligand proteins for pyDock", global_prop, run_remaining_steps)
    oda_time     = launch_step(oda,     "step2_oda_receptor", "optimal docking area (ODA) analysis for the receptor",  global_prop, run_remaining_steps)
    oda_time    += launch_step(oda,     "step3_oda_ligand", "optimal docking area (ODA) analysis for the ligand",  global_prop, run_remaining_steps)
    ftdock_time  = launch_step(ftdock,  "step4_ftdock", "sample docking poses using ftdock (FFT-based algorithm)",  global_prop, run_remaining_steps)
    scoring_time = launch_step(dockser, "step5_dockser", "score docking poses using pyDock",  global_prop, run_remaining_steps)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "makePDB"
    makepdb_time = launch_step(makePDB, "step6_makePDB", "generate PDB files for top scoring docking poses",  global_prop, run_remaining_steps)

    # Read the pydock score ranking and create poses table
    create_poses_table(global_prop, global_paths)

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "oda_filtering"
    oda_filter_time = launch_step(oda_filtering,      "step7_oda_filtering", "filtering with ODA patches",  global_prop, run_remaining_steps)

    # Pass the poses table to the next step
    global_prop["step8_distance_filtering"]["poses_table"] = global_prop["step7_oda_filtering"]["poses_table"]

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "distance_filtering"
    dis_filter_time = launch_step(distance_filtering, "step8_distance_filtering", "filtering with distance between residues",  global_prop, run_remaining_steps)

    # Pass the poses table to the next step
    global_prop["step9_rmsd_filtering"]["poses_table"] = global_prop["step8_distance_filtering"]["poses_table"]

    # Run step or link previous run
    run_remaining_steps = run_remaining_steps or skip_until == "clustering"
    clu_filter_time = launch_step(rmsd_clustering,    "step9_rmsd_filtering", "clustering with RMSD",  global_prop, run_remaining_steps)

    # Pass the poses table to the next step
    global_prop["step10_oda_decoration"]["poses_table"] = global_prop["step9_rmsd_filtering"]["poses_table"]

    # Pass the input ligand and receptor chains to next step
    global_prop["step10_oda_decoration"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["mol"]
    global_prop["step10_oda_decoration"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["mol"]

    # Decorate the receptor and the ligand proteins of each pose with the corresponding ODA values
    decoration_time = launch_step(decorate_with_oda, "step10_oda_decoration", "decorate receptor and ligand proteins of each pose with the corresponding ODA values", global_prop, True)

    # Save the poses table to a csv file
    final_poses_table = clean_poses_table(global_prop["step10_oda_decoration"]["poses_table"])
    final_poses_table.to_csv(os.path.join(output_path, "summary.csv"), index=False)

    # Print timing information to log file
    total_elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow name: Protein protein docking')
    global_log.info('  Output path: %s' % output_path)
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
    global_log.info('  rmsd_filtering: %.2f minutes' % (clu_filter_time/60))
    global_log.info('  oda_decoration: %.2f minutes' % (decoration_time/60))

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
