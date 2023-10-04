#!/usr/bin/env python3

# Importing all the needed libraries
import os
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_pydock.pydock.setup import setup
from biobb_pydock.pydock.ftdock import ftdock
from biobb_pydock.pydock.dockser import dockser
from biobb_pydock.pydock.makePDB import makePDB
from biobb_pydock.pydock.oda import oda

import MDAnalysis as mda
from MDAnalysis.analysis.diffusionmap import DistanceMatrix
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

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

    The input zip file contains docking poses in PDB format. The function will decorate the receptor and ligand with the ODA values.
    Then a threshold for the ODA values will be used to determine if a residue is part of an ODA patch or not. Finally the docking
    poses will be filtered based on how many interface neighbor residues pertain to an ODA patch.
    """

    # Prepare oda_filtering step folder
    oda_poses_paths = prepare_oda_filtering_step(input_receptor_path, input_ligand_path, input_zip_path, properties)

    # Filter docking poses based on ODA patches 

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

def decorate_with_oda_values(pose_path: str, input_receptor_path: str, input_ligand_path: str, oda_poses_folder_path: str, receptor_chain: str, ligand_chain: str):
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

    # Load the docking pose
    pose = mda.Universe(pose_path, format = "PDB", topology_format="PDB")

    # Select the receptor and ligand atoms from the docking pose
    receptor = pose.select_atoms(f"chainid {receptor_chain}")

    # Show the receptor atoms
    print(receptor)

    ligand = pose.select_atoms(f"chainid {ligand_chain}")

    # Load the receptor and ligand with the ODA values 
    receptor_oda = mda.Universe(input_receptor_path, format = "PDB", topology_format="PDB") 
    ligand_oda = mda.Universe(input_ligand_path, format = "PDB", topology_format="PDB")

    # Iterate over all ligand atoms and decorate them with the corresponding ODA value
    for atom, atom_oda in zip(ligand.atoms, ligand_oda.atoms):
        atom.bfactor = atom_oda.bfactor
    
    # Iterate over all receptor atoms and decorate them with the corresponding ODA value
    for atom, atom_oda in zip(receptor.atoms, receptor_oda.atoms):
        atom.bfactor = atom_oda.bfactor

    # Find the name of the docking pose file
    pose_name = Path(pose_path).name

    # Create new file name for the decorated docking pose
    decorated_pose_name = f"oda_{pose_name}"

    # Create new path for the decorated docking pose
    decorated_pose_path = str(Path(oda_poses_folder_path).joinpath(decorated_pose_name))

    # Save the decorated docking pose
    pose.atoms.write(decorated_pose_path)

    return decorated_pose_path

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

    # Prepare clustering step folder 
    poses_ranked_paths, paths_ranking = prepare_clustering_step(input_zip_path, ranking_df, properties)

    # - it could be parallelized with MDAnalysis.analysis.encore.confdistmatrix if needed
    # Calculate the RMSD matrix using only CA ligand atoms without alignment 
    u = mda.Universe(poses_ranked_paths[0], poses_ranked_paths)
    ligand_chain = properties["ligand_chain"]
    matrix = DistanceMatrix(u, select =f'chainID {ligand_chain} and name CA').run()

    # Calculate the linkage matrix
    Z = linkage(matrix.results.dist_matrix, method = "average")

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


def main_wf(configuration_path, receptor_pdb_path, ligand_pdb_path, output_path, previous_output_path, skip_until):
    '''
    Main protein-protein docking workflow. Can be used to sample protein-protein docking poses,
    rank them according to pyDock score and cluster them based on RMSD. The best scoring pose from each cluster 
    is saved into a zip file in PDB format.

    Inputs
    ------

        configuration_path   (str): path to input.yml
        receptor_pdb_path    (str): path to receptor pdb file
        ligand_pdb_path      (str): path to ligand pdb file
        output_path          (str): path to output folder
        previous_output_path (str): path to previous output folder. Used if you want to re-use some steps of a previous run
        skip_until           (str): skip everything until this step. Options: makePDB, clustering, oda_filtering 

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

    # Skip step if requested
    run_remaining_steps = run_remaining_steps or skip_until == "makePDB"
    if run_remaining_steps:

        # STEP 6: Generate PDB files for top scoring docking poses
        global_log.info("step6_makePDB: generate PDB files for top scoring docking poses")
        makePDB(**global_paths["step6_makePDB"], properties=global_prop["step6_makePDB"])
    
    else:

        # Link Step 6 folder in previous output path to output path
        link_previous_steps(global_log, ["step6_makePDB"], output_path, previous_output_path)

    # Keep only top scoring docking poses in the ranking data frame
    rank1 = global_prop["step6_makePDB"]["rank1"]
    rank2 = global_prop["step6_makePDB"]["rank2"]
    top_ranking_df = ranking_df[(ranking_df["RANK"] >= rank1) & (ranking_df["RANK"] <= rank2)] 

    # Skip step if requested
    run_remaining_steps = run_remaining_steps or skip_until == "clustering"
    if run_remaining_steps:
       
        # Add docking name and ligand chain to clustering step properties
        global_prop["step7_clustering"]["docking_name"] = global_prop["step6_makePDB"]["docking_name"]
        global_prop["step7_clustering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]

        # STEP 7: Clustering with RMSD
        global_log.info("step7_clustering: clustering with RMSD")
        rmsd_clustering(**global_paths["step7_clustering"], properties=global_prop["step7_clustering"], ranking_df=top_ranking_df) 

    else:

        # Link Step 7 folder in previous output path to output path
        link_previous_steps(global_log, ["step7_clustering"], output_path, previous_output_path)

    # Add receptor and ligand chains to oda_filtering step properties
    global_prop["step8_oda_filtering"]["receptor_chain"] = global_prop["step1_setup"]["receptor"]["newmol"]
    global_prop["step8_oda_filtering"]["ligand_chain"] = global_prop["step1_setup"]["ligand"]["newmol"]

    # STEP 8: Filtering with ODA patches
    global_log.info("step8_oda_filtering: filtering with ODA patches")
    oda_filtering(**global_paths["step8_oda_filtering"], properties=global_prop["step8_oda_filtering"])

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
                        help="Skip everything until this step. Options: makePDB, clustering, oda_filtering",
                        required=False)

    # Output                        
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path=args.config_path,
            receptor_pdb_path=args.receptor_pdb_path,
            ligand_pdb_path=args.ligand_pdb_path,
            output_path=args.output_path,
            previous_output_path=args.previous_output_path,
            skip_until=args.skip_until)
