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

import MDAnalysis as mda
from MDAnalysis.analysis.diffusionmap import DistanceMatrix
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

def prepare_step5(paths: dict, prop: dict):
    """ Create step5_clustering folder, create top docking poses subfolder and unzip top docking poses into it"""
    
    # Create step5_clustering folder
    fu.create_dir(prop["path"])

    # Create subfolder for top poses
    top_poses_path = fu.create_dir(str(Path(prop["path"]).joinpath("top_poses")))

    # Extract the file with the top docking poses into the subfolder
    pdb_paths = fu.unzip_list(zip_file = paths["input_zip_path"], dest_dir = top_poses_path)

    return pdb_paths, top_poses_path

def prepare_top_poses(pdb_paths, top_poses_path, top_ranking_df, docking_name):
    """ Change pdb file names from conformation number to rank. Return a list with the new file paths and a list with the corresponding rank for each file
    """
    paths_ranking = []
    new_pdb_paths = []

    # Change pdb file name for rank and create a list with the ranking
    for pdb_path in pdb_paths:

        # Get the conformation number from the file name
        conformation_number = Path(pdb_path).stem.strip(f"{docking_name}_")

        # Get the rank from the ranking data frame
        rank = top_ranking_df[top_ranking_df["Conf"] == int(conformation_number)]["RANK"].values[0]

        # Add the rank to the list
        paths_ranking.append(f"{rank}")

        # Create the new file path
        new_pdb_path = str(Path(top_poses_path).joinpath(f"{rank}.pdb"))

        # Rename the file with the new file path
        os.rename(pdb_path, new_pdb_path)

        # Add the new file path to the list
        new_pdb_paths.append(new_pdb_path)
    
    return new_pdb_paths, paths_ranking

def plot_dendogram(linkage_matrix, labels, global_prop):
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
        color_threshold = global_prop["step5_clustering"]["rmsd_threshold"],
        labels = labels
        )
    plt.savefig(str(Path(global_prop["step5_clustering"]["path"]).joinpath("dendrogram.png")))
    plt.close() 

def main_wf(configuration_path):
    '''
    Main protein-protein docking workflow. Can be used to sample protein-protein docking poses,
    rank them according to pyDock score and cluster them based on RMSD. The best scoring pose from each cluster 
    is saved into a zip file in PDB format.

    Inputs
    ------

        configuration_path (str): path to input.yml

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    start_time = time.time()
    
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Sampling protein-protein docking poses
    
    # STEP 1
    global_log.info("step1_setup: setup receptor and ligand proteins for pyDock")
    setup(**global_paths["step1_setup"], properties=global_prop["step1_setup"])
    
    # STEP 2
    global_log.info("step2_ftdock: sample docking poses using ftdock (FFT-based algorithm)")
    ftdock(**global_paths["step2_ftdock"], properties=global_prop["step2_ftdock"])

    # Scoring protein-protein docking poses

    # STEP 3
    global_log.info("step3_dockser: score docking poses using pyDock")
    dockser(**global_paths["step3_dockser"], properties=global_prop["step3_dockser"])

    # Read csv with ranking of docking poses (pydock score)
    ranking_df = pd.read_csv(global_paths["step3_dockser"]["output_ene_path"], sep='\s+', skiprows=[1], header=0)

    # STEP 4
    global_log.info("step4_makePDB: generate PDB files for top scoring docking poses")
    makePDB(**global_paths["step4_makePDB"], properties=global_prop["step4_makePDB"])
    
    # Keep only docking poses between rank1 and rank2
    rank1 = global_prop["step4_makePDB"]["rank1"]
    rank2 = global_prop["step4_makePDB"]["rank2"]
    top_ranking_df = ranking_df[(ranking_df["RANK"] >= rank1) & (ranking_df["RANK"] <= rank2)]

    # Clustering with RMSD

    # STEP 5
    global_log.info("step5_clustering: clustering with RMSD")

    # Prepare step5_clustering folder - unzip top docking poses into subfolder
    conf_pdb_paths, top_poses_path = prepare_step5(global_paths["step5_clustering"], global_prop["step5_clustering"])

    # Prepare top docking poses - change file names from conformation number to rank and create a list with the corresponding rank for each file
    pdb_paths, paths_ranking = prepare_top_poses(conf_pdb_paths, top_poses_path, top_ranking_df, global_prop['step4_makePDB']['docking_name'])

    # Calculate the RMSD matrix using only CA ligand atoms without alignment - it could be parallelized with MDAnalysis.analysis.encore.confdistmatrix if needed
    u = mda.Universe(pdb_paths[0], pdb_paths)
    ligand_chain = global_prop["step1_setup"]["ligand"]["newmol"]
    matrix = DistanceMatrix(u, select =f'chainID {ligand_chain} and name CA').run()

    # Calculate the linkage matrix
    Z = linkage(matrix.results.dist_matrix, method = "average")

    # Get the cluster labels from the linkage matrix
    cluster_labels = fcluster(Z, t = global_prop["step5_clustering"]["rmsd_threshold"], criterion = "distance")

    plot_dendogram(Z, paths_ranking, global_prop)

    # Find the best ranking pose path for each cluster
    top_distinct_poses_paths = []
    for cluster in range(1, max(cluster_labels) + 1):

        # Get the ranking of the poses in the cluster
        cluster_ranking = [paths_ranking[i] for i in np.where(cluster_labels == cluster)[0]]

        # Sort the ranking
        cluster_ranking.sort()

        # Get the best ranking pose path
        top_distinct_poses_paths.append(str(Path(top_poses_path).joinpath(f"{cluster_ranking[0]}.pdb")))

    # Zip the best ranking pose paths into a zip file
    fu.zip_list(zip_file = global_paths["step5_clustering"]["output_zip_path"], file_list = top_distinct_poses_paths)

    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
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

    args = parser.parse_args()

    main_wf(configuration_path=args.config_path)
