##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import re
import time
import glob
import argparse
import csv
import yaml
import json
import numpy as np
from pathlib import Path

import MDAnalysis as mda

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_gromacs.gromacs.gmxselect import gmxselect
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def get_clusters_population(log_path: str, output_path: str, global_log) -> list:
    '''
    Reads the centroids' ID and populations from the clustering log, sorts the clusters by 
    population in descending order and writes the sorted list to a JSON file in the clustering_folder.

    Inputs
    -------

        log_path          : log path of clustering step
        output_path       : path to the output folder where the population of each cluster will be saved in a JSON file
        global_log        : global log object

    Outputs
    -------

        /json_file  with sorted cluster ids and populations

        clusters_population (list<tuple>): list with tuples containing (population, cluster_ID) sorted by population
    '''

    if os.path.exists(log_path):
        file = open(log_path)
    # If not found, return error
    else:
        global_log.error("Clustering log file not found")
        return

    # Read file
    csv_reader = csv.reader(file)

    cluster_ids = []
    populations = []

    # Condition to wait until we reach clustering information in file
    start = False

    # For each row in file
    for row in csv_reader:

        # When we have reached clustering information
        if start: 
            # Split the whole row into columns using '|' as separators
            col = row[0].split('|')

            # If the first column contains something more than whitespace
            if len(col[0].strip()): 
                # Save cluster id and population
                cluster_ids.append(col[0].strip())
                populations.append(col[1].strip().split()[0])

        # We have reached beginning of clustering information, start saving
        if len(row) and row[0].startswith('cl.'):
            start = True
        
    # Close log file
    file.close()
    
    # Convert to integers
    populations = [int(x) for x in populations]
    cluster_ids = [int(x) for x in cluster_ids]

    # Sort clusters by population
    clusters_population = sorted(zip(populations, cluster_ids), reverse=True)

    ordered_clusters = []

    # Save sorted cluster ids and populations in dictionaries and populations in list
    for population,index in clusters_population:
        ordered_clusters.append({
            'cluster': index,
            'population': population
            })
    
    # Save list in JSON file
    with open(os.path.join(output_path, "clusters.json"), 'w') as outfile:
        json.dump(ordered_clusters, outfile)

    outfile.close()

    # If the number of clusters is very large, issue a warning - the user might want to increase the RMSD cutoff
    if len(cluster_ids)>100:
        global_log.warning(f"   Warning: Large number of clusters found. Consider increasing the RMSD cutoff.")
        global_log.warning(f"   Warning: Number of clusters: {len(cluster_ids)}")
        global_log.warning(f"   Warning: Number of clusters with more than 1 member: {len([x for x in populations if x > 1])}")
        global_log.warning(f"   Warning: Large number of clusters might make model extraction very slow.")
    
    # Return list with cluster population and id sorted by population
    return clusters_population

def create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path, output_summary_path = None):
    '''
    Print in log file all available pockets after filtering for each centroid. 
    Include also information for each pocket and population for each centroid if available.
    
    Inputs
    ------

        cluster_names         (list): List with names of centroids
        cluster_populations   (list): List with tuples containing (population, cluster_ID) ordered by population
        cluster_filtered_pockets (dict): Dictionary with filtered pockets for each centroid
        global_paths          (dict): global paths dictionary
        output_path            (str): path to output folder
        output_summary_path    (str): path to output summary file
    
    Output
    ------

        Info dumped to yaml summary file
    '''

    # Find step names
    cavity_analysis_folder = 'step3_cavity_analysis'

    # Find summary file name
    pockets_summary_filename = Path(global_paths[cavity_analysis_folder]['output_summary']).name

    # Dictionary where all available cluster_summary dictionaries will be saved
    global_summary = {}

    # For each centroid
    for cluster_index, cluster_name in enumerate(cluster_names):
        
        # Find list of pocket names for this cluster
        pocket_names = cluster_filtered_pockets[cluster_name]

        # If any pockets are found 
        if len(pocket_names) > 0:

            # Dictionary with information for this cluster
            cluster_summary = {}

            # If clustering was done externally we might not have this information 
            if cluster_populations is not None:

                # Save population of centroid
                cluster_summary.update({'population' : cluster_populations[cluster_index][0]})

            # Path to this cluster's summary with information for all pocket found
            summary_path = os.path.join(output_path, cluster_name, cavity_analysis_folder, pockets_summary_filename)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # For each filtered pocket
            for pocket_name in pocket_names:

                # Save entry with pocket information
                cluster_summary.update({pocket_name : all_pockets_summary[pocket_name]})

            # Save pocket IDs
            cluster_summary.update({'pockets' : pocket_names})

            # Update global_summary
            global_summary.update({cluster_name : cluster_summary})

    # Save global summary in YAML file
    if output_summary_path is None:
        output_summary_path = os.path.join(output_path, 'summary.yml')
    with open(output_summary_path, 'w') as outfile:
        yaml.dump(global_summary, outfile)
    outfile.close()

    return 

def filter_residue_com(input_pockets_zip: str, input_pdb_path: str, output_filter_pockets_zip: str, properties: dict, global_log):
    """
    Function that filters pockets by the distance of their center of mass to a group of residues.

    Inputs
    ------

        input_pockets_zip         (str): path to input pockets zip file
        input_pdb_path            (str): path to input pdb with the pocket model 
        output_filter_pockets_zip (str): path to filtered pockets zip file
        properties               (dict): dictionary with properties for this step
        global_log                (log): global log object
    
    Output
    ------

        filtered_pocket_IDs (list(str)): list with pocket IDs that passed the filter
    """


    # To return and use to create the summary file
    filtered_pocket_IDs = []

    # Create step folder
    fu.create_dir(properties['path'])

    # Extract all pockets in step folder
    pocket_paths = fu.unzip_list(zip_file=input_pockets_zip, dest_dir=properties['path'])

    # If no pockets are found, return
    if len(pocket_paths) == 0:
        global_log.warning("No pockets found after filtering")
        # Erase all the pockets remaining in the step folder
        fu.rm_file_list(file_list=pocket_paths)
        return filtered_pocket_IDs
    
    # Load input pdb
    model_universe = mda.Universe(input_pdb_path)

    # Select the residue of interest, e.g., residue number 42
    residue_selection = model_universe.select_atoms(properties['residue_selection'])

    # Compute the center of mass of the selected residue
    residue_com = np.array(residue_selection.center_of_mass())

    # Save paths to pqr files in another list
    pockets_pqr_paths = []
    
    # For each pocket
    for pocket_path in pocket_paths:

        # If path is a pqr file, append to list
        if pocket_path.endswith('.pqr'):
            pockets_pqr_paths.append(pocket_path)

    # To zip only the filtered pockets
    filtered_pocket_paths = []

    # Iterate over all pqr files
    for pocket_pqr_path in pockets_pqr_paths:

        # Find pocket ID as file name without "_vert.pqr"
        pocket_ID = Path(pocket_pqr_path).stem.replace("_vert", "")

        # Load pocket
        pocket_universe = mda.Universe(pocket_pqr_path)

        # Select all atoms in pocket
        pocket_selection = pocket_universe.select_atoms('all')

        # Compute the center of mass of the pocket
        pocket_com = np.array(pocket_selection.center_of_mass())

        # Compute distance between pocket and residue center of mass using numpy
        distance = np.linalg.norm(pocket_com - residue_com)

        # Save pocket if distance is smaller than threshold
        if distance < properties['distance_threshold']:

            # Save pocket
            filtered_pocket_IDs.append(pocket_ID)
            filtered_pocket_paths.append(pocket_pqr_path)

    # If no pockets are found, return
    if len(filtered_pocket_paths) == 0:
        global_log.warning("No pockets found after filtering")
        # Erase all the pockets remaining in the step folder
        fu.rm_file_list(file_list=pocket_paths)
        return filtered_pocket_IDs
    
    # Zip filtered pockets
    fu.zip_list(zip_file=output_filter_pockets_zip, file_list=filtered_pocket_paths)

    # Erase all the pockets remaining in the step folder
    fu.rm_file_list(file_list=pocket_paths)

    return filtered_pocket_IDs


def main_wf(configuration_path, traj_path, top_path, clustering_path, output_path, output_summary_path):
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and analyzes the cavities of the most representative
    structures. Then filters the cavities according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        configuration_path  (str): path to YAML configuration file
        traj_path           (str): (Optional) path to trajectory file
        top_path            (str): (Optional) path to topology file
        clustering_path     (str): (Optional) path to the folder with the most representative structures in pdb format from an external clustering 
        output_path         (str): (Optional) path to output folder
        output_summary_path (str): (Optional) path to output summary file

    Outputs
    -------

        /output folder

        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties
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
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # If clustering is not given externally -> cluster the input trajectory and return cluster path
    if clustering_path is None:
        
        # Enforce traj_path if provided
        if traj_path is not None:
            global_paths['step1_gmx_cluster']['input_traj_path'] = traj_path

        # Enforce top_path if provided
        if top_path is not None:
            global_paths['step0A_make_ndx']['input_structure_path'] = top_path
            global_paths['step0B_add_rmsd_group']['input_structure_path'] = top_path
            global_paths['step0C_add_output_group']['input_structure_path'] = top_path
            global_paths['step1_gmx_cluster']['input_structure_path'] = top_path

        # STEP 0: Create index file
        global_log.info("step0A_make_ndx: Creation of index file")
        make_ndx(**global_paths['step0A_make_ndx'], properties=global_prop['step0A_make_ndx'])

        # STEP 0: Add RmsdGroup to index file
        global_log.info("step0B_add_rmsd_group: Adding RmsdGroup to index file")
        gmxselect(**global_paths['step0B_add_rmsd_group'], properties=global_prop['step0B_add_rmsd_group'])

        # STEP 0: Add OutputGroup to index file
        global_log.info("step0C_add_output_group: Adding OutputGroup to index file")
        gmxselect(**global_paths['step0C_add_output_group'], properties=global_prop['step0C_add_output_group'])

        # STEP 1: Cluster trajectory with gmx_cluster
        global_log.info("step1_gmx_cluster: Clustering structures from the trajectory")
        gmx_cluster(**global_paths["step1_gmx_cluster"], properties=global_prop["step1_gmx_cluster"])

        # Save centroid IDs and populations in JSON file
        global_log.info( "step1_gmx_cluster: Reading clustering outcome, generating clusters JSON file")
        cluster_populations = get_clusters_population(log_path = global_paths["step1_gmx_cluster"]['output_cluster_log_path'],
                                                      output_path = global_prop["step1_gmx_cluster"]['path'],
                                                      global_log = global_log)

        # Number of clusters: minimum between number of clusters requested and number of clusters obtained
        num_clusters = min(conf.properties['num_clusters'], len(cluster_populations))

        # Cluster names are the cluster IDs
        cluster_names = [str(cluster_populations[i][1]) for i in range(num_clusters)]
    else:

        # Obtain the full sorted list of pdb files from clustering path
        pdb_paths = sorted(glob.glob(os.path.join(clustering_path,"*.pdb")))

        # Population information will not be available in this case
        cluster_populations = None

        # Number of clusters: number of pdb files
        num_clusters = len(pdb_paths)

        # Cluster names are the pdb file names
        cluster_names = [Path(pdb_path).stem for pdb_path in pdb_paths]
    
    global_log.info(f"Number of models to analyze: {num_clusters}")

    # Dictionary to save the filtered pocket IDs for each model
    cluster_filtered_pockets = {}

    for cluster_index, cluster_name in enumerate(cluster_names):

        cluster_prop = conf.get_prop_dic(prefix=cluster_name)
        cluster_paths = conf.get_paths_dic(prefix=cluster_name)

        if clustering_path is None:

            # Update input structure path and model index of step 2
            cluster_paths['step2_extract_models']['input_structure_path'] = global_paths['step1_gmx_cluster']['output_pdb_path']
            cluster_prop['step2_extract_models']['models'] = [cluster_populations[cluster_index][1]]

            # STEP 2: Extract one model from previous clustering
            extract_model(**cluster_paths['step2_extract_models'], properties=cluster_prop['step2_extract_models'])

        else:

            # Update input structure path of step3
            cluster_paths['step3_cavity_analysis']['input_pdb_path'] = pdb_paths[cluster_index]
            cluster_paths['step5_filter_residue_com']['input_pdb_path'] = pdb_paths[cluster_index]
        
        # STEP 3: Cavity analysis
        global_log.info("step3_cavity_analysis: Compute protein cavities using fpocket")
        fpocket_run(**cluster_paths['step3_cavity_analysis'], properties=cluster_prop["step3_cavity_analysis"])

        # STEP 4: Filtering cavities
        global_log.info("step4_filter_cavities: Filter found cavities")
        fpocket_filter(**cluster_paths['step4_filter_cavities'], properties=cluster_prop["step4_filter_cavities"])

        # STEP 5: Filter by pocket center of mass 
        global_log.info("step5_filter_residue_com: Filter cavities by center of mass distance to a group of residues") 
        filtered_pockets_IDs = filter_residue_com(**cluster_paths['step5_filter_residue_com'], properties=cluster_prop["step5_filter_residue_com"], global_log=global_log)

        # Update dictionary with filtered pockets
        cluster_filtered_pockets.update({cluster_name : filtered_pockets_IDs})

    # Create summary with available pockets per cluster 
    global_log.info("    Creating YAML summary file...")
    create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path, output_summary_path)

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)
    
    parser.add_argument('--traj_path', dest='traj_path',
                        help="Path to input trajectory (xtc, trr, cpt, gro, g96, pdb, tng)", 
                        required=False)

    parser.add_argument('--top_path', dest='top_path',
                        help="Path to input structure (gro, pdb)",
                        required=False) 

    parser.add_argument('--clustering_path', dest='clustering_path',
                        help="Input path to representative structures (folder with pdb files)", 
                        required=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    parser.add_argument('--output_summary', dest='output_summary_path',
                        help="Output summary path (default: /output_path/summary.yml)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            traj_path = args.traj_path,
            top_path = args.top_path,
            clustering_path = args.clustering_path,
            output_path = args.output_path,
            output_summary_path = args.output_summary_path)