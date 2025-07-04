#!/shared/work/BiobbWorkflows/envs/biobb_sp_cavity_analysis/bin/python

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import time
import glob
import argparse
import csv
import yaml
import json
import shutil
import numpy as np
from pathlib import Path

import MDAnalysis as mda

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_gromacs.gromacs.gmxselect import gmxselect
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_analysis.ambertools.cpptraj_convert import cpptraj_convert
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def is_gromacs_format(traj_path: str) -> bool:
    """
    Checks if the trajectory is in a GROMACS-compatible format (xtc, trr, cpt, g96, gro, pdb, tng)

    Inputs
    ------

        traj_path: 
            path to the trajectory file
    
    Output
    -------

        bool: True if the trajectory is in a GROMACS-compatible format, False otherwise
    """

    # List of GROMACS-compatible formats
    gromacs_formats = ['.xtc', '.trr', '.cpt', '.g96', '.gro', '.pdb', '.tng']

    # Check if the trajectory file is in a GROMACS-compatible format
    if any([traj_path.endswith(format) for format in gromacs_formats]):
        return True
    else:
        return False

def set_gromacs_path(global_prop: dict, binary_path: str) -> None:
    """
    Set the path to the GROMACS binary for all steps using GROMACS.

    Inputs
    ------

        global_prop: Dictionary containing all the properties of the workflow.
        binary_path: Path to the GROMACS binary.
    """

    list_of_steps = ['step1A_traj_preparation_ndx', 'step1B_add_selection_group', 'step2A_strip_traj', 'step2B_strip_top',
                     'step3A_rmsd_calculation_ndx', 'step3B_add_rmsd_group', 'step3C_add_output_group', 'step4_gmx_cluster']

    for step in list_of_steps:
        global_prop[step]['binary_path'] = binary_path

def get_clusters_population(log_path: str, output_path: str, global_log) -> list:
    '''
    Reads the centroids' ID and populations from the log, sorts the clusters by population 
    in descending order and writes the sorted list to a JSON file in the output folder.

    Inputs
    ------

        log_path          : path to log file of clustering step
        output_path       : path to the output folder where the JSON file will be saved
        global_log        : global log object

    Outputs
    -------

        clusters_population (list<tuple>): list with tuples containing (population, cluster_ID) sorted by population
    '''

    # Open log file if it exists
    if os.path.exists(log_path):
        file = open(log_path)
    else:
        global_log.error("Clustering log file not found")
        return

    # Read file
    csv_reader = csv.reader(file)

    cluster_ids = []
    populations = []

    # Parse the log file
    start = False
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

def create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path):
    '''
    Creates 3 sorted summary files with information for each pocket found and filtered for each model. 
    The summary file is a YAML file.

    
    Inputs
    ------

        cluster_names            (list): List with names of clusters
        cluster_populations      (list): List with tuples containing (population, cluster_ID) ordered by population
        cluster_filtered_pockets (dict): Dictionary with filtered pocket's IDs for each cluster
        global_paths             (dict): global paths dictionary
        output_path               (str): path to output folder
    
    Output
    ------

        Info dumped to yaml summary files
    '''

    # Find step names
    cavity_analysis_folder = 'step6_cavity_analysis'

    # Find unfiltered summary file name
    pockets_summary_filename = Path(global_paths[cavity_analysis_folder]['output_summary']).name

    # Dictionary where all available cluster_summary dictionaries will be saved
    global_summary = {}

    # For each cluster
    for cluster_index, cluster_name in enumerate(cluster_names):
        
        # Find list of filtered pocket IDs for this cluster
        filtered_pocket_names = cluster_filtered_pockets[cluster_name]

        # If any pockets are found 
        if len(filtered_pocket_names) > 0:

            # Dictionary with information for this cluster
            cluster_summary = {}

            # If clustering was done externally we might not have this information 
            if cluster_populations is not None:
                # Save population of cluster
                cluster_summary.update({'population' : cluster_populations[cluster_index][0]})

            # Path to this cluster's summary with information for all pocket found
            summary_path = os.path.join(output_path, cluster_name, cavity_analysis_folder, pockets_summary_filename)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # For each filtered pocket
            for pocket_name in filtered_pocket_names:

                # Save entry with pocket information
                cluster_summary.update({pocket_name : all_pockets_summary[pocket_name]})

            # Save pocket IDs
            cluster_summary.update({'pockets' : filtered_pocket_names})

            # Update global_summary
            global_summary.update({cluster_name : cluster_summary})

    # Sort models by 3 criteria (volume, druggability score, score)
    sorted_pockets_by_volume, sorted_pockets_by_drug_score, sorted_pockets_by_score = sort_summary(global_summary)
    
    # Create file names for sorted summary files
    volume_summary_path = os.path.join(output_path, f"summary_by_volume.yml")
    drug_score_summary_path = os.path.join(output_path, f"summary_by_drug_score.yml")
    score_summary_path = os.path.join(output_path, f"summary_by_score.yml")

    # Write the sorted pockets by volume to a YAML file
    with open(volume_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_volume, f, sort_keys = False)

    # Write the sorted pockets by druggability score to a YAML file
    with open(drug_score_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_drug_score, f, sort_keys = False)
    
    # Write the sorted pockets by score to a YAML file
    with open(score_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_score, f, sort_keys = False)

    return 

def filter_residue_com(input_pockets_zip: str, input_pdb_path: str, output_filter_pockets_zip: str, properties: dict, global_log):
    """
    Function that filters pockets by the distance of their center of mass to a group of residues.

    Inputs
    ------

        input_pockets_zip         (str): path to input pockets zip file
        input_pdb_path            (str): path to input pdb with the pocket model (pdb of the receptor)
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

    # Check if input pockets zip file exists
    if not os.path.exists(input_pockets_zip):
        global_log.warning("Input pockets zip file not found, previous step didn't find any pockets or failed")
        return filtered_pocket_IDs

    # Check if step should run
    if not properties['run_step']:

        # Copy input pockets zip file to output filtered pockets zip file
        shutil.copyfile(input_pockets_zip, output_filter_pockets_zip) 

        # Find list of filtered pocket IDs 
        filtered_pocket_IDs = get_pockets_IDs(output_filter_pockets_zip, properties, global_log)

        # Log warning
        global_log.warning("    Skipping step because run_step = False")

        return filtered_pocket_IDs
    
    # Extract all pockets in step folder
    pocket_paths = fu.unzip_list(zip_file=input_pockets_zip, dest_dir=properties['path'])

    # If no pockets are found, return
    if len(pocket_paths) == 0:
        global_log.warning("No pockets found after filtering in previous step")
        return filtered_pocket_IDs
    
    # Load input pdb
    model_universe = mda.Universe(input_pdb_path)

    # Select the residues of interest, e.g. residue number 42
    residue_selection = model_universe.select_atoms(properties['residue_selection'])

    # Compute the center of mass of the selected residues
    residue_com = np.array(residue_selection.center_of_mass())

    # Save paths to pqr files in another list
    pockets_pqr_paths = []
    
    # For each pocket
    for pocket_path in pocket_paths:
        # If path is a pqr file, append to list
        if pocket_path.endswith('.pqr'):
            pockets_pqr_paths.append(pocket_path)

    # Save pockets that pass the distance filter in another list
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
        
        # Erase all the pockets remaining in the step folder
        fu.rm_file_list(file_list=pocket_paths)

        # Log warning
        global_log.warning("No pockets found after filtering")

        return filtered_pocket_IDs
    
    # Zip filtered pockets
    fu.zip_list(zip_file=output_filter_pockets_zip, file_list=filtered_pocket_paths)

    # Erase all the pockets remaining in the step folder
    fu.rm_file_list(file_list=pocket_paths)

    return filtered_pocket_IDs

def sort_summary(pockets_summary: dict):
    """
    Function that reads the dictionary with all models and pockets and sorts the models by:

        1. Volume of the largest pocket in that model
        2. Druggability score of the best pocket in that model
        3. Score of the best pocket in that model
    
    It returns the 3 dictionaries with sorted models and their pockets.

    Inputs
    ------

        pockets_summary (dict): dictionary with all models and pockets
    
    Outputs
    -------

        sorted_pockets_by_volume     (dict): dictionary with models sorted by volume
        sorted_pockets_by_drug_score (dict): dictionary with models sorted by druggability score
        sorted_pockets_by_score      (dict): dictionary with models sorted by score
    
    """

    # Sort the pockets by volume
    sorted_pockets_by_volume = dict(sorted(pockets_summary.items(), key = lambda x: largest_volume(x[1]), reverse = True))

    # Sort the pockets by druggability score
    sorted_pockets_by_drug_score = dict(sorted(pockets_summary.items(), key = lambda x: highest_drug_score(x[1]), reverse = True))

    # Sort the pockets by score
    sorted_pockets_by_score = dict(sorted(pockets_summary.items(), key = lambda x: highest_score(x[1]), reverse = True))
    
    return sorted_pockets_by_volume, sorted_pockets_by_drug_score, sorted_pockets_by_score

def largest_volume(model: dict):
    """
    Function to sort the pockets by volume
    """

    # Find the largest volume
    largest_volume = 0
    for pocket in model["pockets"]:
        if model[pocket]["volume"] > largest_volume:
            largest_volume = model[pocket]["volume"]

    return largest_volume

def highest_drug_score(model: dict):
    """
    Function to sort the pockets by druggability score
    """

    # Find the highest druggability score
    highest_drug_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["druggability_score"] > highest_drug_score:
            highest_drug_score = model[pocket]["druggability_score"]

    return highest_drug_score

def highest_score(model: dict):
    """
    Function to sort the pockets by score
    """

    # Find the highest score
    highest_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["score"] > highest_score:
            highest_score = model[pocket]["score"]

    return highest_score

def get_pockets_IDs(input_pockets_zip: str, properties: dict, global_log):
    """
    Function that retrieves all the pocket IDs from the pockets zip file.

    Inputs
    ------

        input_pockets_zip (str): path to input pockets zip file
        properties       (dict): dictionary with properties for this step
        global_log        (log): global log object
    
    Output
    ------

        filtered_pocket_IDs (list(str)): list with all pocket IDs found in the pockets zip file
    """

    # To return and use to create the summary file
    filtered_pocket_IDs = []

    # Check if input pockets zip file exists
    if not os.path.exists(input_pockets_zip):
        global_log.warning("Input pockets zip file not found, last step didn't find any pockets or failed")
        return filtered_pocket_IDs

    # Extract all pockets in step folder
    pocket_paths = fu.unzip_list(zip_file=input_pockets_zip, dest_dir=properties['path'])

    # If no pockets are found, return
    if len(pocket_paths) == 0:
        global_log.warning("No pockets found after filtering in last step")
        return filtered_pocket_IDs 
    
    # Save paths to pqr files in another list
    pockets_pqr_paths = []

    # For each pocket
    for pocket_path in pocket_paths:
        # If path is a pqr file, append to list
        if pocket_path.endswith('.pqr'):
            pockets_pqr_paths.append(pocket_path)
    
    # Iterate over all pqr files
    for pocket_pqr_path in pockets_pqr_paths:

        # Find pocket ID as file name without "_vert.pqr"
        pocket_ID = Path(pocket_pqr_path).stem.replace("_vert", "")
        filtered_pocket_IDs.append(pocket_ID)
    
    # Erase all the pockets remaining in the step folder
    fu.rm_file_list(file_list=pocket_paths)

    return filtered_pocket_IDs

def check_arguments(global_log, traj_path, top_path, clustering_path):
    """
    Check the arguments provided by the user
    """

    # If the user doesn't provide traj_path and top_path or clustering_path 
    if (None in [traj_path, top_path]) and clustering_path is None:

        global_log.error("ERROR: traj_path and top_path or clustering_path must be provided")
        raise SystemExit

    # If the user provides traj_path and top_path and clustering_path -> exit
    if (None not in [traj_path, top_path]) and clustering_path is not None:
        global_log.error("ERROR: traj_path, top_path and clustering_path are provided, provide either traj_path and top_path or clustering_path")
        raise SystemExit

    # If the user provides traj_path and not top_path -> exit
    if traj_path is not None and top_path is None:
        global_log.error("ERROR: top_path must be provided if traj_path is provided")
        raise SystemExit
    
    # If the user provides top_path and not traj_path -> exit
    if top_path is not None and traj_path is None:
        global_log.error("ERROR: traj_path must be provided if top_path is provided")
        raise SystemExit
    
    # If the user provides traj_path and it doesn't exist -> exit
    if traj_path is not None and not os.path.exists(traj_path):
        global_log.error("ERROR: traj_path doesn't exist")
        raise SystemExit

    # If the user provides top_path and it doesn't exist -> exit
    if top_path is not None and not os.path.exists(top_path):
        global_log.error("ERROR: top_path doesn't exist")
        raise SystemExit
    
    # If the user provides clustering_path and it doesn't exist -> exit
    if clustering_path is not None and not os.path.exists(clustering_path):
        global_log.error("ERROR: clustering_path doesn't exist")
        raise SystemExit

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
# Global properties (common for all steps)
binary_path: gmx                                                                                # GROMACS binary path
working_dir_path: output                                                                        # Workflow default output directory
can_write_console_log: True                                                                     # Output log to console
remove_tmp: True                                                                                # Remove temporal files
restart: True                                                                                   # Do not execute steps if output files are already created
num_clusters: 20                                                                                # Number of most populated clusters to extract from the trajectory and analyze 
                                                                                                # with fpocket, if representative structures are given instead of a traj, this
                                                                                                # number is ignored

# Step 0: Convert from Amber to Gromacs compatible format
# Optional step (will be executed if the trajectory is not in a Gromacs-compatible format)
step0_convert_amber_traj:
  tool: cpptraj_convert
  paths:
    input_traj_path: /path/to/trajectory.dcd                                        # Amber compatible trajectory file
    input_top_path: /path/to/topology.pdb                                           # topology file
    output_cpptraj_path: trajectory.xtc
  properties:
    mask: "all-atoms"                                                               # Any Amber atom selection syntax
    format: "xtc"

# Step 1: Create index file to select some atoms from the trajectory
# Optional step (activate from command line with --prepare_traj)
step1A_traj_preparation_ndx:
  tool: make_ndx 
  paths:
    input_structure_path: dependency/step0_convert_amber_traj/input_top_path
    output_ndx_path: index.ndx
  properties:
    selection: "System"      

step1B_add_selection_group:                                       
  tool: gmxselect
  paths:
    input_structure_path: dependency/step0_convert_amber_traj/input_top_path
    input_ndx_path: dependency/step1A_traj_preparation_ndx/output_ndx_path
    output_ndx_path: index_selection.ndx
  properties:
    selection: '"Selection" resnr 1 to 196'                                       # Gromacs selection syntax
    append: True        

# Step 2: Extract requested atoms from the Gromacs compatible trajectory and topology
# Optional step (activate from command line with --prepare_traj)
step2A_strip_traj:
  tool: gmx_trjconv_trj
  paths: 
    input_traj_path: dependency/step0_convert_amber_traj/output_cpptraj_path
    input_top_path: dependency/step0_convert_amber_traj/input_top_path
    input_index_path: dependency/step1B_add_selection_group/output_ndx_path
    output_traj_path: trajectory.xtc
  properties:
    selection: "Selection" 
    start: 0
    end: 100 
    dt: 1

step2B_strip_top:                       
  tool: gmx_trjconv_str
  paths: 
    input_structure_path: dependency/step0_convert_amber_traj/input_top_path
    input_top_path: dependency/step0_convert_amber_traj/input_top_path
    input_index_path: dependency/step1B_add_selection_group/output_ndx_path
    output_str_path: stripped_topology.pdb
  properties:
    selection: "Selection"   

# Step 3: Create index file to select the atoms for the RMSD calculation
step3A_rmsd_calculation_ndx:
  tool: make_ndx 
  paths:
    input_structure_path: dependency/step2B_strip_top/output_str_path
    output_ndx_path: index.ndx
  properties:
    selection: "System"      

step3B_add_rmsd_group:                                       
  tool: gmxselect
  paths:
    input_structure_path: dependency/step2B_strip_top/output_str_path
    input_ndx_path: dependency/step3A_rmsd_calculation_ndx/output_ndx_path
    output_ndx_path: index_rmsd.ndx
  properties:
    selection: '"RmsdGroup" resnr 181 to 296' 
    append: True        

step3C_add_output_group:                                       
  tool: gmxselect
  paths:
    input_structure_path: dependency/step2B_strip_top/output_str_path
    input_ndx_path: dependency/step3B_add_rmsd_group/output_ndx_path
    output_ndx_path: index_rmsd_output.ndx
  properties:
    selection: '"OutputGroup" group "Protein"'
    append: True

# Steps 4-5: Cluster trajectory and extract centroids pdb
step4_gmx_cluster:
  tool: gmx_cluster
  paths:
    input_traj_path: dependency/step2A_strip_traj/output_traj_path
    input_structure_path: dependency/step2B_strip_top/output_str_path
    input_index_path: dependency/step3C_add_output_group/output_ndx_path
    output_pdb_path: output.cluster.pdb
    output_cluster_log_path: output.cluster.log
    output_rmsd_cluster_xpm_path: output.rmsd-clust.xpm
    output_rmsd_dist_xvg_path: output.rmsd-dist.xvg
  properties:
    fit_selection: RmsdGroup       
    output_selection: OutputGroup        
    dista: False
    method: linkage                      
    cutoff: 0.10                    # nm (RMSD cut-off)
    nofit: True                     # Wether to use the RmsdGroups to fit the traj before computing the RMSD or not

step5_extract_models:
  tool: extract_model
  paths:
    input_structure_path: dependency/step4_gmx_cluster/output_pdb_path
    output_structure_path: cluster.pdb     
  properties:

# Step 6-8: Cavity analysis with fpocket on centroids + filtering
step6_cavity_analysis:
  tool: fpocket_run
  paths:
    input_pdb_path: dependency/step5_extract_models/output_structure_path
    output_pockets_zip: all_pockets.zip
    output_summary: summary.json
  properties:
    min_radius: 3
    max_radius: 6
    num_spheres: 35
    sort_by: druggability_score

step7_filter_cavities:
  tool: fpocket_filter
  paths:
    input_pockets_zip: dependency/step6_cavity_analysis/output_pockets_zip
    input_summary: dependency/step6_cavity_analysis/output_summary
    output_filter_pockets_zip: filtered_pockets.zip
  properties:
    score: [0.4, 1]
    druggability_score: [0.4, 1]
    volume: [200, 5000]

step8_filter_residue_com:
  paths: 
    input_pockets_zip: dependency/step7_filter_cavities/output_filter_pockets_zip
    input_pdb_path: dependency/step5_extract_models/output_structure_path
    output_filter_pockets_zip: filtered_pockets.zip
  properties:
    residue_selection: "resid 31 or resid 21"      # MDAnalysis selection string
    distance_threshold: 8                          # Distance threshold in Angstroms (6-8 are reasonable values if the residue/s are part of the pocket)
    run_step: False                                 # Run step or not
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
        
# Main workflow   
def main_wf(configuration_path: str,  
            traj_path: str, 
            top_path: str, 
            clustering_path: str, 
            distance_threshold: float, 
            prepare_traj: bool, 
            filtering_selection: str, 
            output_path: str
    ):
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and analyzes the cavities of the most representative
    structures. Then filters the cavities according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        configuration_path: 
            path to YAML configuration file
        traj_path:  
            path to trajectory file
        top_path:  
            path to topology file
        clustering_path:  
            path to the folder with the most representative structures in pdb format from an external clustering 
        prepare_traj:  
            flag to prepare the trajectory for clustering
        filtering_selection:  
            residue selection to filter pockets by distance to center of mass
        distance_threshold:  
            distance threshold to filter pockets by distance to center of mass
        output_path:  
            path to output folder

    Outputs
    -------

        /output folder

        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Default configuration file
    default_config = False
    if configuration_path is None:
        default_config = True
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
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Enforce gromacs binary path for all steps using gromacs
    if conf.properties.get('binary_path'):
        global_log.info(f"Using GROMACS binary path: {conf.properties['binary_path']}")
        set_gromacs_path(global_prop, conf.properties['binary_path'])

    # Check arguments
    check_arguments(global_log, traj_path, top_path, clustering_path)

    # If clustering is not given externally -> cluster the input trajectory
    if clustering_path is None:

        # If the trajectory is not in a GROMACS-compatible format, convert it
        if not is_gromacs_format(traj_path):

            # Enforce traj and top paths
            global_paths['step0_convert_amber_traj']['input_traj_path'] = traj_path
            global_paths['step0_convert_amber_traj']['input_top_path'] = top_path

            # STEP 0: Convert the trajectory to xtc format
            global_log.info("step0_convert_amber_traj: Converting AMBER trajectory to xtc format")
            cpptraj_convert(**global_paths['step0_convert_amber_traj'], properties=global_prop['step0_convert_amber_traj'])

            # Change the traj path to the converted version of the trajectory
            traj_path = global_paths['step0_convert_amber_traj']['output_cpptraj_path']
        
        # If the trajectory is in a GROMACS-compatible format, pass the path to the next step
        else:

            # Enforce traj path to next step
            global_paths['step2A_strip_traj']['input_traj_path'] = traj_path

        # If the user wants to prepare the trajectory before clustering
        if prepare_traj:

            # Enforce traj and top paths
            global_paths['step1A_traj_preparation_ndx']['input_structure_path'] = top_path
            global_paths['step1B_add_selection_group']['input_structure_path'] = top_path
            global_paths['step2A_strip_traj']['input_top_path'] = top_path
            global_paths['step2B_strip_top']['input_structure_path'] = top_path
            global_paths['step2B_strip_top']['input_top_path'] = top_path

            # STEP 1A: Create index file for atoms selection
            global_log.info("step1A_traj_preparation_ndx: Creation of index file for atoms selection")
            make_ndx(**global_paths['step1A_traj_preparation_ndx'], properties=global_prop['step1A_traj_preparation_ndx'])

            # STEP 1B: Add selection group to index file
            global_log.info("step1B_add_selection_group: Adding selection group to index file")
            gmxselect(**global_paths['step1B_add_selection_group'], properties=global_prop['step1B_add_selection_group'])

            # STEP 2A: Strip trajectory
            global_log.info("step2A_strip_traj: Stripping trajectory, keeping just selected atoms")
            gmx_trjconv_trj(**global_paths['step2A_strip_traj'], properties=global_prop['step2A_strip_traj'])

            # STEP 2B: Strip topology
            global_log.info("step2B_strip_top: Stripping topology, keeping just selected atoms")
            gmx_trjconv_str(**global_paths['step2B_strip_top'], properties=global_prop['step2B_strip_top'])

        # Trajectory and topology are already prepared (fitted, with no strange atoms and in GROMACS-compatible format)
        else:
            
            # Enforce traj and top paths. 
            global_paths['step4_gmx_cluster']['input_traj_path'] = traj_path
            global_paths['step3A_rmsd_calculation_ndx']['input_structure_path'] = top_path
            global_paths['step3B_add_rmsd_group']['input_structure_path'] = top_path
            global_paths['step3C_add_output_group']['input_structure_path'] = top_path
            global_paths['step4_gmx_cluster']['input_structure_path'] = top_path

        # STEP 3A: Create index file for rmsd calculation
        global_log.info("step3A_rmsd_calculation_ndx: Creation of index file")
        make_ndx(**global_paths['step3A_rmsd_calculation_ndx'], properties=global_prop['step3A_rmsd_calculation_ndx'])

        # STEP 3B: Add rmsd group to index file
        global_log.info("step3B_add_rmsd_group: Adding RmsdGroup to index file")
        gmxselect(**global_paths['step3B_add_rmsd_group'], properties=global_prop['step3B_add_rmsd_group'])

        # STEP 3C: Add output group to index file
        global_log.info("step3C_add_output_group: Adding OutputGroup to index file")
        gmxselect(**global_paths['step3C_add_output_group'], properties=global_prop['step3C_add_output_group'])

        # STEP 4: Cluster trajectory with gmx_cluster
        global_log.info("step4_gmx_cluster: Clustering structures from the trajectory")
        gmx_cluster(**global_paths["step4_gmx_cluster"], properties=global_prop["step4_gmx_cluster"])

        # Save centroid IDs and populations in JSON file
        global_log.info( "step4_gmx_cluster: Reading clustering outcome, generating clusters JSON file")
        cluster_populations = get_clusters_population(log_path = global_paths["step4_gmx_cluster"]['output_cluster_log_path'],
                                                      output_path = global_prop["step4_gmx_cluster"]['path'],
                                                      global_log = global_log)

        # Number of clusters: minimum between number of clusters requested and number of clusters obtained
        num_clusters = min(conf.properties['num_clusters'], len(cluster_populations))

        # Cluster names are the cluster IDs
        cluster_names = [str(cluster_populations[i][1]) for i in range(num_clusters)]
    
    # If clustering is given externally
    else:

        # Obtain the full sorted list of pdb files from clustering path
        # If the clustering path is a file, we assume it is a single pdb file
        if os.path.isfile(clustering_path):
            global_log.info("External clustering file provided")
            pdb_paths = [clustering_path]
        else:
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

    # For representative structure (model)
    for cluster_index, cluster_name in enumerate(cluster_names):

        # Create sub folder for the model
        cluster_prop = conf.get_prop_dic(prefix=cluster_name)
        cluster_paths = conf.get_paths_dic(prefix=cluster_name)

        # Enforce filtering selection
        if filtering_selection is not None:
            cluster_prop['step8_filter_residue_com']['residue_selection'] = filtering_selection 
            
        # If clustering was done here, extract the model from the clustering results
        if clustering_path is None:

            # Update input structures path and model index
            cluster_paths['step5_extract_models']['input_structure_path'] = global_paths['step4_gmx_cluster']['output_pdb_path']
            cluster_prop['step5_extract_models']['models'] = [cluster_populations[cluster_index][1]]

            # STEP 5: Extract one model from the input structures path
            extract_model(**cluster_paths['step5_extract_models'], properties=cluster_prop['step5_extract_models'])

        # If clustering was done externally, just update the input pdb path
        else:

            cluster_paths['step6_cavity_analysis']['input_pdb_path'] = pdb_paths[cluster_index]
            cluster_paths['step8_filter_residue_com']['input_pdb_path'] = pdb_paths[cluster_index]
        
        # STEP 6: Cavity analysis
        global_log.info("step6_cavity_analysis: Compute protein cavities using fpocket")
        fpocket_run(**cluster_paths['step6_cavity_analysis'], properties=cluster_prop["step6_cavity_analysis"])

        # STEP 7: Filtering cavities
        global_log.info("step7_filter_cavities: Filter found cavities")
        fpocket_filter(**cluster_paths['step7_filter_cavities'], properties=cluster_prop["step7_filter_cavities"])
        
        # Enforce distance threshold
        if distance_threshold:
            cluster_prop['step8_filter_residue_com']['distance_threshold'] = distance_threshold

        # STEP 8: Filter by pocket center of mass 
        global_log.info("step8_filter_residue_com: Filter cavities by center of mass distance to a group of residues") 
        filtered_pockets_IDs = filter_residue_com(**cluster_paths['step8_filter_residue_com'], properties=cluster_prop["step8_filter_residue_com"], global_log=global_log)

        # Update dictionary with filtered pockets
        cluster_filtered_pockets.update({cluster_name : filtered_pockets_IDs})

        # Save model pdb file in sub folder
        model_subfolder = os.path.join(output_path, cluster_name)
        shutil.copyfile(cluster_paths['step6_cavity_analysis']['input_pdb_path'], os.path.join(model_subfolder, 'model.pdb'))

    # Create summary with available pockets per cluster 
    global_log.info("    Creating YAML summary file...")
    create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path)

    if default_config:
        # Move the default configuration file to the output path
        shutil.move(configuration_path, os.path.join(output_path, 'config.yml'))
        configuration_path = os.path.join(output_path, 'config.yml')
        
    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow name: Cavity analysis')
    global_log.info('  Output path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=False)
    
    parser.add_argument('--traj_path', dest='traj_path',
                        help="Path to input trajectory (GROMACS or AMBER formats)", 
                        required=False)

    parser.add_argument('--top_path', dest='top_path',
                        help="Path to input structure (gro, pdb)",
                        required=False) 

    parser.add_argument('--clustering_path', dest='clustering_path',
                        help="Input path to representative structures (folder with pdb files)", 
                        required=False)

    parser.add_argument('--prepare_traj', action='store_true',
                        help="Prepare trajectory for clustering (activate if you need to delete atoms or you don't have a gromacs compatible trajectory file)",
                        required=False)

    parser.add_argument('--filtering_selection', dest='filtering_selection',
                        help="Residue selection to filter pockets by distance to center of mass",
                        required=False)
    
    parser.add_argument('--distance_threshold', dest='distance_threshold', type=float,
                        help="Distance threshold to filter pockets by distance to center of mass",
                        required=False)
    
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            traj_path = args.traj_path,
            top_path = args.top_path,
            clustering_path = args.clustering_path,
            prepare_traj = args.prepare_traj,
            filtering_selection = args.filtering_selection,
            distance_threshold = args.distance_threshold,
            output_path = args.output_path)