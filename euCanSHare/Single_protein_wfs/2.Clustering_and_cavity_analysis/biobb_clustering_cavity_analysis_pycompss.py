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

from pycompss.api.api import compss_barrier
from biobb_adapters.pycompss.biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_adapters.pycompss.biobb_gromacs.gromacs.gmxselect import gmxselect
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_adapters.pycompss.biobb_structure_utils.utils.extract_model import extract_model
from biobb_adapters.pycompss.biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_adapters.pycompss.biobb_vs.fpocket.fpocket_filter import fpocket_filter

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
    cavity_analysis_folder = 'step3_cavity_analysis'

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

def check_arguments(global_log, global_paths, traj_path, top_path, clustering_path):
    """
    Check the arguments provided by the user
    """

    # If the user doesn't provide traj_path and top_path or clustering_path 
    if (None in [traj_path, top_path]) and clustering_path is None:

        # And the traj and topology are not given through the configuration file either -> exit
        config_traj_path = global_paths['step1_gmx_cluster'].get('input_traj_path', None)
        config_top_path = global_paths['step0A_make_ndx'].get('input_structure_path', None)

        if (None in [config_traj_path, config_top_path]):
            global_log.error("ERROR: Either clustering_path or both traj_path and top_path must be provided")
            raise SystemExit
        
        # If the user provides traj_path and top_path through the configuration file -> check if they exist
        if not os.path.exists(config_traj_path):
            global_log.error("ERROR: traj_path doesn't exist")
            raise SystemExit

        if not os.path.exists(config_top_path):
            global_log.error("ERROR: top_path doesn't exist")
            raise SystemExit

    # If the user provides traj_path and top_path and clustering_path -> exit
    if (None not in [traj_path, top_path]) and clustering_path is not None:
        global_log.error("ERROR: Both traj_path and top_path and clustering_path are provided, clustering_path will be used")
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


def main_wf(configuration_path, traj_path, top_path, clustering_path, output_path):
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

    # Check arguments
    check_arguments(global_log, global_paths, traj_path, top_path, clustering_path)

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

    # For each model of the receptor
    for cluster_index, cluster_name in enumerate(cluster_names):

        # Create sub folder for the model
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

        # NOTE: this step should be decorated accordingly - not a standard biobb
        # STEP 5: Filter by pocket center of mass 
        global_log.info("step5_filter_residue_com: Filter cavities by center of mass distance to a group of residues") 
        filtered_pockets_IDs = filter_residue_com(**cluster_paths['step5_filter_residue_com'], properties=cluster_prop["step5_filter_residue_com"], global_log=global_log)

        # Update dictionary with filtered pockets
        cluster_filtered_pockets.update({cluster_name : filtered_pockets_IDs})

        # Save model pdb file in sub folder
        model_subfolder = os.path.join(output_path, cluster_name)
        shutil.copyfile(cluster_paths['step3_cavity_analysis']['input_pdb_path'], os.path.join(model_subfolder, 'model.pdb'))

    # NOTE: here all processes should wait for the others to finish, but we cannot use files as not all of them will find pockets after filtering - no output file
    compss_barrier()

    # Create summary with available pockets per cluster 
    global_log.info("    Creating YAML summary file...")
    create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path)

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
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            traj_path = args.traj_path,
            top_path = args.top_path,
            clustering_path = args.clustering_path,
            output_path = args.output_path)