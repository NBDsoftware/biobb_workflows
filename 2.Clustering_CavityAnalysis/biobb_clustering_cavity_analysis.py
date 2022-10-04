##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import time
import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path, PurePath

from pkg_resources import resource_string
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter



def addModelSuffix(original_path, i):
    '''
    Adds '_modeli' to original_path before file extension. For example:

    original_path = Path('/home/user/ClusteringWF/output/step2_extract_models/output.pdb')
    i = 0

    new_path = '/home/user/ClusteringWF/output/step2_extract_models/output_model0.pdb'

    Inputs
        original_path  (Path class from pathlib):  path to original file including filename
        i  (int):  integer 
    '''

    filename = original_path.name
    directory = original_path.parent

    filename_parts = filename.split(".")

    new_path = str(directory) + "/" + filename_parts[0] + "_model" + str(i) + "." + filename_parts[1]

    return new_path


def moveFileIfItExists(origin_src, dest_src):
    '''
    Checks the existence of 'origin_src' file and moves it to 'dest_src' if it exists

    Inputs
        origin_src  (str):  path to original file including filename
        dest_src    (str):  path to destination file including filename
    '''

    if (os.path.exists(origin_src)):
        shutil.move(origin_src, dest_src)
    
    return



def saveClusters(log_file, json_file, clustering_folder):
    '''
    Reads the clusters' id and population from log_file, sorts the clusters by 
    population in descending order and writes the sorted list to a JSON file.

    Inputs
        log_file          (str): log file from clustering (reads from here)
        json_file         (str): path for JSON file with sorted cluster ids and populations (writes here) 
        clustering_folder (str): str with name of clustering step folder - just in case Restart 
                                 has been set to True and files have been moved

    Outputs:
        clusters_iterator (list<tuple>): list with cluster information ordered by population
        (implicitly) JSON file with sorted cluster ids and populations
    '''

    # Open log file - try to find it in working dir or in step dir
    if (os.path.exists(log_file)):
        file = open(log_file)
    elif (os.path.exists(clustering_folder + "/output.cluster.log")):
        file = open(clustering_folder + "/output.cluster.log")
    else:
        print("Error: clustering log file not found, try restart: false")
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
    clusters_iterator = sorted(zip(populations, cluster_ids), reverse=True)

    ordered_clusters = []

    # Save sorted cluster ids and populations in dictionaries and populations in list
    for population,index in clusters_iterator:
        ordered_clusters.append({
            'cluster': index,
            'population': population
            })
    
    # Save list in JSON file
    with open(json_file, 'w') as outfile:
        json.dump(ordered_clusters, outfile)

    outfile.close()
    
    # Return list with cluster population and id sorted by population
    return clusters_iterator



def main(args):
    
    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.config_path )
    wdir = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Saving names of the steps included in config.yml
    steps_names = list(global_paths.keys())

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

# STEP 1: Clustering 

    step_dir = wdir + "/" + steps_names[0] + "/"
    
    # Write next action to global log
    global_log.info(steps_names[0] + ": Clustering structures from the trajectory")

    # Action: Initialize GMXCluster and call launch method
    gmx_cluster(**global_paths[steps_names[0]], properties=global_prop[steps_names[0]])

    # Write next action to global log
    global_log.info(steps_names[0] + ": Reading clustering outcome, generating clusters JSON file")

    # Action: save Clusters' id and population in JSON file
    clusters = saveClusters(log_file = steps_names[0] + "_cluster.log", json_file = step_dir + "clusters.json", clustering_folder = step_dir)

    # Write next action to global log
    global_log.info(steps_names[0] + ": Moving clustering output files to folder")

    # Action: move gmx_cluster files into step1_gmx_cluster folder
    moveFileIfItExists(steps_names[0] + "_cluster.log", step_dir + "output.cluster.log")
    moveFileIfItExists(steps_names[0] + "_rmsd-clust.xpm", step_dir + "output.rmsd-clust.xpm")
    moveFileIfItExists(steps_names[0] + "_rmsd-dist.xvg", step_dir + "output.rmsd-dist.xvg")

# STEP 2: Extract model/s

    # Write next action to global log
    global_log.info(steps_names[1] + ": Extracting models of N most populated clusters")

    # Action: Extract the most representative PDB models
    # Properties and paths of step
    prop = global_prop[steps_names[1]]
    paths = global_paths[steps_names[1]]

    # The models that will be extracted are min(number_clusters, models2extract) -> one model per found cluster at most
    modelsToExtract = prop['models_to_extract']
    modelsToExtract = min(modelsToExtract,len(clusters))

    # Default path for output
    default_output = Path(paths['output_structure_path'])

    for i in range(modelsToExtract):

        # Add 'models' keyword and value to properties
        prop.update({'models':[clusters[i][1]]})

        # Modify the output name according to model index
        output_path = addModelSuffix(default_output, i)

        # Update path dictionary
        paths.update({'output_structure_path':output_path})

        # Extract the model 
        extract_model(**paths, properties=prop)
    
# STEP 3: Cavity analysis

    # Write next action to global log
    global_log.info(steps_names[2] + ": Compute protein cavities using fpocket")

    # Action: Search for cavities for each model using fpocket
    # Properties and paths of step
    prop = global_prop[steps_names[2]]
    paths = global_paths[steps_names[2]]

    # From the extracted models use only some, defined in input.yml, cannot be greater than modelsToExtract
    modelsToUse = prop['models_to_use']
    modelsToUse = min(modelsToExtract, modelsToUse)

    # Default paths from input.yml
    default_input = Path(paths['input_pdb_path'])
    default_zip = Path(paths['output_pockets_zip'])
    default_summary = Path(paths['output_summary'])

    for i in range(modelsToUse):

        # Modify the path names according to model index
        input_path = addModelSuffix(default_input, i)
        zip_path = addModelSuffix(default_zip, i)
        summary_path = addModelSuffix(default_summary, i)

        # Update path dictionary
        paths.update({'input_pdb_path':input_path})
        paths.update({'output_pockets_zip':zip_path})
        paths.update({'output_summary':summary_path})

        fpocket_run(**paths, properties=prop)

# STEP 4: Filtering cavities

    # Write next action to global log
    global_log.info(steps_names[3] + ": Filter found cavities")

    # Action: Filter the set of cavities for each model
    # Properties and paths of step
    prop = global_prop[steps_names[3]]
    paths = global_paths[steps_names[3]]

    # Default paths from input.yml
    default_zip_input = Path(paths['input_pockets_zip'])
    default_summary_input = Path(paths['input_summary'])
    default_output = Path(paths['output_filter_pockets_zip'])

    for i in range(modelsToUse):

        # Modify the path names according to model index
        zip_input = addModelSuffix(default_zip_input, i)
        summary_input = addModelSuffix(default_summary_input, i)
        output_path = addModelSuffix(default_output, i)

        # Update path dictionary
        paths.update({'input_pockets_zip':zip_input})
        paths.update({'input_summary':summary_input})
        paths.update({'output_filter_pockets_zip':output_path})

        fpocket_filter(**paths, properties=prop)

    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % args.config_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)")
    
    args = parser.parse_args()

    main(args)

