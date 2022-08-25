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
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def saveClusters(log_file, json_file):
    '''
    Reads the clusters' id and population from log_file, sorts the clusters by 
    population in descending order and writes the sorted list to a JSON file.

    Inputs
        log_file  (str):  log file from clustering (reads from here)
        json_file (str): path for JSON file with sorted cluster ids and populations (writes here) 

    Outputs:
        clusters_iterator (list<tuple>): list with cluster information ordered by population
        (implicitly) JSON file with sorted cluster ids and populations
    '''

    # NOTE: if restart is set to True and the first step was already executed we find an error as the log_file is not found in the working directory

    # Open log file
    file = open(log_file)
    
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

def main(config, system=None):
    
    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(config, system)
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

    # STEP 0: Clustering 

    step_dir = wdir + "/" + steps_names[0] + "/"
    
    # Write next action to global log
    global_log.info(steps_names[0] + ": Clustering structures from the trajectory")

    # NOTE: would the input trajectory be dry in general? Is it necessary?

    # Action: Initialize GMXCluster and call launch method
    # NOTE: perhaps hierarchical clustering is more adequate - study which is the most adequate method
    #       one can compare the numbe of frames with the number of clusters and cluster again until the num clusters << num of frames
    gmx_cluster(**global_paths[steps_names[0]], properties=global_prop[steps_names[0]])

    # Write next action to global log
    global_log.info(steps_names[0] + ": Reading clustering outcome, generating clusters JSON file")

    # Action: save Clusters' id and population in JSON file
    clusters=saveClusters(steps_names[0] + "_cluster.log", step_dir + "clusters.json")

    # Write next action to global log
    global_log.info(steps_names[0] + ": Moving clustering output files to folder")

    # Action: move gmx_cluster files into step1_gmx_cluster folder
    # NOTE: if these files depend on cluster method we'll get an error - or if the first step was executed and have been moved already
    shutil.move(steps_names[0] + "_cluster.log", step_dir + "output.cluster.log")
    shutil.move(steps_names[0] + "_rmsd-clust.xpm", step_dir + "output.rmsd-clust.xpm")
    shutil.move(steps_names[0] + "_rmsd-dist.xvg", step_dir + "output.rmsd-dist.xvg")

    # STEP 1: Extract model/s

    # Write next action to global log
    # NOTE: this N should be defined at config.yml or flag level
    global_log.info(steps_names[1] + ": Extracting models of N most populated clusters")

    # Action: Extract the most representative PDB models
    # Properties and paths of step
    prop = global_prop[steps_names[1]]
    paths = global_paths[steps_names[1]]

    # Consider only the models with highest population, define how many in input.yml file
    models2extract = prop['models_to_extract']

    for i in range(models2extract):

        # Add 'models' keyword and value to properties
        prop.update({'models':[clusters[i][1]]})
        
        # Modify the output name in paths
        model_path = Path(paths['output_structure_path']).parent
        model_path = str(model_path) + "/model" + str(i) + ".pdb" 
        paths.update({'output_structure_path':model_path})

        # Extract the model 
        extract_model(**paths, properties=prop)

    # NOTE: the log file could include here a warning if there are too many clusters or a 
    # representative cluster is missed - RMSD criteria could be included?
    
    # STEP 2: Cavity analysis

    # Write next action to global log
    global_log.info(steps_names[2] + ": Compute protein cavities using fpocket")

    # Action: Search for cavities for each model using fpocket
    # Properties and paths of step
    prop = global_prop[steps_names[2]]
    paths = global_paths[steps_names[2]]

    # From the extracted models use only some, define in input.yml
    models2use = prop['models_to_use']

    # NOTE: check models2use < or = than models2extract
    for i in range(models2use):

        fpocket_run(**paths, properties=prop)

    # STEP 3: Filtering cavities

    # Write next action to global log
    global_log.info(steps_names[3] + ": Filter found cavities")

    # Action: Filter the set of cavities for each model
    # Properties and paths of step
    prop = global_prop[steps_names[3]]
    paths = global_paths[steps_names[3]]

    for i in range(models2use):

        fpocket_filter(**paths, properties=prop)

    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    if system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)

