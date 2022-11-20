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
from pathlib import Path

from pkg_resources import resource_string
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_gromacs.gromacs.make_ndx import make_ndx
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
    ------

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
    ------

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
    -------

        log_file          (str): log file from clustering (reads from here)
        json_file         (str): path for JSON file with sorted cluster ids and populations (writes here) 
        clustering_folder (str): str with name of clustering step folder - just in case Restart 
                                 has been set to True and files have been moved

    Outputs
    -------

        clusters_iterator (list<tuple>): list with cluster information ordered by population
        (implicitly) JSON file with sorted cluster ids and populations
    '''

    # Open log file - try to find it in working dir or in step dir
    if os.path.exists(log_file):
        file = open(log_file)
    elif os.path.exists(clustering_folder + "/output.cluster.log"):
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

def validateStep(*output_paths):
    '''
    Check all output files exist and are not empty
    
    Inputs
    ------
        *output_paths (str): variable number of paths to output file/s

    Output
    ------
        validation_result (bool): result of validation
    '''

    # Initialize value 
    validation_result = True

    # Check existence of files
    for path in output_paths:
        validation_result = validation_result and os.path.exists(path)

    # Check files are not empty if they exist
    if (validation_result):

        for path in output_paths:
            file_not_empty = os.stat(path).st_size > 0
            validation_result = validation_result and file_not_empty

    return validation_result

def main(args):
    
    start_time = time.time()

    # Set default value for 'to_do' arg
    if args.to_do is None:
        args.to_do = 'all'

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.config_path )
    wdir = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=wdir, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Saving names of the steps included in config.yml
    step_names = list(global_paths.keys())

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

# STEP 0: Index file creation

    # Properties and paths of step
    prop_ndx = global_prop[step_names[0]]
    paths_ndx = global_paths[step_names[0]]
    
    # Only if ndx file is not given already with '--ndx-file'
    if args.ndx_path is None: 

        # Write next action to global log
        global_log.info(step_names[0] + ": Creation of index file for clustering" + "\n")

        # Action: Initialize MakeNdx and call launch method
        make_ndx(**paths_ndx, properties=prop_ndx)

        # Validate step
        if not validateStep(paths_ndx["output_ndx_path"]):
            global_log.info("    ERROR: No index file was created. Check atom selection string in input file")
            return 0

    # Check if this should be the final step
    if args.to_do == 'ndx':
        global_log.info("Index file created.")
        return 0

# STEP 1: Clustering 

    step_dir = os.path.join(wdir, step_names[1])
    
    # Write next action to global log
    global_log.info(step_names[1] + ": Clustering structures from the trajectory" + "\n")

    # Properties and paths of step
    prop_clust = global_prop[step_names[1]]
    paths_clust = global_paths[step_names[1]]

    if args.ndx_path is None:
        # If index path was not provided, it was created and paths are correct
        global_log.info("   Using created index file")
    else:
        # If index path was provided, change the corresponding paths
        global_log.info("   Using external index file")
        paths_clust.update({"input_index_path" : args.ndx_path})
        paths_clust.update({"input_index_path" : paths_ndx["input_structure_path"]})

    # Action: Initialize GMXCluster and call launch method
    gmx_cluster(**paths_clust, properties=prop_clust)

    # NOTE: The clustering is always done using all the atoms for the RMSD calculation? Is it possible to use a sub-set with gmx cluster? Maybe including ttclust project? or mdtraj?
    # Using gromos the RMSD calculation is done with respect to a certain group... But we have to align previously? the idea is to align with respect to one group and compute rmsd with respect to another.
    
    # Validate step
    if not validateStep(paths_clust["output_pdb_path"]):
        global_log.info("    ERROR: No PDB file was created.")
        return 0

    # Write next action to global log
    global_log.info(step_names[1] + ": Reading clustering outcome, generating clusters JSON file")

    # Action: save Clusters' id and population in JSON file
    clusters = saveClusters(log_file = step_names[1] + "_cluster.log", json_file = os.path.join(step_dir, "clusters.json"), clustering_folder = step_dir)

    # Write next action to global log
    global_log.info(step_names[1] + ": Moving clustering output files to folder")

    # Action: move gmx_cluster files into step1_gmx_cluster folder
    moveFileIfItExists(step_names[1] + "_cluster.log", os.path.join(step_dir, "output.cluster.log"))
    moveFileIfItExists(step_names[1] + "_rmsd-clust.xpm", os.path.join(step_dir, "output.rmsd-clust.xpm"))
    moveFileIfItExists(step_names[1] + "_rmsd-dist.xvg", os.path.join(step_dir, "output.rmsd-dist.xvg"))

# STEP 2: Extract model/s

    # Write next action to global log
    global_log.info(step_names[2] + ": Extracting models of N most populated clusters")

    # Action: Extract the most representative PDB models
    # Properties and paths of step
    prop = global_prop[step_names[2]]
    paths = global_paths[step_names[2]]

    # The models that will be extracted are min(number_clusters, models2extract) -> one model per found cluster at most
    modelsToExtract = prop['models_to_extract']
    modelsToExtract = min(modelsToExtract,len(clusters))

    # Default path for output
    default_output = Path(paths['output_structure_path'])

    for i in range(modelsToExtract):

        # Update 'models' index
        prop.update({'models':[clusters[i][1]]})

        # Modify the output name - according to model index
        output_path = addModelSuffix(default_output, i)

        # Update paths dictionary
        paths.update({'output_structure_path':output_path})

        # Extract the model 
        extract_model(**paths, properties=prop)

        # Validate step
        if not validateStep(paths["output_structure_path"]):
            global_log.info("    ERROR: No model PDB was extracted.")
            return 0
    
    # Check if this should be the final step
    if args.to_do == 'cluster':
        global_log.info("Clustering completed.")
        return 0

# STEP 3: Cavity analysis

    # Write next action to global log
    global_log.info(step_names[3] + ": Compute protein cavities using fpocket")

    # Action: Search for cavities for each model using fpocket
    # Properties and paths of step
    prop = global_prop[step_names[3]]
    paths = global_paths[step_names[3]]

    # From the extracted models use only some, defined in input.yml, cannot be greater than modelsToExtract
    modelsToUse = prop['models_to_use']
    modelsToUse = min(modelsToExtract, modelsToUse)

    # Default paths from input.yml
    default_input = Path(paths['input_pdb_path'])
    default_zip = Path(paths['output_pockets_zip'])
    default_summary = Path(paths['output_summary'])

    for i in range(modelsToUse):

        # Modify the path names - according to model index
        input_path = addModelSuffix(default_input, i)
        zip_path = addModelSuffix(default_zip, i)
        summary_path = addModelSuffix(default_summary, i)

        # Update path dictionary
        paths.update({'input_pdb_path':input_path})
        paths.update({'output_pockets_zip':zip_path})
        paths.update({'output_summary':summary_path})

        fpocket_run(**paths, properties=prop)

        # Validate step
        if not validateStep(paths["output_pockets_zip"], paths["output_summary"]):
            global_log.info("    ERROR: no output from fpocket.")
            return 0

    # Check if this should be the final step
    if args.to_do == 'cavity':
        global_log.info("Cavity analysis completed.")
        return 0

# STEP 4: Filtering cavities

    # Write next action to global log
    global_log.info(step_names[4] + ": Filter found cavities")

    # Action: Filter the set of cavities for each model
    # Properties and paths of step
    prop = global_prop[step_names[4]]
    paths = global_paths[step_names[4]]

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
    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)
    
    parser.add_argument('--ndx-file', dest='ndx_path',
                        help="(Opt, default: None) External index file path to be used. Select groups for 'fit_selection' and 'output_selection' contained in index file.", 
                        required=False)

    # Execute workflow until 'to_do' step -> all executes all steps (all is default)
    parser.add_argument('--until', dest='to_do', 
                        help="(Opt, default: all) Extent of the pipeline to execute (ndx, cluster, cavity, all)", 
                        required=False)

    args = parser.parse_args()

    main(args)

