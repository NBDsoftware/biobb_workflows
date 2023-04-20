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
import json
import shutil
from pathlib import Path

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_gromacs.gromacs.gmxselect import gmxselect
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter


def add_suffix(original_path, suffix):
    '''
    Adds suffix to original_path before file extension. For example:

    Inputs
    ------

        original_path  (Path class):  path to original file including filename
        suffix                (str):  suffix string 

    Output
    ------

        new_path   (str): new path

    original_path = Path('/home/user/ClusteringWF/output/step2_extract_models/output.pdb')
    suffix = '0'

    new_path = '/home/user/ClusteringWF/output/step2_extract_models/output_0.pdb'
    '''
    # Identify original file name, extension and path
    filename = original_path.stem
    fileExtension = original_path.suffix
    directory = original_path.parent

    # Create new file name
    new_fileName = filename + '_' + suffix + fileExtension

    # Create new path
    new_path = os.path.join(str(directory), new_fileName)

    return new_path

def add_suffix_to_paths(all_paths, suffix, *keywords):
    '''
    Goes over paths corresponding to keywords in all_paths dictionary, adds 
    suffix to each filename.
    
    Inputs
    ------

        all_paths  (dict): dictionary with old paths of step
        suffix      (str): string corresponding to suffix
        keywords    (str): keywords of output/input paths that will be modified

                    /path/to/file/filename.txt  ---> /path/to/file/filename_suffix.txt

    Output
    ------


        (implicit) all_paths  (dict): dictionary with paths corresponding to "keywords" modified
    '''

    # For all keys passed in keywords, modify path 
    for key in keywords:

        # Original path
        original_path = Path(all_paths[key])

        # Add suffix to path
        newpath = add_suffix(original_path, suffix)

        # Update paths dictionary
        all_paths.update({key : newpath})

    return

def move_files(dest_src: str, *origin_src):
    '''
    Checks the existence of all paths in origin_src and move them to 'dest_src' if they exist

    Inputs
    ------

        dest_src    :  path to destination folder
        *origin_src  :  paths to original file including filename
    '''

    # For each path in origin_src
    for path in origin_src:
        # If path exists
        if os.path.exists(path):
            # Create new path
            new_path = os.path.join(dest_src, os.path.basename(path))
            # Move path to new path
            shutil.move(path, new_path)

def get_clusters_population(log_name: str, clustering_folder: str, global_log) -> list:
    '''
    Reads the centroids' ID and populations from the clustering log, sorts the clusters by 
    population in descending order and writes the sorted list to a JSON file in the clustering_folder.

    Inputs
    -------

        log_name          : log file name from clustering 
        clustering_folder : path to clustering step folder

    Outputs
    -------

        /json_file  with sorted cluster ids and populations

        clusters_population (list<tuple>): list with tuples containing (population, cluster_ID) sorted by population
    '''
    # NOTE: this is needed because gmx_cluster doesn't include all the required output paths in the constructor and dictionaries
    # Try to find log in working dir
    if os.path.exists(log_name):
        file = open(log_name)
    # Try to find log in clustering folder
    elif os.path.exists(os.path.join(clustering_folder, log_name)):
        file = open(os.path.join(clustering_folder, log_name))
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
    with open(os.path.join(clustering_folder, "clusters.json"), 'w') as outfile:
        json.dump(ordered_clusters, outfile)

    outfile.close()
    
    # Return list with cluster population and id sorted by population
    return clusters_population

def merge_files(source_file_path: str, destination_file_path: str) -> None:
    '''
    Append the contents of "source_file_path" to "destination_file_path"

    Inputs
    ------

        source_file_path       : path to source file
        destination_file_path  : path to destination path
    '''

    # Open the source file in read mode 
    source_file = open(source_file_path, 'r')

    # Open destination file in append mode
    destination_file = open(destination_file_path, 'a+')

    # Append contents of source to destination
    destination_file.write(source_file.read())

    # Close both
    source_file.close()
    destination_file.close()

def create_summary(pockets_path, default_summary_path, models_population, global_log):
    '''
    Print in log file all available pockets for each model found in the input folder for the pocket selection step. 
    Include also information for each pocket and population for each model.
    
    Inputs
    ------

        pockets_path         (str): Path to pockets folder 
        default_summary_path (str): Default path to cavity summary with information for each pocket
        models_population   (list): List with tuples containing (population, cluster_ID) ordered by population
        global_log        (Logger): Object of class Logger (dumps info to global log file of this workflow)
    
    Output
    ------

        Info dumped to log file

        global_summary (dict):  dictionary with models and pockets ordered by modelname. See example:
        
            {
                modelname1 : {        
                    population : 143, 
                    pockets : [1 , 2], 
                    pocket1 : {info pocket 1}, 
                    pocket2 : {info pocket 2}
                    }, 
                modelname2 : {
                    ...
                } 
            }
    '''
    
    global_log.info("    Available models and pockets after filtering:")

    # Pattern that can be used to extract ID from string with pocket name: (str) pocket3 -> (int) 3
    pocketID_pattern = r'pocket(\d+)'

    # Pattern that can be used to extract model ID from string with model name
    # This ID reflects the ordering of the population, not the original ID of the clustering step
    modelID_pattern = r'cluster_(\d+)'

    # Dictionary where all available model_summary dictionaries will be saved
    global_summary = {}

    # Convert string to Path class
    pockets_path = Path(pockets_path)

    # Pattern for pocket filtering log file names
    logNamePattern = pockets_path.name + "_log*.out"

    # Find all log files matching pattern
    logList = sorted(pockets_path.rglob(logNamePattern))

    # Iterate through all log files -> find model name and available pockets if any
    for log in logList:

        # Find model name NOTE: hardcoded name of step and file, should be changed
        modelName = find_matching_line(pattern=r'step3_cavity_analysis/all_pockets_(\S+).zip', filepath=str(log))

        # Find pockets  NOTE: is there a more robust way of finding models and pockets after filtering? 
        pocket_lines = find_matching_lines(pattern=r'pocket\d+$', filepath=str(log))

        # If some model and pocket/s are found in log
        if None not in (modelName, pocket_lines):

            # Dictionary where pockets, model population and information of each pocket will be saved
            model_summary = {}

            # If clustering was done externally we might not have this information 
            if models_population is not None:

                # Search for the model ID in model name
                match = re.search(modelID_pattern, modelName)
                model_ID = int(match[1])

                # Save population of model (model = centroid of cluster)
                model_summary.update({'population' : models_population[model_ID][0]})

            # Path to this model's summary with information for all pocket found
            summary_path = add_suffix(original_path=Path(default_summary_path), suffix=modelName)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # list with pocket IDs for this model
            pocket_IDs = []

            # For each filtered pocket: pocket1, pocket4 ...
            for pocket_line in pocket_lines:

                # Strip whitespace if any
                pocket_line = pocket_line.strip()

                # Save entry with pocket information
                model_summary.update({pocket_line : all_pockets_summary[pocket_line]})

                # Search for the pocket ID in pocket name
                match = re.search(pocketID_pattern, pocket_line)

                # Append pocket ID to list
                pocket_IDs.append(int(match[1]))

            # Save pocket IDs
            model_summary.update({'pockets' : pocket_IDs})

            # Update global_summary
            global_summary.update({modelName : model_summary})

    # Print in log file with readable format
    pretty_summary = json.dumps(global_summary, indent=2)

    global_log.info(pretty_summary)

    return global_summary

def find_matching_line(pattern, filepath):
    '''
    Finds first line in file containing a given pattern

    Inputs
    ------

        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------
    
        line (str): line matching the pattern or None if no line matches the pattern
    '''

    # Open log file
    file = open(filepath, 'r')
    
    # Read all lines
    lines = file.readlines()

    # Search matching string in each line
    for line in lines:

        match = re.search(pattern, line)

        if match is not None:

            # Close file
            file.close()

            return match[1]
    
    file.close()

    return None

def find_matching_lines(pattern, filepath):
    '''
    Finds all lines in file containing a given pattern

    Inputs
    ------

        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------

        lines (list(str)): lines matching the pattern or None if no line matches the pattern
    '''

    # Open log file
    file = open(filepath, 'r')
    
    # Read all lines
    lines = file.readlines()

    # List to save matching lines
    matchingLines = []

    # Search matching string in each line
    for line in lines:

        match = re.search(pattern, line)

        if match is not None:

            matchingLines.append(match[0])

    # Close file
    file.close()

    # If no lines contained the pattern, return None
    if len(matchingLines) == 0:
        matchingLines = None
    
    return matchingLines


def main_wf(configuration_path, clustering_path = None):
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and analyzes the cavities of the most representative
    structures. Then filters the cavities according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        configuration_path (str): path to input.yml 
        clustering_path    (str): (Optional) path to the folder with the most representative structures in pdb format from an external clustering (*.pdb)

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties

        global_summary  (dict): dictionary with models and pockets. See example:
        
            {
                modelname1 : {
                    population : 143,  
                    pockets: [1, 2],
                    pocket1 : {info pocket 1}, 
                    pocket2 : {info pocket 2}
                    }, 
                modelname2 : {
                    ...
                } 
            }
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # If clustering is not given externally -> cluster the input trajectory and return cluster path
    if clustering_path is None:
    
        # STEP 0: Create index file with FitGroup for gmx_cluster
        global_log.info("step0_gmx_select_fit: Creation of index file for clustering 'fit_selection'")
        gmxselect(**global_paths['step0_gmx_select_fit'], properties=global_prop['step0_gmx_select_fit'])

        # STEP 0: Create index file with OutputGroup for gmx_cluster
        global_log.info("step0_gmx_select_output: Creation of index file for clustering 'output_selection'")
        gmxselect(**global_paths['step0_gmx_select_output'], properties=global_prop['step0_gmx_select_output'])

        # Merge both index files 
        merge_files(source_file_path = global_paths['step0_gmx_select_fit']["output_ndx_path"], 
                    destination_file_path = global_paths['step0_gmx_select_output']["output_ndx_path"])

        # STEP 1: Clustering trajectory with gmx_cluster
        global_log.info("step1_gmx_cluster: Clustering structures from the trajectory")
        gmx_cluster(**global_paths["step1_gmx_cluster"], properties=global_prop["step1_gmx_cluster"])

        # Save Model's (centroids of clusters) ID and population in JSON file
        global_log.info( "step1_gmx_cluster: Reading clustering outcome, generating clusters JSON file")
        models_population = get_clusters_population(log_name = "step1_gmx_cluster_cluster.log",
                                                    clustering_folder = global_prop["step1_gmx_cluster"]['path'],
                                                    global_log = global_log)

        # NOTE: this is needed because gmx_cluster doesn't include all the output paths in the constructor and dictionaries - thus some are left behind and not copied inside the step folder
        move_files(global_prop["step1_gmx_cluster"]['path'], "step1_gmx_cluster_cluster.log", "step1_gmx_cluster_rmsd-clust.xpm", "step1_gmx_cluster_rmsd-dist.xvg")

        # Number of models to extract
        num_models_to_extract = min(prop_models['models_to_extract'], len(models_population))

        # STEP 2: Extract the most representative PDB models
        prop_models = global_prop["step2_extract_models"]
        global_log.info(f"step2_extract_models: Extracting models of the {num_models_to_extract} most populated clusters")

        # Extract the most populated models from the pdb with all centroids
        for i in range(num_models_to_extract):

            # Copy original paths - to avoid accumulation of suffixes
            paths_models = global_paths["step2_extract_models"].copy()

            # Modify the path names according to model index
            add_suffix_to_paths(paths_models, str(i), "output_structure_path")

            # Update 'models' property (ID)
            prop_models.update({'models':[models_population[i][1]]})

            # Extract the model 
            extract_model(**paths_models, properties=prop_models)
        
        # Obtain the full sorted list of pdb files from previous step path
        pdb_paths = sorted(glob.glob(os.path.join(prop_models["path"],"*.pdb")))

    # If clustering is given externally 
    else:

        # Obtain the full sorted list of pdb files from clustering path
        pdb_paths = sorted(glob.glob(os.path.join(clustering_path,"*.pdb")))

        # Population information will not be available in this case
        models_population = None

    # STEP 3: Cavity analysis
    global_log.info("step3_cavity_analysis: Compute protein cavities using fpocket")

    # Number of models to use
    models_to_use = min(global_prop["step3_cavity_analysis"]['models_to_use'], len(pdb_paths))

    # Analyze cavities of each pdb 
    for pdb_path in pdb_paths[:models_to_use]:

        # Copy original paths
        paths_fpocket = global_paths["step3_cavity_analysis"].copy()

        # Modify the path names according to pdb name
        add_suffix_to_paths(paths_fpocket, Path(pdb_path).stem, "output_pockets_zip", "output_summary")

        # Update input pdb path
        paths_fpocket.update({'input_pdb_path': pdb_path})

        # Run cavity analysis
        fpocket_run(**paths_fpocket, properties=global_prop["step3_cavity_analysis"])

    # STEP 4: Filtering cavities
    global_log.info("step4_filter_cavities: Filter found cavities")

    # Filter the cavities found for each pdb
    for pdb_path in pdb_paths[:models_to_use]:

        # Copy original paths
        paths_filter = global_paths["step4_filter_cavities"].copy()

        # Modify the path names according to pdb name
        add_suffix_to_paths(paths_filter, Path(pdb_path).stem, "input_pockets_zip", "input_summary", "output_filter_pockets_zip")

        # Filter cavities
        fpocket_filter(**paths_filter, properties=global_prop["step4_filter_cavities"])

    # Create and print available models with their pockets and populations after filtering
    global_summary = create_summary(pockets_path = global_prop["step4_filter_cavities"]["path"], 
                                    default_summary_path = global_paths['step3_cavity_analysis']['output_summary'],
                                    models_population = models_population, 
                                    global_log = global_log)

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

    return global_paths, global_prop, global_summary

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)

    parser.add_argument('--input-clust', dest='clust_path',
                        help="Input path to representative structures (pdb files, all will be used)", 
                        required=False)

    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            clustering_path = args.clust_path)

