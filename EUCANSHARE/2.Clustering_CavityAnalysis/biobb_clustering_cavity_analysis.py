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


def addSuffix(original_path, suffix):
    '''
    Adds suffix to original_path before file extension. For example:

    Inputs
    ------

        suffix                (str):  suffix string 
        original_path  (Path class):  path to original file including filename

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
    Reads the centroids' id from the clustering and populations from log_file, sorts the clusters by 
    population in descending order and writes the sorted list to a JSON file.

    Inputs
    -------

        log_file          (str): log file from clustering (reads from here)
        json_file         (str): path for JSON file with sorted cluster ids and populations (writes here) 
        clustering_folder (str): str with name of clustering step folder - just in case Restart 
                                 has been set to True and files have been moved

    Outputs
    -------

        clusters_iterator (list<tuple>): list with tuples containing (population, cluster_ID) sorted by population
        (implicitly) JSON file with sorted cluster ids and populations
    '''

    # Open log file - try to find it in working dir or in step dir
    if os.path.exists(log_file):
        file = open(log_file)
    elif os.path.exists(clustering_folder + "/output.cluster.log"):
        file = open(clustering_folder + "/output.cluster.log")
    else:
        print("Error: clustering log file not found, try 'restart: false'")
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

def joinFiles(source_file_path, destination_file_path):
    '''
    Append the contents of "source_file_path" to "destination_file_path"

    Inputs
    ------

        source_file_path       (str): path to source file
        destination_file_path  (str): path to destination path
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

    return

def addSuffixToPaths(all_paths, suffix, *keywords):
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

        all_paths  (dict): dictionary with paths with those corresponding to "keywords" modified
    '''

    # For all keys passed in keywords, modify path 
    for key in keywords:

        # Original path
        original_path = Path(all_paths[key])

        # Add suffix to path
        newpath = addSuffix(original_path, suffix)

        # Update paths dictionary
        all_paths.update({key : newpath})

    return all_paths

def checkInput(trajectory_path, topology_path, clustering_path, global_log):
    '''
    Checks either a trajectory and a topology path or, alternatively, a clustering path are given.
    If it's not the case, print an error.

    Inputs
    ------

        trajectory_path (str): path to trajectory
        topology_path   (str): path to topology
        clustering_path (str): path to clustering
        global_log   (logger): used to print the error

    Output
    ------

        check (bool): True if inputs OK, False if inputs are missing.
    '''

    check = False

    # If these paths were provided
    if None not in (trajectory_path,topology_path):
        # And they exist
        if os.path.exists(trajectory_path) and os.path.exists(topology_path):
            # Check is True
            check = True

    # If this path was provided
    if clustering_path is not None:
        # And it exists
        if os.path.isdir(clustering_path):
            # And it contains files with the .pdb extension
            if len(glob.glob(os.path.join(clustering_path,"*.pdb"))) > 0:
                # Check is True
                check = True

    # If check is still False
    if check is False:
        # Provide warning to user and exit
        global_log.info("ERROR: Input check failed. Make sure the following are suitable trajectory and topology paths or a suitable clustering path with pdb files in it.")
        global_log.info("       trajectory_path: {}".format(trajectory_path))
        global_log.info("       topology_path: {}".format(topology_path))
        global_log.info("       clustering_path: {}".format(clustering_path))

    return check

def createSummary(pockets_path, default_summary_path, models_population, global_log):
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
        global_summary (dict): dictionary with models and pockets. See example:
        
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
    logList = pockets_path.rglob(logNamePattern)

    # Iterate through all log files -> find model name and available pockets if any
    for log in logList:

        # Find model name NOTE: hardcoded name of step and file, should be changed
        modelName = findMatchingLine(pattern=r'step3_cavity_analysis/all_pockets_(\S+).zip', filepath=str(log))

        # Find pockets  NOTE: is there a more robust way of finding models and pockets after filtering? 
        pocket_lines = findMatchingLines(pattern=r'pocket\d+$', filepath=str(log))

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
            summary_path = addSuffix(original_path=Path(default_summary_path), suffix=modelName)

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

def findMatchingLine(pattern, filepath):
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

def findMatchingLines(pattern, filepath):
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


def main_wf(configuration_path, trajectory_path = None, topology_path = None, clustering_path = None, last_step = None):
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and analyzes the cavities of the most representative
    structures. Then filters the cavities according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        configuration_path (str): path to input.yml 
        trajectory_path    (str): path to the trajectory that will be clustered, alternatively one can give the clustering path. (accepted formats: xtc, trr, cpt, gro, g96, pdb or tng)
        topology_path      (str): path to the topology of the trajectory in trajectory_path (accepted formats: tpr, gro, g96, pdb or brk)
        clustering_path    (str): path to the folder with the most representative structures in pdb format from an external clustering (*.pdb)
        last_step          (str): last step of the workflow to execute ('ndx', 'cluster', 'cavity', 'all')

    Outputs
    -------

        /output folder
        working_dir_path (str): path to working directory 
        clusters_path    (str): path to clustering results (most representative clusters)
        pockets_path     (str): path to filtered pockets (pockets that passed the filter)
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

    # Set default value for 'last_step' arg
    if last_step is None:
        last_step = 'all'

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

    # Check input files are given: either trajectory + topology or clustering has to be given
    if not checkInput(trajectory_path, topology_path, clustering_path, global_log):
        return None, None

    # If clustering is not given externally -> cluster the input trajectory and return cluster path
    if clustering_path is None:
    
    # STEP 0: Fit Index file creation

        # Properties and paths of step
        prop_selFit = global_prop['step0_gmx_select_fit']
        paths_selFit= global_paths['step0_gmx_select_fit']

        # If topology is given through CLI -> prioritize over input.yml
        if topology_path is not None:
            paths_selFit.update({'input_structure_path' : topology_path})
        
        # Write next action to global log
        global_log.info("step0_gmx_select_fit: Creation of index file for clustering 'fit_selection'")

        # Create index file with 'FitGroup'
        gmxselect(**paths_selFit, properties=prop_selFit)

        # Validate step
        if not validateStep(paths_selFit["output_ndx_path"]):
            global_log.info("    ERROR: No index file was created. Check gmx_select_fit atom selection string")
            return None, None

    # STEP 0: Output Index file creation

        # Properties and paths of step
        prop_selOut = global_prop['step0_gmx_select_output']
        paths_selOut = global_paths['step0_gmx_select_output']

        # If topology is given through CLI -> prioritize over input.yml
        if topology_path is not None:
            paths_selOut.update({'input_structure_path' : topology_path})
        
        # Write next action to global log
        global_log.info("step0_gmx_select_output: Creation of index file for clustering 'output_selection'")

        # Create index file with 'OutputGroup'
        gmxselect(**paths_selOut, properties=prop_selOut)

        # Validate step
        if not validateStep(paths_selOut["output_ndx_path"]):
            global_log.info("    ERROR: No index file was created. Check gmx_select_output atom selection string")
            return None, None

        # Add FitGroup to index_Fit_Output.ndx 
        joinFiles(source_file_path = paths_selFit["output_ndx_path"], destination_file_path = paths_selOut["output_ndx_path"])

        # Check if this should be the final step
        if last_step == 'ndx':
            global_log.info("Index files created.")
            return None, None

    # STEP 1: Clustering with gmx_cluster

        step_dir = os.path.join(conf.get_working_dir_path(), "step1_gmx_cluster")
        
        # Write next action to global log
        global_log.info("step1_gmx_cluster: Clustering structures from the trajectory")

        # Properties and paths of step
        prop_clust = global_prop["step1_gmx_cluster"]
        paths_clust = global_paths["step1_gmx_cluster"]

        # If trajectory path is given through CLI -> prioritize over input.yml
        if trajectory_path is not None:
            paths_clust.update({'input_traj_path' : trajectory_path})
        
        # If topology is given through CLI -> prioritize over input.yml
        if topology_path is not None:
            paths_clust.update({'input_structure_path' : topology_path})

        # NOTE: The clustering is done using all the atoms for the RMSD calculation. Is it possible to use a sub-set with gmx cluster? Maybe including ttclust project? or mdtraj?
        # Using gromos method the RMSD calculation is done with respect to a certain group... But it also aligns wrt this group... 
        # The idea would be to align with respect to one group and compute rmsd with respect to another, this cannot be easily implemented for now -> external clustering can be provided

        # Cluster trajectory
        gmx_cluster(**paths_clust, properties=prop_clust)

        # Save path of clustering results
        clusters_path = prop_clust["path"]

        # Validate step
        if not validateStep(paths_clust["output_pdb_path"]):
            global_log.info("    ERROR: No PDB file was created.")
            return clusters_path, None

        # Write next action to global log
        global_log.info( "step1_gmx_cluster: Reading clustering outcome, generating clusters JSON file")

        # Save Model's (centroids of clusters) ID and population in JSON file
        models_population = saveClusters(log_file = "step1_gmx_cluster_cluster.log", 
                                json_file = os.path.join(step_dir, "clusters.json"), 
                                clustering_folder = step_dir)

        # Write next action to global log
        global_log.info("step1_gmx_cluster: Moving clustering output files to folder")

        # Action: move gmx_cluster files into step1_gmx_cluster folder NOTE: check if this is still needed
        moveFileIfItExists("step1_gmx_cluster_cluster.log", os.path.join(step_dir, "output.cluster.log"))
        moveFileIfItExists("step1_gmx_cluster_rmsd-clust.xpm", os.path.join(step_dir, "output.rmsd-clust.xpm"))
        moveFileIfItExists("step1_gmx_cluster_rmsd-dist.xvg", os.path.join(step_dir, "output.rmsd-dist.xvg"))

    # STEP 2: Extract the most representative PDB models

        # Write next action to global log
        global_log.info("step2_extract_models: Extracting models of N most populated clusters")

        # Properties and paths of step
        prop_models = global_prop["step2_extract_models"]
        paths_models = global_paths["step2_extract_models"]

        # The models that will be extracted are: min(models2extract, number_clusters)
        modelsToExtract = min(prop_models['models_to_extract'], len(models_population))

        # Extract the most populated models from the pdb with all centroids
        for i in range(modelsToExtract):

            # Copy original paths
            paths_models = global_paths["step2_extract_models"].copy()

            # Modify the path names according to model index
            addSuffixToPaths(paths_models, str(i),
                            "output_structure_path")

            # Update 'models' property
            prop_models.update({'models':[models_population[i][1]]})

            # Extract the model 
            extract_model(**paths_models, properties=prop_models)

            # Validate step
            if not validateStep(paths_models["output_structure_path"]):
                global_log.info("    ERROR: No model PDB was extracted.")
                return clusters_path, None
        
        # Obtain the full sorted list of pdb files from previous step path
        pdb_paths = sorted(glob.glob(os.path.join(prop_models["path"],"*.pdb")))

        # Save cluster path
        clusters_path = prop_models["path"]

        # Check if this should be the final step
        if last_step == 'cluster':
            global_log.info("Clustering completed.")
            return clusters_path, None

    # If clustering is given externally 
    else:

        # Path with clustering was given
        clusters_path = clustering_path

        # Obtain the full sorted list of pdb files from given path
        pdb_paths = sorted(glob.glob(os.path.join(clustering_path,"*.pdb")))

        # Population information will not be available NOTE: it would be nice to give this as well as 
        # part of the external clustering results
        models_population = None

# STEP 3: Cavity analysis

    # Write next action to global log
    global_log.info("step3_cavity_analysis: Compute protein cavities using fpocket")

    # Properties of step
    prop_fpocket = global_prop["step3_cavity_analysis"]

    # Number of models to use
    modelsToUse = min(prop_fpocket['models_to_use'], len(pdb_paths))

    # Keep only the ones that will be used
    pdb_paths = pdb_paths[:modelsToUse]

    # Analyze cavities of each pdb 
    for pdb_path in pdb_paths:

        # Copy original paths
        paths_fpocket = global_paths["step3_cavity_analysis"].copy()

        # Modify the path names according to pdb name
        addSuffixToPaths(paths_fpocket,
                         Path(pdb_path).stem,
                        "output_pockets_zip",
                        "output_summary")

        # Update input pdb path
        paths_fpocket.update({'input_pdb_path': pdb_path})

        # Run cavity analysis
        fpocket_run(**paths_fpocket, properties=prop_fpocket)

        # Validate step
        if not validateStep(paths_fpocket["output_pockets_zip"], paths_fpocket["output_summary"]):
            global_log.info("    ERROR: no output from fpocket.")
            return None, None

    pockets_path = prop_fpocket["path"]

    # Check if this should be the final step
    if last_step == 'cavity':
        global_log.info("Cavity analysis completed.")
        return clusters_path, pockets_path

# STEP 4: Filtering cavities

    # Write next action to global log
    global_log.info("step4_filter_cavities: Filter found cavities")

    # Properties and paths of step
    prop_filter = global_prop["step4_filter_cavities"]
    paths_filter = global_paths["step4_filter_cavities"]

    # Filter the cavities found for each pdb
    for pdb_path in pdb_paths:

        # Copy original paths
        paths_filter = global_paths["step4_filter_cavities"].copy()

        # Modify the path names according to pdb name
        addSuffixToPaths(paths_filter,
                         Path(pdb_path).stem,
                        "input_pockets_zip",
                        "input_summary",
                        "output_filter_pockets_zip")

        # Filter cavities
        fpocket_filter(**paths_filter, properties=prop_filter)
    
    # Find path to all filtered pockets
    pockets_path = prop_filter["path"]

    # Create and print available models with their pockets and populations after filtering
    global_summary = createSummary(pockets_path = pockets_path, 
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

    return conf.get_working_dir_path(), clusters_path, pockets_path, global_summary

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)", 
                        required=True)

    parser.add_argument('--input-traj', dest='input_traj',
                        help="Input trajectory path, prioritized over input.yml paths (xtc, trr, cpt, gro, g96, pdb or tng format)", 
                        required=False)

    parser.add_argument('--input-top', dest='input_topology',
                        help="Input topology/structure path, prioritized over input.yml paths (tpr, gro, g96, pdb or brk format)", 
                        required=False)

    # NOTE: add option to give external clustering result to wf - as this clustering has limitations
    parser.add_argument('--input-clust', dest='clust_path',
                        help="Input path to representative structures (pdb files, all will be used)", 
                        required=False)

    # Execute workflow until 'last_step' -> all executes all steps (all is default)
    parser.add_argument('--until', dest='last_step', 
                        help="(Opt, default: all) Extent of the pipeline to execute (ndx, cluster, cavity, all)", 
                        required=False)

    args = parser.parse_args()

_,_,_,_ = main_wf(configuration_path = args.config_path, 
                    trajectory_path = args.input_traj, 
                    topology_path = args.input_topology, 
                    clustering_path = args.clust_path, 
                    last_step = args.last_step)

