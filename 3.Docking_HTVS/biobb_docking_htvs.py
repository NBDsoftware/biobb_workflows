#!/usr/bin/env python3

# Importing all the needed libraries
import os
import re
import time
import argparse
import itertools 
from pathlib import Path
import multiprocessing as mp
from biobb_io.api.drugbank import drugbank
from biobb_io.api.ideal_sdf import ideal_sdf
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_vs.utils.box import box
from biobb_chemistry.babelm.babel_convert import babel_convert
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens
from biobb_vs.vina.autodock_vina_run import autodock_vina_run

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

def findMatchingStr(pattern, filepath):
    '''
    Finds str in file corresponding to a given pattern

    Inputs
    ------
        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output
    ------
        match_str (str): string matching the pattern or None if no piece matches the pattern
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

            return match.group(1)

    file.close()

    return None

def findNumberInString(string):
    '''
    Finds and returns first integer number in string. If no integer is found a 0 is returned.

    Inputs
    ------
        string (str): string to search int in
    
    Output
    ------
        number (int): integer found in string or 0 if no integer is found.
    '''

    # Pattern corresponding to a digit of one or more characters
    digitPattern = r'\d+'

    # Search for the first match in string
    number = re.search(digitPattern, string)
    
    if (number == None):
        # If there is no match, return a 0
        number = 0
    else:
        # If there is a match, convert match object to string
        number = number[0]

    return int(number)

def printAvailablePockets(pockets_path, global_log):
    '''
    Print in log file all available pockets for each model found in the input folder for the pocket selection step
    
    Inputs
    ------
        pockets_path  (str): Path to pockets folder including file name
        global_log (Logger): Object of class Logger (dumps info to global log file of this workflow)

    Output
    ------
        Info dumped to log file
    '''
    
    # Convert string to Path class
    pocketsLogPath = Path(pockets_path)

    # Path of pocket filtering step
    stepPath = pocketsLogPath.parent

    # Pattern for pocket filtering log files names
    logNamePattern = stepPath.name + "_log*.out"

    # Find all log files matching pattern
    logList = stepPath.rglob(logNamePattern)

    # Iterate through all log files -> print the centroid name + available pockets in the global log
    for log in logList:

        # Find pockets 
        pockets = findMatchingLines(pattern=r'pocket\d+$', filepath=str(log))

        # Find name of centroid
        centroidName = findMatchingStr(pattern=r'/all_pockets_(\S+).zip', filepath=str(log))

        # Print to global log
        global_log.info("   Model: {}".format(centroidName))

        for pocket in pockets:
            global_log.info("        {}".format(pocket))

    return

def findAffinityInLog(log_path):
    '''
    Find best binding affinity among different poses
    
    Inputs
    ------
        log_path (str): paths of log file with results

    Outputs
    -------
        affinity (float): best affinity among poses  
    '''

    # Find affinity of best pose from log
    affinity = float(findMatchingStr(pattern=r'\s+1\s+(\S+)\s+', filepath=log_path))

    return affinity

def printTopLigands(affinities_list, properties, global_log):
    '''
    Print top ligands 

    Inputs
    ------
        affinities_list  (list): list with tuples -> (affinity, ligand identifier)
        properties       (dict): properties of step 8
        global_log     (logger): global log
    
    Output
    ------

        affinities_list        (list): list with top ligand tuples -> (affinity, ligand identifier) ordered by affinity
    '''
    # NOTE: perhaps print csv with data?

    # Sort list according to affinity
    affinities_list = sorted(affinities_list)

    # Find number of ligands to print in global log
    numLigandsToPrint = min(properties['number_top_ligands'], len(affinities_list)) 

    # Exclude the rest
    affinities_list = affinities_list[:numLigandsToPrint]

    # Print best ligands and their affinities in log file
    for item in affinities_list:

        ligand_affinity, ligand_identifier = item

        global_log.info("    Affinity: {} for ligand {}".format(ligand_affinity, ligand_identifier))

    return affinities_list

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
    
    # Erase files if they are empty
    if (validation_result == False):
        removeFiles(*output_paths)

    return validation_result

def readLigandLine(ligand_line):
    '''
    Read line from ligand library. Accepts the following formats as valid:

    Format 1:

    "ligand1_id"

    Format 2:

    "ligand1_id  name_ligand1"

    Where ligand_id is either a PDB code, Drug Bank ID or SMILES and name_ligand is a string with the ligand name
    
    Inputs
    ------
        ligand_line (str): the line from the ligand library

    Output
    ------
        ligand_ID   (str): ligand identifier (or None if there is no ID)
        ligand_Name (str): ligand name (or None if there is no name)
    '''

    # Divide line by white space
    ligand_line = ligand_line.split()

    # If ligand_line was not empty or just whitespace
    if len(ligand_line) > 0:

        # Save ligand ID
        ligand_ID = ligand_line[0]

        # If there is a name
        if len(ligand_line)>1:
            
            # Save ligand name
            ligand_Name = ligand_line[1]
        
        else:
            # There is no name
            ligand_Name = None
    else:
        # There is no ID or name
        ligand_ID = None
        ligand_Name = None


    return ligand_ID, ligand_Name

def identifyFormat(ligand_ID):
    '''
    Guess the format of the identifier
    WARNING: there is currently no way of distinguishing 3 char SMILES from 3-letter PDB Ligand ID.
    
    Inputs
    ------
        ligand_ID (str): ligand identifier (SMILES, PDB ID or Drug Bank ID)

    Output
    ------
        ID_format (str): ligand_ID format ("PDB", "DB" or "SMILES")
    '''

 # Drug Bank format starts with 'DB'
    if (ligand_ID[0:2] == 'DB'):

        ID_format = 'DB'
    
 # PBD identifier
    elif (len(ligand_ID) == 3):

        ID_format = 'PDB'
    
 # Otherwise it should be a SMILES code 
    else:

        ID_format = 'SMILES'

    return ID_format

def sourceLigand(ligand_ID, ligand_Name, ligand_index, paths, prop):
    '''
    Sources ligand in sdf format from ligand identifier. 
    The accepted formats for the identifier are: SMILES, PDB ID and Drug Bank ID
    
    Inputs
    ------
        ligand_ID       (str): ligand identifier (SMILES, PDB ID or Drug Bank ID)
        ligand_Name     (str): ligand name
        ligand_index    (int): ligand index / counter
        paths          (dict): paths of step
        prop           (dict): properties of step
        
 
    Output
    ------
        output_path      (str): name of sdf file where ligand was saved
        successful_step (bool): success of step, False if output file was not created or is empty
                                Protection against corrupt ligand IDs
        ligand sdf file 
    '''

    # Identify format of ligand ID
    ID_format = identifyFormat(ligand_ID)

    # Modify the default output path 
    # According to ligand ID (unless SMILES is ID) or ligand Name if there is any
    addLigandSuffixToPaths(paths, ligand_ID, ligand_Name, ligand_index,
                            'output_sdf_path', 'output_path')

 # If ligand_ID is a PDB entry ID
    if (ID_format == 'PDB'):

        # Update properties 
        prop.update({'ligand_code': ligand_ID})                     

        # Action: retrieve ligand
        try:
            ideal_sdf(**paths, properties=prop)
        except:
            return False

 # If ligand_ID is a Drug Bank ID
    elif (ID_format == 'DB'):

        # Update properties
        prop.update({'drugbank_id': ligand_ID})

        # Action: retrieve ligand 
        try:
            drugbank(**paths, properties=prop)
        except:
            return False
        
 # If ligand_ID is a SMILES code
    elif (ID_format == 'SMILES'):
        
        # write the SMILES code in a file to give it as input to Open Babel 
        smiles_path = writeSMILES(SMILES = ligand_ID, 
                                 ligand_index = ligand_index, 
                                 step_path = prop['path'])

        # Update properties
        prop.update({'input_format': 'smiles'})
        prop.update({'output_format': 'sdf'})
        prop.update({'coordinates': 3})

        # Update paths
        paths.update({'input_path': smiles_path})

        # Action: format conversion using Open Babel (SMILES -> sdf)
        try:
            babel_convert(**paths, properties=prop)
        except:
            return False

        # Erase tmp SMILES file
        removeFiles(smiles_path)

    # Check output exists and is not empty (to skip corrupt ligand identifiers)
    successful_step = validateStep(paths['output_sdf_path']) or validateStep(paths['output_path'])

    return successful_step

def removeFiles(*file_path):
    '''
    Removes files in '*file_path' if they exist

    Inputs
    ------
        file_path  (str): variable number of paths to files including filename
    '''
    for path in file_path:
        if os.path.exists(path):
            os.remove(path)
    
    return

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

def addLigandSuffixToPaths(all_paths, ligand_ID, ligand_Name, ligand_index, *keywords):
    '''
    Modifies all_paths according to ligand_index, ligand_ID and ligand_Name if any. 
    
    Inputs
    ------

        all_paths      (dict):  dictionary with all paths
        ligand_ID       (str):  drug bank ID 
        ligand_Name     (str):  name of ligand
        ligand_index    (int):  ligand counter
        keywords    (str): keywords of output/input paths that will be modified
    
    Output
    ------

        suffix      (str): suffix used to modify paths
        (Implicit) all_paths paths are modified

    For example:

    path in all_paths = '/home/user/DockingWF/step4_drugbank/drug.sdf'
    ligand_ID = DB00530
    ligand_Name = None
    ligand_index = 4

    new_path = '/home/user/DockingWF/step4_drugbank/drug_4_DB00530.sdf'

    If 'ligand_ID' is a SMILES code then adds '_i_SMILES' to original_path before file extension
    '''
    
    # Identify ID format
    ID_format = identifyFormat(ligand_ID)

    # Add name information if any
    if ligand_Name is None:
        ligandName = ""
    else:
        ligandName = "_" + ligandName
    
    # SMILES code has problematic characters and is too long - just put "SMILES" in name
    if (ID_format == 'SMILES'):

        suffix = str(ligand_index) + ligandName + "_SMILES"
        
    # If format is PDB or DB - print code in file name
    else:
        suffix = str(ligand_index) + ligandName + "_" + str(ligand_ID) 
    
    # For all keys passed in keywords, modify path 
    for key in keywords:

        # Original path
        original_path = Path(all_paths[key])

        # Add suffix to path
        newpath = addSuffix(original_path, suffix)

        # Update paths dictionary
        all_paths.update({key : newpath})

    return suffix

def writeSMILES(SMILES, ligand_index, step_path):
    '''
    Writes a SMILES code into a file in step_path. The name of the file will
    be "ligand_name".smiles or ligand_"ligand_index".smiles

    Inputs
    ------
        SMILES              (str):  SMILES code
        ligand_index        (int):  ligand index (counter for the ligands)
        step_path    (Path class):  path of step where file will be written
    '''
    # Name of file
    smiles_filename = "ligand_" + str(ligand_index) + ".smiles"

    # Path including name
    smiles_path = os.path.join(str(step_path), smiles_filename) 

    # Save SMILES in tmp file inside step_path
    smiles_tmp_file = open(smiles_path, "w")
    smiles_tmp_file.write(SMILES)
    smiles_tmp_file.close()

    return smiles_path

def parallel_docking(ligands_queue, affinities_list, global_prop, global_paths):
    '''
    Function that will be used by different processes to dock all ligands in the input library

    Inputs
    ------
        ligands_queue    (Queue from multiprocessing): queue shared between processes with tuples containing 
                                                        the ligand index and corresponding line from input library
        affinities_list                        (list): results list containing the affinity for each ligand and ligand details 
    '''
    
    # Properties of steps
    prop_ligandLib  = global_prop["step4_source_lig"]
    prop_ligConv  = global_prop["step5_babel_prepare_lig"]
    prop_autodock  = global_prop["step6_autodock_vina_run"]
    prop_poseConv  = global_prop["step7_babel_prepare_pose"]

    # NOTE: if library is very big avoid prints per ligand, only warnings and errors - decide what to print and use lock to access same global_log

    while True:

        # Get item from shared queue
        queue_item = ligands_queue.get()

        # Each valid item has a ligand index and ligand information
        ligand_index, ligand_line = queue_item

        # If ligand information is None, we have reached the end, exit
        if ligand_line == None:
            return

        # Read ligand information
        ligand_ID, ligand_Name = readLigandLine(ligand_line)

        # Skip empty lines
        if ligand_ID is not None:

        # STEP 4: Source ligand in sdf format

            # Copy original paths to modify them
            paths_ligandLib = global_paths["step4_source_lig"].copy()

            # Add ligand index to log file names
            prop_ligandLib.update({'prefix' : str(ligand_index)})

            step4_successful = sourceLigand(ligand_ID, ligand_Name, ligand_index, 
                                                paths_ligandLib, prop_ligandLib)

            # Skip corrupt ligand IDs 
            if (step4_successful):
            
            # STEP 5: Convert ligand from sdf format to pdbqt format

                # Copy original paths to modify them
                paths_ligConv = global_paths["step5_babel_prepare_lig"].copy()
                
                # Modify paths according to ligand info
                ligand_identifier = addLigandSuffixToPaths(paths_ligConv, ligand_ID, ligand_Name, ligand_index,
                                        'input_path','output_path')
                
                # Add ligand index to log file names
                prop_ligConv.update({'prefix' : str(ligand_index)})
                
                # Action: format conversion using Open Babel
                babel_convert(**paths_ligConv, properties=prop_ligConv)


            # STEP 6: Autodock vina

                # Copy original paths to modify them
                paths_autodock = global_paths["step6_autodock_vina_run"].copy()            

                # Modify paths according to ligand info
                addLigandSuffixToPaths(paths_autodock, ligand_ID, ligand_Name, ligand_index,
                                        'output_pdbqt_path','output_log_path', 'input_ligand_pdbqt_path')

                # Add ligand index to log file names
                prop_autodock.update({'prefix' : str(ligand_index)})

                # Action: Run autodock vina
                autodock_vina_run(**paths_autodock, properties=prop_autodock)
                
                # Find best affinity among different poses - the rest can be checked from step7 output 
                affinity = findAffinityInLog(log_path = paths_autodock['output_log_path'])

                # Append result
                affinities_list.append((affinity, ligand_identifier))

            # STEP 7: Convert poses to PDB 

                # Copy original paths to modify them
                paths_poseConv = global_paths["step7_babel_prepare_pose"].copy()    

                # Modify paths according to ligand info
                addLigandSuffixToPaths(paths_poseConv, ligand_ID, ligand_Name, ligand_index,
                                        'input_path','output_path')

                # Add ligand index to log file names
                prop_poseConv.update({'prefix' : str(ligand_index)})

                # Action: Convert pose to PDB
                babel_convert(**paths_poseConv, properties=prop_poseConv)

                # Remove unnecessary files 
                removeFiles(paths_ligandLib['output_sdf_path'], 
                            paths_ligandLib['output_path'],
                            paths_ligConv['output_path'],
                            paths_autodock['output_pdbqt_path'])



def main_wf(configuration_path, ligand_lib_path, last_step = None, input_pockets_path = None, 
              pocket_ID = None, pocket_residues_path = None, input_structure_path = None):

    # Added last two arguments for ensemble docking :) NOTE: improve docs

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
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Set number of CPU cores
    
    # If we are in an HPC environment:
    if "SLURM_JOB_ID" in os.environ:
        num_cores = int(os.getenv('SLURM_CPUS_PER_TASK'))
    
    # If not, check num of CPUs
    else:
        num_cores = mp.cpu_count()
    
    # Initialize Manager with a list and a Queue
    manager = mp.Manager()
    affinities_list = manager.list()
    ligands_queue = manager.Queue(num_cores)

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

    # If "pocket_residues_path" provided, skip step 1, box will be formed using pocket residues instead of pocket .pqr file
    if pocket_residues_path is None:
    
    # STEP 1: Pocket selection from filtered list 

        # Write next action to global log
        global_log.info("step1_fpocket_select: Extract pocket cavity")

        # Properties and paths of step
        props = global_prop["step1_fpocket_select"]
        paths = global_paths["step1_fpocket_select"]

        # If model, pockets and pocket ID are provided through arguments -> prioritize over input.yml (ensemble docking)
        if None not in (input_pockets_path, pocket_ID, input_structure_path):

            paths.update({'input_pockets_zip' : input_pockets_path})
            props.update({'pocket' : pocket_ID})

        # Action: pocket selection
        fpocket_select(**paths, properties=props)

        # Validate pocket selection (make sure chosen pocket exists)
        if not validateStep(paths["output_pocket_pdb"],paths["output_pocket_pqr"]):

            # If output is not valid, print error message and available options in log file
            global_log.info("    ERROR: fpocket_select failed to select a pocket from input file.")
            global_log.info("           Check the selected pocket exists:")

            # Print available pockets in input folder
            printAvailablePockets(paths["input_pockets_zip"], global_log)

            return 

# STEP 2: Generate box around selected cavity

    # Write next action to global log
    global_log.info("step2_box: Generating cavity box")

    # Properties and paths of step
    props_box = global_prop["step2_box"]
    paths_box = global_paths["step2_box"]
    
    # If model and pocket_residues_path are provided through arguments -> prioritize over input.yml (ensemble docking with residues defining pocket)
    if None not in (pocket_residues_path, input_structure_path):

        paths_box.update({'input_pdb_path' : pocket_residues_path})

    # Action: box creation
    box(**paths_box, properties=props_box)

# STEP 3: Prepare target protein for docking 
    
    # Write next action to global log     
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    
    # Properties and paths of step
    props_addH = global_prop["step3_str_check_add_hydrogens"]
    paths_addH = global_paths["step3_str_check_add_hydrogens"]

    # If model, pockets and pocket ID are provided through arguments -> prioritize over input.yml (ensemble docking)
    if None not in (input_pockets_path, pocket_ID, input_structure_path):

        paths_addH.update({'input_structure_path' : input_structure_path})

    # Action: Preparation of target protein
    str_check_add_hydrogens(**paths_addH, properties=props_addH) 

    # Check if this should be the final step
    if last_step == 'prepare':
        global_log.info("Receptor preparation completed.")
        return 

# STEPS 4-5-6-7: For each ligand in library: obtain molecule, prepare ligand, run docking, prepare poses

    # Create directory for step 4. 
    # If SMILES are given, the code will try to save a txt file in this path with the SMILES before step4 is executed
    if not os.path.exists(global_prop['step4_source_lig']['path']):
        os.makedirs(global_prop['step4_source_lig']['path'])

    # Empty pool of workers
    pool = []

    # Start as many workers as available cores
    for i in range(num_cores):
        process = mp.Process(target = parallel_docking, args = (ligands_queue, affinities_list, global_prop.copy(), global_paths.copy()))
        process.start()
        pool.append(process)
    
    # Open ligand library, ligand_lib is an iterable obj
    with open(ligand_lib_path) as ligand_lib:

        # Add as many None as processes, they will be used as an exit condition
        extended_ligand_lib = itertools.chain(ligand_lib, (None,)*num_cores)

        # Fill queue with lines from ligand library
        for ligand_line in enumerate(extended_ligand_lib):

            ligands_queue.put(ligand_line)

    # Join all workers
    for process in pool:

        process.join()

# STEP 8: Find top ligands (according to lowest affinity)
    
    # Write action to global log
    global_log.info("step8_show_top_ligands: print identifiers of top ligands ranked by lowest affinity")  

    bestAffinities =  printTopLigands(affinities_list, global_prop['step8_show_top_ligands'], global_log)

    # Print timing info only if we are not calling docking_htvs from another script
    if None in (input_pockets_path, pocket_ID, input_structure_path):
        
        # Timing information
        elapsed_time = time.time() - start_time
        global_log.info('')
        global_log.info('')
        global_log.info('Execution successful: ')
        global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
        global_log.info('  Config File: %s' % configuration_path)
        global_log.info('  Ligand library: %s' % ligand_lib_path)
        global_log.info('')
        global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
        global_log.info('')

    return conf.get_working_dir_path(), bestAffinities

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('--lig-lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (Ligand PDB code, SMILES or Drug Bank ID) per line.",
                        required=True)
    
    # Execute workflow until 'to_do' step -> all executes all steps (all is default)
    parser.add_argument('--until', dest='to_do', 
                        help="(Opt) Extent of the pipeline to execute (preparation, all)", 
                        required=False)

    args = parser.parse_args()

    _,_= main_wf(args.config_path, args.ligand_lib, args.to_do)


