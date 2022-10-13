#!/usr/bin/env python3

# Importing all the needed libraries
from mailbox import NoSuchMailboxError
import os
import re
import time
import argparse
from pathlib import Path, PurePath
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

    Inputs:
        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output:
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

def findMatchingLine(pattern, filepath):
    '''
    Finds first line in file containing a given pattern

    Inputs:
        pattern  (regex pattern): regular expression pattern to search in file lines
        filepath           (str): file path to search in
    
    Output:
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

            return line
    
    return None

def findNumberInString(string):
    '''
    Finds and returns first integer number in string. If no integer is found a 0 is returned.

    Inputs:
        string (str): string to search int in
    
    Output:
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
    
    Inputs:
        pockets_path  (str): Path to pockets folder including file name
        global_log (Logger): Object of class Logger (dumps info to global log file of this workflow)
    Output:
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

    # Iterate through all log files -> print the model ID + available pockets in the global log
    for log in logList:

        # Find model from file name
        logName = log.name
        modelID = findNumberInString(logName.strip(stepPath.name))

        # Find pockets 
        pockets = findMatchingLines(pattern=r'pocket\d+$', filepath=str(log))

        # Print to global log
        global_log.info("   Model {}".format(modelID))

        for pocket in pockets:
            global_log.info("        {}".format(pocket))

    return

def findTopLigands(paths, properties, ligands, global_log):
    '''
    Find top scoring ligands according to binding affinity calculated in step 6. 
    The number of ligands considered is determined in the properties of step 8
    
    Inputs:
        paths         (dict): paths of step 8 
        properties    (dict): properties of step 8
        ligands       (list): list of ligand identifiers
        global_log  (Logger):

    Output:
        bestLigandAffinities      (list): list with lowest affinities 
        bestLigandIDs             (list): list with corresponding ligand identifiers
    '''

# 1. Create 2 lists: one with binding affinities and another with ligandID_indices 

    affinities = []
    ligandID_indices = []

    # Default path for docking log, convert from string to Path class
    dockingLogsPath = Path(paths["docking_logs"])
    
    # Step path for docking step
    stepPath = dockingLogsPath.parent

    # "output_vina*log"
    logNamePattern = dockingLogsPath.stem + "*" + dockingLogsPath.suffix

    # List all files matching the pattern
    logList = stepPath.rglob(logNamePattern)

    # NOTE: extend to many ligands - might be 2 or 3 characters, take number 

    # For each docked ligand
    for log in logList:

        # Find index of ligand from log name
        logName = log.name
        ligandID_index = findNumberInString(logName.strip(dockingLogsPath.stem))

        # Find affinity of best pose from log
        line = findMatchingLine(pattern=r'   1 ', filepath=str(log))
        
        # Separate line
        lineParts = line.split()

        # Save just the affinity
        affinity = float(lineParts[1])

        # Save both: best affinity and index
        affinities.append(affinity)
        ligandID_indices.append(ligandID_index)

# 2. Find min and its index, save in another list, remove min. Repeat "number_top_ligands" times
    
    # Read how many ligands to include in the top
    numTopLigands = properties["number_top_ligands"]

    bestLigandAffinities = []
    bestLigandIDs = []

    # Number of top ligands should not be larger than docked ligands
    for i in range(min(numTopLigands, len(affinities))):

        # Find minimum affinity
        bestLigandAffinity = min(affinities)

        # Index of minimum affinity
        bestAffinityIndex = affinities.index(bestLigandAffinity)

        # Corresponding ligand ID index - as given in log file name
        bestLigandIDIndex = ligandID_indices[bestAffinityIndex]

        # Save affinity
        bestLigandAffinities.append(bestLigandAffinity)

        # Save corresponding ligandID
        bestLigandIDs.append(ligands[bestLigandIDIndex])

        # Remove affinity
        affinities.remove(bestLigandAffinity)

        # Remove ligandID index
        ligandID_indices.remove(bestLigandIDIndex)

    # Print best ligands and their affinities in log file
    for i in range(len(bestLigandAffinities)):

        global_log.info("    Affinity: {} for ligand {}".format(bestLigandAffinities[i], bestLigandIDs[i]))

    return

def validateStep(*output_paths):
    '''
    Check all output files exist and are not empty
    
    Inputs:
        *output_paths (str): variable number of paths to output file/s

    Output:
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

def readLigandLibFile(ligand_lib_path):
    '''
    Read all ligand identifiers from ligand library file
    
    Inputs:
        ligand_lib_path (str): path pointing to ligand library, including file name

    Output:
        ligand_lib (list(str)): list of ligand identifiers 
    '''

    # NOTE: How big will the library be? maybe we should process it by pieces...

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Strip whitespace
        ligand_lib = [line.rstrip() for line in ligand_lib]

    return ligand_lib


def identifyFormat(ligand_ID):
    '''
    Guess the format of the identifier
    WARNING: there is currently no way of distinguishing 3 char SMILES from 3-letter PDB Ligand ID.
    
    Inputs:
        ligand_ID (str): ligand identifier (SMILES, PDB ID or Drug Bank ID)

    Output:
        ID_format (str): ligand_ID format ("PDB", "DB" or "SMILES")
    '''

 # Drug Bank format starts with 'DB'
    if (ligand_ID[0:2] == 'DB'):

        ID_format = 'DB'
    
 # PBD identifier
    # NOTE: should they be just letters? not numbers?
    elif (len(ligand_ID) == 3):

        ID_format = 'PDB'
    
 # Otherwise it should be a SMILES code 
    else:

        ID_format = 'SMILES'

    return ID_format

def sourceLigand(ligand_ID, default_output_path, paths, prop, ligand_index):
    '''
    Sources ligand in sdf format from ligand identifier. 
    The accepted formats for the identifier are: SMILES, PDB ID and Drug Bank ID
    
    Inputs:
        ligand_ID                  (str): ligand identifier (SMILES, PDB ID or Drug Bank ID)
        default_output_path (Path class): path chosen by default for output
        paths                     (dict): paths of step
        prop                      (dict): properties of step
        ligand_index               (int): index of ligand
 
    Output:
        output_path      (str): name of sdf file where ligand was saved
        successful_step (bool): success of step, False if output file was not created or is empty
                                Protection against corrupt ligand IDs
        ligand sdf file 
    '''
    # Identify format of ligand ID
    ID_format = identifyFormat(ligand_ID)

    # If format is SMILES, we will need to write the SMILES code in a file to give it as input to Open Babel (SMILES -> sdf)
    smiles_filename = "ligand_" + str(ligand_index) + ".smiles"
    step_directory = default_output_path.parent
    smiles_path = str(step_directory) + "/" + smiles_filename

 # If ligand_ID is a PDB entry ID
    if (ID_format == 'PDB'):

        # Update properties 
        prop.update({'ligand_code': ligand_ID})
        prop.update({'api_id': 'pdbe'})                 # NOTE: pdbe or pdb?

        # Modify the output path - according to drug ID
        output_path = addLigandIDSuffix(default_output_path, ligand_ID, ligand_index)

        # Update paths dictionary
        paths.update({'output_sdf_path': output_path})                             

        # Action: retrieve ligand
        try:
            ideal_sdf(**paths, properties=prop)
        except:
            return output_path, False

 # If ligand_ID is a Drug Bank ID
    elif (ID_format == 'DB'):

        # Update properties
        prop.update({'drugbank_id': ligand_ID})

        # Modify the output path - according to drug ID
        output_path = addLigandIDSuffix(default_output_path, ligand_ID, ligand_index)

        # Update paths dictionary
        paths.update({'output_sdf_path': output_path})

        # Action: retrieve ligand 
        try:
            drugbank(**paths, properties=prop)
        except:
            return output_path, False
        

 # If ligand_ID is a SMILES code
    elif (ID_format == 'SMILES'):

        # Update properties
        prop.update({'input_format': 'smiles'})
        prop.update({'output_format': 'sdf'})
        prop.update({'coordinates': 3})

        # If step4 directory has not been created yet -> create dir 
        if not os.path.exists(step_directory):
            os.makedirs(step_directory)

        # Save SMILES in tmp file inside step_directory
        smiles_tmp_file = open(smiles_path, "w")
        smiles_tmp_file.write(ligand_ID)
        smiles_tmp_file.close()

        # Modify the output path - according to drug ID
        output_path = addLigandIDSuffix(default_output_path, ligand_ID, ligand_index)

        # Update paths
        paths.update({'input_path': smiles_path})
        paths.update({'output_path': output_path})

        # Action: format conversion using Open Babel (SMILES -> sdf)
        try:
            babel_convert(**paths, properties=prop)
        except:
            return output_path, False

        # Erase tmp file
        removeFiles(smiles_path)
    
    # Check output exists and is not empty (to skip corrupt ligand identifiers)
    successful_step = validateStep(output_path)

    return output_path, successful_step

def removeFiles(*file_path):
    '''
    Removes files in '*file_path' if they exist

    Inputs
        file_path  (str): variable number of paths to files including filename
    '''
    for path in file_path:
        if os.path.exists(path):
            os.remove(path)
    
    return

def addLigandIDSuffix(original_path, ligand_ID, ligand_index):
    '''
    Adds 'ligand_ID' to original_path before file extension. For example:

    original_path = Path('/home/user/DockingWF/step4_drugbank/drug.sdf')
    ligand_ID = DB00530

    new_path = '/home/user/DockingWF/step4_drugbank/drug_DB00530.sdf'

    If 'ligand_ID' is a SMILES code then adds '_i_SMILES' to original_path before file extension

    Inputs
        original_path  (Path class):  path to original file including filename
        ligand_ID             (str):  drug bank ID 
    '''
    
    ID_format = identifyFormat(ligand_ID)

    filename = original_path.name
    directory = original_path.parent

    filename_parts = filename.split(".")

    # If format is PDB or DB - print code in file name
    if (ID_format != 'SMILES'):
        new_path = str(directory) + "/" + filename_parts[0] + "_" + str(ligand_index) + "_" + str(ligand_ID) + "." + filename_parts[1]
    
    # SMILES code has problematic characters and is too long - just put SMILES instead
    else:
        new_path = str(directory) + "/" + filename_parts[0] + "_" + str(ligand_index) + "_SMILES" + "." + filename_parts[1]

    return new_path


def main(args):

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.config_path)
    
    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    
    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 


# STEP 1: Pocket selection from filtered list 

    # Write next action to global log
    global_log.info("step1_fpocket_select: Extract pocket cavity")

    # Properties and paths of step
    props = global_prop["step1_fpocket_select"]
    paths = global_paths["step1_fpocket_select"]

    # Action: pocket selection
    fpocket_select(**paths, properties=props)

    # Validate pocket selection (make sure chosen pocket exists)
    if not validateStep(paths["output_pocket_pdb"],paths["output_pocket_pqr"]):

        # If output is not valid, print error message and available options in log file
        global_log.info("    ERROR: fpocket_select failed to select a pocket from input file.")
        global_log.info("           Check the selected pocket in the properties of this step exists.")

        # Print available pockets in input folder
        printAvailablePockets(paths["input_pockets_zip"], global_log)

        # Exit 
        return

# STEP 2: Generate box around selected cavity

    # NOTE: should we adjust offset depending on largest molecule in database?

    # Write next action to global log
    global_log.info("step2_box: Generating cavity box")
    
    # Action: box creation
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])


# STEP 3: Prepare target protein for docking 

    # NOTE: try to fix case in which residues start again - for more complicated proteins pre-processed with maestro or else
    
    # Write next action to global log     
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    
    # Action: Preparation of target protein
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 


# STEP 4-5-6-7: For each ligand in library: obtain molecule, prepare ligand, run docking, prepare poses

    # Write next action to global log
    global_log.info("step4_source_lig: Extracting small molecules from library")
    
    # Properties and paths of next steps
    prop_ligandLib  = global_prop["step4_source_lig"]
    paths_ligandLib = global_paths["step4_source_lig"]

    prop_ligConv  = global_prop["step5_babel_prepare_lig"]
    paths_ligConv = global_paths["step5_babel_prepare_lig"]

    prop_autodock  = global_prop["step6_autodock_vina_run"]
    paths_autodock = global_paths["step6_autodock_vina_run"]

    prop_poseConv  = global_prop["step7_babel_prepare_pose"]
    paths_poseConv = global_paths["step7_babel_prepare_pose"]

    # Load drug list
    if (args.ligand_lib):
        # If file with drug library is given
        ligand_list = readLigandLibFile(args.ligand_lib)
    else:
        # If no file is given, use list from input
        ligand_list = prop_ligandLib['ligand_list']

    # Total number of drug identifiers
    num_ligandIDs = len(ligand_list)

    # If num_ligandIDs > 20 then save minimum number of files
    drugLibIsLarge = num_ligandIDs > 20

    # Default paths that will be modified according to drug ID
    default_ligandLib_output = Path(paths_ligandLib['output_sdf_path'])
    default_ligConv_output = Path(paths_ligConv['output_path'])
    default_autodock_pdbqt = Path(paths_autodock['output_pdbqt_path'])
    default_autodock_log = Path(paths_autodock['output_log_path'])

    default_poseConv_input = Path(paths_poseConv['input_path'])
    default_poseConv_output = Path(paths_poseConv['output_path'])

    ligandID_index = 0

    for ligand in ligand_list:

    # STEP 4: Source ligand in sdf format

        drugbank_output_path, step4_successful = sourceLigand(ligand, default_ligandLib_output, 
                                                               paths_ligandLib, prop_ligandLib, ligandID_index)

        # Skip corrupt ligand IDs 
        if (step4_successful):
        
        # STEP 5: Convert from sdf format to pdbqt format

            # Write next action to global log
            global_log.info("step5_babel_prepare_lig: Preparing small molecule (ligand) for docking")

            # Modify the output path - according to drug ID
            ligConv_output_path = addLigandIDSuffix(default_ligConv_output, ligand, ligandID_index)
        
            # Update paths dictionary
            paths_ligConv.update({'input_path':drugbank_output_path})
            paths_ligConv.update({'output_path':ligConv_output_path})

            # Action: format conversion using Open Babel
            babel_convert(**paths_ligConv, properties=prop_ligConv)


        # STEP 6: Autodock vina

            # Write action to global log
            global_log.info("step6_autodock_vina_run: Running the docking")

            # Modify pdbqt and log paths - according to drug ID
            autodock_pdbqt_path = addLigandIDSuffix(default_autodock_pdbqt, ligand, ligandID_index)
            autodock_log_path = addLigandIDSuffix(default_autodock_log, ligand, ligandID_index)

            # Update paths in dictionary  
            paths_autodock.update({'input_ligand_pdbqt_path': ligConv_output_path})
            paths_autodock.update({'output_pdbqt_path': autodock_pdbqt_path})
            paths_autodock.update({'output_log_path': autodock_log_path})

            # Action: Run autodock vina
            autodock_vina_run(**paths_autodock, properties=prop_autodock)


        # STEP 7: Convert poses to PDB

            # Write action to global log
            global_log.info("step7_babel_prepare_pose: Converting ligand pose to PDB format")  

            # Modify input and output paths - according to drug ID
            poseConv_input_path = addLigandIDSuffix(default_poseConv_input, ligand, ligandID_index)
            poseConv_output_path = addLigandIDSuffix(default_poseConv_output, ligand, ligandID_index)

            # Update paths in dictionary
            paths_poseConv.update({'input_path': poseConv_input_path})
            paths_poseConv.update({'output_path': poseConv_output_path})

            # Action: Convert pose to PDB
            babel_convert(**paths_poseConv, properties=prop_poseConv)

            # If drug library is large remove unnecessary files 
            if (drugLibIsLarge):
                removeFiles(drugbank_output_path, ligConv_output_path, autodock_pdbqt_path)
        
        else:

            global_log.info("    WARNING: failed to source ligand {} with ligand ID: {}".format(ligandID_index, ligand))
            global_log.info("            Skipping to next ligand")
        
        # Increase ligand ID count
        ligandID_index += 1

# STEP 8: Find top ligands (according to lowest affinity)
    
    # Write action to global log
    global_log.info("step8_show_top_ligands: print identifiers of top ligands ranked by lowest affinity")  

    # Properties and paths of step
    props = global_prop["step8_show_top_ligands"]
    paths = global_paths["step8_show_top_ligands"]

    # Find and print top ligands in log file
    findTopLigands(paths, props, ligand_list, global_log)

# Timing information
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % args.config_path)
    global_log.info('  Ligand library: %s' % args.ligand_lib)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)")

    parser.add_argument('--lig-lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (Ligand PDB code, SMILES or Drug Bank ID) per line.")
    

    args = parser.parse_args()

    main(args)
