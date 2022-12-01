#!/usr/bin/env python3

# Importing all the needed libraries
import os
import re
import time
import argparse
from pathlib import Path
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

def findTopLigands(paths, properties, ligand_IDs, ligand_Names):
    '''
    Find top scoring ligands according to binding affinity calculated in step 6. 
    The number of ligands considered is determined in the properties of step 8 or 
    the maximum number of ligands
    
    Inputs
    ------
        paths         (dict): paths of step 8 
        properties    (dict): properties of step 8
        ligand_IDs    (list): list of ligand identifiers
        ligand_Names  (list): list of ligand names
        global_log  (Logger):

    Outputs
    -------
        bestLigandAffinities  (list): list with lowest affinities 
        bestLigandIDs         (list): list with corresponding ligand identifiers/names
    '''

# 1. Create 2 lists: one with binding affinities and another with ligandID_indices 

    affinities = []
    ligand_indices = []

    # Default path for docking log, convert from string to Path class
    dockingLogsPath = Path(paths["docking_logs"])
    
    # Step path for docking step
    stepPath = dockingLogsPath.parent

    # "output_vina*log"
    logNamePattern = dockingLogsPath.stem + "*" + dockingLogsPath.suffix

    # List all files matching the pattern
    logList = stepPath.rglob(logNamePattern)

    # For each docked ligand
    for log in logList:

        # Find index of ligand from log name
        logName = log.name
        ligand_index = findNumberInString(logName.strip(dockingLogsPath.stem))

        # Find affinity of best pose from log
        affinity = findMatchingStr(pattern=r'\s+1\s+(\S+)\s+', filepath=str(log))
        
        # If autodock was successful and we find an affinity in log
        if affinity is not None:
        
            # Save both: best affinity and index
            affinities.append(float(affinity))
            ligand_indices.append(ligand_index)

# 2. Find min and its index, save in another list, remove min. Repeat "number_top_ligands" times
    
    # Read how many ligands to include in the top
    numTopLigands = min(properties["number_top_ligands"], len(affinities)) 

    bestLigandAffinities = []
    bestLigandIDs = []
    bestLigandNames = []

    # Number of top ligands should not be larger than docked ligands
    for i in range(numTopLigands):

        # Find minimum affinity
        bestLigandAffinity = min(affinities)

        # Index of minimum affinity
        bestAffinityIndex = affinities.index(bestLigandAffinity)

        # Corresponding ligand ID index - as given in log file name
        bestLigandIndex = ligand_indices[bestAffinityIndex]

        # Save affinity
        bestLigandAffinities.append(bestLigandAffinity)

        # Save corresponding ligandID
        bestLigandIDs.append(ligand_IDs[bestLigandIndex])

        # Save corresponding ligand Name
        bestLigandNames.append(ligand_Names[bestLigandIndex])

        # Remove affinity
        affinities.remove(bestLigandAffinity)

        # Remove ligand index
        ligand_indices.remove(bestLigandIndex)

    return bestLigandAffinities, bestLigandIDs, bestLigandNames

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

def readLigandLibFile(ligand_lib_path):
    '''
    Read all ligand identifiers from ligand library file. 
    The expected format is one of the following:

    Format 1:

    ligand1_id
    ligand2_id
    .
    .
    .

    Format 2:

    ligand1_id  name_ligand1
    ligand2_id  name_ligand2
    .
    .
    .

    Where ligand_id is either a PDB code, Drug Bank ID or SMILES and name_ligand is a string with the ligand name
    
    Inputs
    ------
        ligand_lib_path (str): path pointing to ligand library, including file name

    Output
    ------
        ligand_IDs   (list(str)): list of ligand identifiers
        ligand_Names (list(str)): list with ligand names if any 
    '''

    # NOTE: How big will the library be? maybe we should process it by pieces...

    ligand_IDs = []
    ligand_Names = []

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Process every line
        for line in ligand_lib:

            # Divide line by white space
            line = line.split()

            # Append ligand ID to ligand_IDs list
            ligand_IDs.append(line[0])

            # If there is a name, append it to ligand_Names
            if len(line)>1:

                ligand_Names.append(line[1])
            
            else:
                
                ligand_Names.append(None)

    return ligand_IDs, ligand_Names

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
    successful_step = validateStep(paths['output_path']) or validateStep(paths['output_sdf_path'])

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

        all_paths  (dict): all new paths

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
        ligand_Name = ""
    else:
        ligand_Name = "_" + ligand_Name
    
    # SMILES code has problematic characters and is too long - just put "SMILES" in name
    if (ID_format == 'SMILES'):

        suffix = str(ligand_index) + ligand_Name + "_SMILES"
        
    # If format is PDB or DB - print code in file name
    else:
        suffix = str(ligand_index) + ligand_Name + "_" + str(ligand_ID) 
    
    # For all keys passed in keywords, modify path 
    for key in keywords:

        # Original path
        original_path = Path(all_paths[key])

        # Add suffix to path
        newpath = addSuffix(original_path, suffix)

        # Update paths dictionary
        all_paths.update({key : newpath})

    return all_paths

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

    # If step directory has not been created yet -> create dir 
    if not os.path.exists(step_path):
        os.makedirs(step_path)

    # Save SMILES in tmp file inside step_path
    smiles_tmp_file = open(smiles_path, "w")
    smiles_tmp_file.write(SMILES)
    smiles_tmp_file.close()

    return smiles_path

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

    # Launching the actions of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

    # If "pocket_residues_path" provided, skip step 1, box will be formed using pocket residues
    if pocket_residues_path is None:
    # STEP 1: Pocket selection from filtered list 

        # Write next action to global log
        global_log.info("step1_fpocket_select: Extract pocket cavity")

        # Properties and paths of step
        props = global_prop["step1_fpocket_select"]
        paths = global_paths["step1_fpocket_select"]

        # If model, pockets and pocket ID are provided through arguments -> prioritize over input.yml (ensemble docking)
        if None not in (input_pockets_path, pocket_ID):

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
    if pocket_residues_path is not None:

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
    if input_structure_path is not None:

        paths_addH.update({'input_structure_path' : input_structure_path})

    # Action: Preparation of target protein
    str_check_add_hydrogens(**paths_addH, properties=props_addH) 

    # Check if this should be the final step
    if last_step == 'prepare':
        global_log.info("Receptor preparation completed.")
        return 

# STEP 4-5-6-7: For each ligand in library: obtain molecule, prepare ligand, run docking, prepare poses
    
    # Properties of next steps
    prop_ligandLib  = global_prop["step4_source_lig"]
    prop_ligConv  = global_prop["step5_babel_prepare_lig"]
    prop_autodock  = global_prop["step6_autodock_vina_run"]
    prop_poseConv  = global_prop["step7_babel_prepare_pose"]

    # Load drug list
    if (ligand_lib_path):
        # If file with drug library is given -> prioritize over list in input.yml
        ligand_IDs, ligand_Names = readLigandLibFile(ligand_lib_path)
    else:
        # If no file is given, use list from input 
        ligand_IDs = prop_ligandLib['ligand_list']
        ligand_Names = [None]*len(ligand_IDs)

    # NOTE: if library is very big avoid prints per ligand, only warnings and errors - decide what to print and use lock to access same global_log

    for ligand_index in range(len(ligand_IDs)):

        ligand_ID = ligand_IDs[ligand_index]
        ligand_Name = ligand_Names[ligand_index]

    # STEP 4: Source ligand in sdf format

        # Write next action to global log
        global_log.info("step4_source_lig: Extracting small molecule from library")

        # Copy original paths to modify them
        paths_ligandLib = global_paths["step4_source_lig"].copy()

        step4_successful = sourceLigand(ligand_ID, ligand_Name, ligand_index, 
                                            paths_ligandLib, prop_ligandLib)

        # Skip corrupt ligand IDs 
        if (step4_successful):
        
        # STEP 5: Convert from sdf format to pdbqt format

            # Write next action to global log
            global_log.info("step5_babel_prepare_lig: Preparing small molecule (ligand) for docking")

            # Copy original paths to modify them
            paths_ligConv = global_paths["step5_babel_prepare_lig"].copy()
            
            # Modify paths according to ligand info
            addLigandSuffixToPaths(paths_ligConv, ligand_ID, ligand_Name, ligand_index,
                                    'input_path','output_path')

            # Action: format conversion using Open Babel
            babel_convert(**paths_ligConv, properties=prop_ligConv)


        # STEP 6: Autodock vina

            # Write action to global log
            global_log.info("step6_autodock_vina_run: Running the docking")

            # Copy original paths to modify them
            paths_autodock = global_paths["step6_autodock_vina_run"].copy()            

            # Modify paths according to ligand info
            addLigandSuffixToPaths(paths_autodock, ligand_ID, ligand_Name, ligand_index,
                                    'output_pdbqt_path','output_log_path', 'input_ligand_pdbqt_path')

            # Action: Run autodock vina
            autodock_vina_run(**paths_autodock, properties=prop_autodock)


        # STEP 7: Convert poses to PDB

            # Write action to global log
            global_log.info("step7_babel_prepare_pose: Converting ligand pose to PDB format")  

            # Copy original paths to modify them
            paths_poseConv = global_paths["step7_babel_prepare_pose"].copy()    

            # Modify paths according to ligand info
            addLigandSuffixToPaths(paths_poseConv, ligand_ID, ligand_Name, ligand_index,
                                    'input_path','output_path')

            # Action: Convert pose to PDB
            babel_convert(**paths_poseConv, properties=prop_poseConv)

            # Remove unnecessary files 
            removeFiles(paths_ligandLib['output_sdf_path'], 
                        paths_ligandLib['output_path'],
                        paths_ligConv['output_path'],
                        paths_autodock['output_pdbqt_path'])
        
        # If step 4 was not successful
        else:

            global_log.info("    WARNING: failed to source ligand {} with ligand ID: {}".format(ligand_index, ligand_ID))
            global_log.info("            Skipping to next ligand")

# STEP 8: Find top ligands (according to lowest affinity)
    
    # Write action to global log
    global_log.info("step8_show_top_ligands: print identifiers of top ligands ranked by lowest affinity")  

    # Find and print top ligands in log file
    bestAffinities, bestLigandIDs, bestLigandNames = findTopLigands(paths = global_paths["step8_show_top_ligands"], 
                                                            properties = global_prop["step8_show_top_ligands"], 
                                                            ligand_IDs = ligand_IDs, ligand_Names = ligand_Names)

    # Print more info only if we are not calling docking_htvs from another script
    if (input_pockets_path == None) and (pocket_ID == None) and (input_structure_path == None):

        # Print best ligands and their affinities in log file
        for i in range(len(bestAffinities)):

            global_log.info("    Affinity: {} for ligand {}".format(bestAffinities[i], bestLigandIDs[i]))
        
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

    return conf.get_working_dir_path(), bestAffinities, bestLigandIDs, bestLigandNames

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")
    
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('--lig-lib', dest='ligand_lib',
                        help="Path to file with ligand library. The file should contain one ligand identifier (Ligand PDB code, SMILES or Drug Bank ID) per line.",
                        required=False)
    
    # Execute workflow until 'to_do' step -> all executes all steps (all is default)
    parser.add_argument('--until', dest='to_do', 
                        help="(Opt) Extent of the pipeline to execute (preparation, all)", 
                        required=False)

    args = parser.parse_args()

    _,_,_,_= main_wf(args.config_path, args.ligand_lib, args.to_do)
