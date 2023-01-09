##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import sys
import yaml
import glob
import time
import logging
import argparse
import numpy as np
import mdtraj as md
from scipy import stats
from glob import glob
from pathlib import PurePath, Path
from zipfile import ZipFile

# Import biobb
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.grompp_mdrun import grompp_mdrun

def checkFilesExist(*file_path):
    '''
    Returns true if all files exist

    Inputs
    ------

        file_path  (str): variable number of paths to files including filename
    
    Output
    ------

        exist (bool): True if they all exist. False if any of them doesnt exist
    '''

    exist = True

    for path in file_path:

        exist = exist and os.path.isfile(path)

        if not os.path.isfile(path):

            logging.warning(" {} file was not found".format(path))
            
    return exist

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

def createFolder(dir_path):
    '''
    Creates folder dir_path if it doesn't exist

    Inputs
    ------

        dir_path  (str): path including folder name to create
    '''

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    
    return

def checkOutputExistence(output_path, file_name):
    """
    Check the existence of a file called file_name in output_path. 
    If necessary, creates a new name adding a number suffix -> new_file_name
    Returns the complete path to the new file name output_path/new_file_name

    Inputs
    ----------

        output_path    (str): path to output folder where the file will be written
        file_name      (str): name of file
    
    Returns
    -------
        
        file_path    (str): path to new file
    """

    # file name as Path
    file_path = PurePath(file_name)

    # Construct generic pattern for file, example: "rave_analysis*.log"
    file_name_pattern = file_path.stem + "*" + file_path.suffix

    # List all files matching the pattern in the output path
    files_list = list(glob(os.path.join(output_path, file_name_pattern)))

    # If no file matches the pattern
    if len(files_list) == 0:    

        # No need to modify the original name
        file_path = os.path.join(output_path, file_name)

    # If N files match the pattern
    else:

        # We modify the filename with N
        new_file_name = file_path.stem + str(len(files_list)) + file_path.suffix

        file_path = os.path.join(output_path, new_file_name)

    return file_path

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

def launchContinuationMD():
    """
    Launch MD using last accepted structure, last accepted checkpoint and without generating new velocities.

    Inputs
    ------

    (Implicitly)
    Paths as global variables
    MD simulation settings as global dictionary
    
    Output
    ------
    
    lastAcceptedStep = False

    (Implicitly)
    Ouput files of short MD simulation
        
    """

    logging.info('+++++++++ Continuation of previous step...')

    # Modify properties
    md_prop.update({'continuation' : 'yes'})
    md_prop.update({'gen-vel' : 'no'})
    _ = md_prop.pop('gen-temp', None)

    # Execute a short MD simulation, restarting from last checkpoint
    grompp_mdrun(input_gro_path=lastAcceptedGRO,
            input_top_zip_path=args.topology_path,
            input_cpt_path=lastAcceptedCPT,
            output_tpr_path=step_tpr_path,
            output_trr_path=step_trr_path,
            output_gro_path=step_gro_path,
            output_edr_path=step_edr_path,
            output_log_path=step_log_path,
            output_xtc_path=step_xtc_path,
            output_cpt_path=step_cpt_path,
            properties=md_prop)

    lastStepAccepted = False

    return lastStepAccepted

def launchNewVelocitiesMD():
    """
    Launch MD using last accepted structure and generating new velocities.

    Inputs
    ------

    (Implicitly)
    Paths as global variables
    MD simulation settings as global dictionary
    
    Output
    ------

    (Implicitly)
    Ouput files of short MD simulation
        
    """

    logging.info('+++++++++ Generating new velocities...')

    # Modify properties
    md_prop.update({'continuation' : 'no'})
    md_prop.update({'gen-vel' : 'yes'})
    md_prop.update({'gen-temp' : temperature})

    # Execute a short MD simulation, generating new velocities
    grompp_mdrun(input_gro_path=lastAcceptedGRO,
            input_top_zip_path=args.topology_path,
            output_tpr_path= step_tpr_path,
            output_trr_path=step_trr_path,
            output_gro_path=step_gro_path,
            output_edr_path=step_edr_path,
            output_log_path=step_log_path,
            output_xtc_path=step_xtc_path,
            output_cpt_path=step_cpt_path,
            properties=md_prop)
    return

def computeCVslope():
    """
    Analyze user-defined CV. The CV is defined as the distance between two atom selections from
    the YAML configuration. Fits a line to the evolution of the CV and returns the CV average value and the slope of the fitted line

    Inputs
    ------

    (Implicitly)
    Trajectory and topology
    Atomic selections in CONFIG dictionary
    
    Output
    ------

    CV_mean: average value of CV during short MD
    CV_slope: slope of the fitted line to the CV evolution 
    """

    logging.info('+++++++++ Analyzing CV...')

    # Read trajectory and topology
    step_traj = md.load(step_xtc_path, top=step_gro_path)

    # Find the COM of each group selection for each frame (nm)
    com_group1 = md.compute_center_of_mass(traj = step_traj, select = CONFIG['sumd']['colvar']['group1_selection'])
    com_group2 = md.compute_center_of_mass(traj = step_traj, select = CONFIG['sumd']['colvar']['group2_selection'])

    # Find the vector joining the COMs for each frame
    distance_vectors = com_group1 - com_group2

    # Find the norm of this vector for each frame (nm)
    CV = np.linalg.norm(distance_vectors, axis=1)

    logging.info('+++++++++ CV: {}'.format(CV))

    # Time array (index of frames)
    time_au = np.arange(0, len(CV))

    # Mean value of the CV
    CV_mean = np.mean(CV)

    logging.info('+++++++++ CV mean: {}'.format(CV_mean))

    # Compute linear regression 
    res = stats.linregress(time_au, CV)

    # Find slope in nm/frame
    CV_slope = res.slope

    # Find slope m in nm/ns
    CV_slope = CV_slope/traj_freq*1000

    logging.info('+++++++++ CV slope: {}'.format(CV_slope))

    return CV_mean, CV_slope

def acceptLastStep():
    """
    Accept last short MD simulation:

        - Add short trajectory to total trajectory
        - Set the step checkpoint file as the lastAcceptedCPT
        - Set the step gro file as the lastAcceptedGRO

    Inputs
    ------

    (Implicitly)
    Total trajectory path
    Step trajectory path
    Step gro path
    Step cpt path
    Last accepted gro path
    Last accepted cpt path

    Output
    ------

    (Implicitly)
    Files renamed and saved
    """

    logging.info('+++++++++ Accepting step!')

    # Save xtc in final trajectory
    if os.path.exists(trajPath):
        
        # remove previous zip trajectory bundle if it exists
        removeFiles(zipPath)

        # Create a ZipFile object
        zipObject = ZipFile(zipPath, 'w')
        
        # Add traj to zip file
        zipObject.write(trajPath)
        zipObject.write(step_xtc_path)

        # Close zip file
        zipObject.close()

        prop = { 'concatenate': True }

        # Concatenate
        trjcat(input_trj_zip_path=zipPath,
                output_trj_path=trajPath,
                properties=prop)

    else:

        # Change name
        os.rename(step_xtc_path, trajPath)

    # Save checkpoint
    os.rename(step_cpt_path, lastAcceptedCPT)

    # Save GRO file
    os.rename(step_gro_path, lastAcceptedGRO)

    return

def main_wf():
    '''
    Main workflow that takes as input the GROMACs topology and structure files together with a certain configuration as global variables to
    perform Supervised Molecular Dynamics. 

    The main idea is to run short MDs and track the evolution of a certain CV defined by the user. If the CV value approaches the 
    target value during the short MD simulation, then the final state of the system is accepted and used to perform the next short MD.
    However, if the CV value departs from the target value, then the next simulation re-starts from the last accepted state generating 
    new velocities to try to explore a different direction in conformational space.

    Inputs
    ------

    System structure (.gro)
    gromacs topology (.zip of .top)
    Configuration (different global variables)

    Outputs
    -------
    '''

    # Start timer
    start_time = time.time()

    # Initialize total step counter 
    step_counter = 0

    # While all conditions True, keep running steps
    while (step_counter < maxNumSteps) and notWithinThreshold:

        logging.info('++++++ STEP {}'.format(step_counter))

        if lastStepAccepted:

            # If last step was accepted, continue MD with same velocities 
            lastStepAccepted = launchContinuationMD()

        else:
            
            # Otherwise, restart from last accepted structure with new velocities 
            launchNewVelocitiesMD()

        # Use MDTraj to analyze a certain CV (distance between 2 user-defined groups) 
        CV_mean, CV_slope = computeCVslope()

        # Check if distance to target is smaller than CV threshold
        if abs(CV_mean - CV_target) < CV_threshold:

            # We are within threshold!
            notWithinThreshold = False

            logging.info('+++++++++ CV is within threshold of target!')
        
        # Check slope is higher than CV slope threshold
        if abs(CV_slope) > slopeThreshold:

            # Accept new structure if CV_mean ----> CV_target and m > 0
            if (CV_mean < CV_target) and (CV_slope > 0):
                
                # Accept last step
                acceptLastStep()
 
                lastStepAccepted = True

            # Accept new structure if CV_target <---- CV_mean and m < 0
            elif (CV_mean > CV_target) and (CV_slope < 0): 
                
                # Accept last step 
                acceptLastStep()

                lastStepAccepted = True

        # Increase total counter
        step_counter += 1
    
    # Print timing information to the log file
    elapsed_time = time.time() - start_time

    logging.info('+++++++++ Elapsed time: {} minutes'.format(round(elapsed_time/60, 1)))

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('-input', dest='input_path',
                        help="Path to input structure (.gro)", 
                        required=True)

    parser.add_argument('-topology', dest='topology_path',
                        help="Path to input topology (.zip)", 
                        required=True)

    parser.add_argument('-config', dest='configuration',
                        help="Path to configuration file (YAML)", 
                        required=True)

    parser.add_argument('-output', dest='output_path',
                        help="Path where results will be dumped (./output by default)", 
                        required=False)

    args = parser.parse_args()

    
    if args.output_path:
        # If output path is given
        OUTPUT_PATH = args.output_path
    else:
        # If not -> default output path
        OUTPUT_PATH = "output"

    # Create folder if necessary
    createFolder(OUTPUT_PATH)

    # Path to global log file
    log_path = checkOutputExistence(OUTPUT_PATH, "suMD.log")

    # Determine logging level, format and handlers to print to log file and console
    logging.basicConfig( 
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)   
        ])

    # Print log title
    logging.info("SuMD LOG FILE\n")

    # Parse configuration file   
    if checkFilesExist(args.configuration):
        # If configuration file exists -> open configuration file
        with open(args.configuration) as config_file:

            # Load configuration file as dictionary
            CONFIG = yaml.load(config_file, Loader = yaml.FullLoader)
    else:
        # Issue error
        logging.error("Configuration file not found!")
        logging.error("Check the configuration file exists: {}".format(args.configuration))
        sys.exit()

    # Read common MD simulation properties from input configuration
    md_prop = CONFIG['md']['properties'].copy()

    # Find some constants
    CV_target = CONFIG['sumd']['colvar']['target']
    CV_threshold = CONFIG['sumd']['colvar']['threshold']
    maxNumSteps = CONFIG['sumd']['num_steps']['max_total']
    temperature = CONFIG['md']['properties']['mdp']['gen-temp']
    slopeThreshold = CONFIG['sumd']['slope_threshold']
    trajName = CONFIG['sumd']['files']['trajectory']
    
    # Trajectory saving period in ps: time step (ps) * saving frequency (time steps)
    traj_freq = md_prop['mdp']['dt']*md_prop['mdp']['nstxout-compressed']

    # Path to trajectory file: if there are previous trajectories, add number to new trajectory to avoid overwriting
    trajPath = checkOutputExistence(OUTPUT_PATH, trajName)

    # Path to zip file 
    zipPath = checkOutputExistence(OUTPUT_PATH, 'trajectory_bundle.zip') 

    # Names of step temporal files (short MD results)
    step_tpr_path = os.path.join(OUTPUT_PATH, "step.tpr")
    step_trr_path = os.path.join(OUTPUT_PATH, "step.trr")
    step_gro_path = os.path.join(OUTPUT_PATH, "step.gro")
    step_edr_path = os.path.join(OUTPUT_PATH, "step.edr")
    step_log_path = os.path.join(OUTPUT_PATH, "step.log")
    step_xtc_path = os.path.join(OUTPUT_PATH, "step.xtc")
    step_cpt_path = os.path.join(OUTPUT_PATH, "step.cpt")

    # Working directory
    workingPath = os.getcwd()

    # log.err and log.out temporal file paths
    tmpErrsPath = Path(workingPath).rglob('log*.err')
    tmpOutsPath = Path(workingPath).rglob('log*.out')

    # Paths for last accepted gro and cpt files
    lastAcceptedCPT = os.path.join(OUTPUT_PATH, "lastAcceptedStep.cpt")
    lastAcceptedGRO = os.path.join(OUTPUT_PATH, "lastAcceptedStep.gro")

    # Save input GRO file as lastAcceptedGRO
    os.rename(args.input_path, lastAcceptedGRO)

    # Initialize 'within threshold' condition
    notWithinThreshold = True

    # Initialize 'last step accepted' condition
    lastStepAccepted = False

    main_wf(args)

    # NOTE: test with MPI + GPU - improve performance
    
    # NOTE: add limit to the number of failed steps from a given structure -> generate new starting conditions

    # NOTE: test with longer simulation times

    # NOTE: analyze the trajectory

    # NOTE: add improvements from papers: longer MD once within threshold

