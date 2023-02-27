##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import sys
import yaml
import time
import shutil
import logging
import argparse
import numpy as np
import mdtraj as md
from glob import glob
from scipy import stats
from pathlib import PurePath
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

def analyzeCV(MD_paths, Global_properties):
    """
    Analyze user-defined CV. The CV is defined as the distance between two atom selections from
    the YAML configuration. Fits a line to the evolution of the CV and returns the CV average value and the slope of the fitted line

    Inputs
    ------

    MD_paths          (dict): paths of input and output files of short MD
    Global_properties (dict): Atomic selections in Global_properties dictionary
    
    Output
    ------

    CV_mean  (float): average value of CV during short MD
    CV_slope (float): slope of the fitted line to the CV evolution 
    """

    logging.info('+++++++++ Analyzing CV...')

    # Read trajectory and topology
    step_traj = md.load(MD_paths["step_xtc_path"], top=MD_paths["step_gro_path"])

    # Find the COM of each group selection for each frame (nm)
    com_group1 = md.compute_center_of_mass(traj = step_traj, select = Global_properties['sumd']['colvar']['group1_selection'])
    com_group2 = md.compute_center_of_mass(traj = step_traj, select = Global_properties['sumd']['colvar']['group2_selection'])

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

    # Trajectory saving period in ps: time step (ps) * saving frequency (time steps)
    traj_freq = Global_properties['md']['properties']['mdp']['dt']*Global_properties['md']['properties']['mdp']['nstxout-compressed']

    # Find slope m in nm/ns
    CV_slope = CV_slope/traj_freq*1000

    logging.info('+++++++++ CV slope: {}'.format(CV_slope))

    return CV_mean, CV_slope

def acceptLastStep(MD_paths):
    """
    Accept last short MD simulation:

        - Add short trajectory to total trajectory
        - Set the step checkpoint file as the lastAcceptedCPT
        - Set the step gro file as the lastAcceptedGRO

    Inputs
    ------

    MD_paths   (dict): paths of input and output files of short MD

    Output
    ------

    (Implicitly)
    Files renamed and saved
    """

    logging.info('+++++++++ Accepting step!')

    # Save xtc in final trajectory
    if os.path.exists(MD_paths["traj_path"]):
        
        # remove previous zip trajectory bundle if it exists
        removeFiles(MD_paths["zip_path"])

        # Create a ZipFile object
        zipObject = ZipFile(MD_paths["zip_path"], 'w')
        
        # Add traj to zip file
        zipObject.write(MD_paths["traj_path"])
        zipObject.write(MD_paths["step_xtc_path"])

        # Close zip file
        zipObject.close()

        prop = { 'concatenate': True }

        # Concatenate
        trjcat(input_trj_zip_path=MD_paths["zip_path"],
                output_trj_path=MD_paths["traj_path"],
                properties=prop)

    else:

        # Change name
        os.rename(MD_paths["step_xtc_path"], MD_paths["traj_path"])

    # Save checkpoint
    os.rename(MD_paths["step_cpt_path"], MD_paths["last_cpt_path"])

    # Save GRO file
    os.rename(MD_paths["step_gro_path"], MD_paths["last_gro_path"])

    return

def removeTmpLogs(MD_paths):
    """
    Remove temporal log files created by biobb_gromacs

    Inputs
    ------

    MD_paths   (dict): paths of input and output files of short MD
    """

    errLogList = glob(MD_paths["std_err_path"])
    outLogList = glob(MD_paths["std_out_path"])

    for path in errLogList:
        removeFiles(path)
    
    for path in outLogList:
        removeFiles(path)

    return


def main_wf(args):
    '''
    Main workflow that takes as input the GROMACs topology and structure files together with a certain configuration as global variables to
    perform Supervised Molecular Dynamics. 

    The main idea is to run short MDs and track the evolution of a certain CV defined by the user. If the CV value approaches the 
    target value during the short MD simulation, then the final state of the system is accepted and used to perform the next short MD.
    However, if the CV gets away from the target value, then the next simulation re-starts from the last accepted state generating 
    new velocities to try to explore a different direction in conformational space.

    Inputs
    ------

    System structure (.gro)
    gromacs topology (.zip of .top)
    Configuration (different global variables)

    '''

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
            Global_properties = yaml.load(config_file, Loader = yaml.FullLoader)
    else:
        # Issue error
        logging.error("Configuration file not found!")
        logging.error("Check the configuration file exists: {}".format(args.configuration))
        sys.exit()

    # Dictionary with common MD simulation properties from input configuration
    MD_properties = Global_properties['md']['properties'].copy()

    # Working directory
    workingPath = os.getcwd()

    # Dictionary with paths
    MD_paths = {}

    # Path to trajectory file: if there are previous trajectories, add number to new trajectory to avoid overwriting
    MD_paths.update({"traj_path" : checkOutputExistence(OUTPUT_PATH, Global_properties['sumd']['files']['trajectory'])}) 
    # Path to zip file 
    MD_paths.update({"zip_path" : checkOutputExistence(OUTPUT_PATH, 'trajectory_bundle.zip')}) 
    # Names of step temporal files (short MD results)
    MD_paths.update({"step_tpr_path" : os.path.join(OUTPUT_PATH, "step.tpr")})
    MD_paths.update({"step_trr_path" : os.path.join(OUTPUT_PATH, "step.trr")})
    MD_paths.update({"step_gro_path" : os.path.join(OUTPUT_PATH, "step.gro")})
    MD_paths.update({"step_edr_path" : os.path.join(OUTPUT_PATH, "step.edr")})
    MD_paths.update({"step_log_path" : os.path.join(OUTPUT_PATH, "step.log")})
    MD_paths.update({"step_xtc_path" : os.path.join(OUTPUT_PATH, "step.xtc")})
    MD_paths.update({"step_cpt_path" : os.path.join(OUTPUT_PATH, "step.cpt")})
    # log.err and log.out temporal file paths
    MD_paths.update({"std_err_path" : os.path.join(workingPath, "log*.err")})
    MD_paths.update({"std_out_path" : os.path.join(workingPath, "log*.out")})
    # Paths for last accepted gro and cpt files
    MD_paths.update({"last_cpt_path" : os.path.join(OUTPUT_PATH, "lastAcceptedStep.cpt")})
    MD_paths.update({"last_gro_path" : os.path.join(OUTPUT_PATH, "lastAcceptedStep.gro")})
    # Path to topology
    MD_paths.update({"topology" : args.topology_path})

    # Save input GRO file as lastAcceptedGRO
    shutil.copyfile(args.input_path, MD_paths["last_gro_path"])

    # Initialize 'within threshold' condition
    notWithinThreshold = True

    # Initialize 'last step accepted' condition
    lastStepAccepted = False

    # Start timer
    start_time = time.time()

    # Initialize total step counter 
    step_counter = 0

    # Find max number of steps (short MDs) to do
    maxNumSteps = Global_properties['sumd']['num_steps']['max_total']

    # Find target value of CV
    CV_target = Global_properties['sumd']['colvar']['target']

    # Find threshold in the CV within which we start normal MD 
    CV_threshold = Global_properties['sumd']['colvar']['threshold']

    # Find threshold of slope fitting CV evolution
    slopeThreshold = Global_properties['sumd']['slope_threshold']

    # Find temperature of velocity generation
    T = Global_properties["md"]["properties"]["mdp"]["gen-temp"]

    # While all conditions True, keep running steps
    while (step_counter < maxNumSteps) and notWithinThreshold:

        logging.info('++++++ STEP {}'.format(step_counter))

        if lastStepAccepted:

            # If last step was accepted, continue MD with same velocities 
            logging.info('+++++++++ Continuation of previous step...')

            # Modify properties
            MD_properties['mdp'].update({'continuation' : 'yes'})
            MD_properties['mdp'].update({'gen-vel' : 'no'})
            _ = MD_properties['mdp'].pop('gen-temp', None)

            # Execute a short MD simulation, restarting from last checkpoint
            grompp_mdrun(input_gro_path=MD_paths["last_gro_path"],
                    input_top_zip_path=MD_paths["topology"],
                    input_cpt_path=MD_paths["last_cpt_path"],
                    output_tpr_path=MD_paths["step_tpr_path"],
                    output_trr_path=MD_paths["step_trr_path"],
                    output_gro_path=MD_paths["step_gro_path"],
                    output_edr_path=MD_paths["step_edr_path"],
                    output_log_path=MD_paths["step_log_path"],
                    output_xtc_path=MD_paths["step_xtc_path"],
                    output_cpt_path=MD_paths["step_cpt_path"],
                    properties=MD_properties)

            lastStepAccepted = False

        else:
            
            # Otherwise, restart from last accepted structure with new velocities 
            logging.info('+++++++++ Generating new velocities...')

            # Modify properties
            MD_properties['mdp'].update({'continuation' : 'no'})
            MD_properties['mdp'].update({'gen-vel' : 'yes'})
            MD_properties['mdp'].update({'gen-temp' : T})

            # Execute a short MD simulation, generating new velocities
            grompp_mdrun(input_gro_path=MD_paths["last_gro_path"],
                    input_top_zip_path=MD_paths["topology"],
                    output_tpr_path=MD_paths["step_tpr_path"],
                    output_trr_path=MD_paths["step_trr_path"],
                    output_gro_path=MD_paths["step_gro_path"],
                    output_edr_path=MD_paths["step_edr_path"],
                    output_log_path=MD_paths["step_log_path"],
                    output_xtc_path=MD_paths["step_xtc_path"],
                    output_cpt_path=MD_paths["step_cpt_path"],
                    properties=MD_properties)

        # Use MDTraj to analyze a certain CV (distance between 2 user-defined groups) 
        CV_mean, CV_slope = analyzeCV(MD_paths, Global_properties)

        # Check if distance to target is smaller than CV threshold
        if abs(CV_mean - CV_target) < CV_threshold:

            # We are within threshold!
            notWithinThreshold = False

            logging.info('+++++++++ CV is within threshold of target! :)')
        
        # Check slope is higher than CV slope threshold
        if abs(CV_slope) > slopeThreshold:

            # Accept new structure if CV_mean < CV_target and m > 0
            if (CV_mean < CV_target) and (CV_slope > 0):
                
                # Accept last step
                acceptLastStep(MD_paths)
 
                lastStepAccepted = True

            # Accept new structure if CV_mean > CV_target and m < 0
            elif (CV_mean > CV_target) and (CV_slope < 0): 
                
                # Accept last step 
                acceptLastStep(MD_paths)

                lastStepAccepted = True

        # Increase total counter
        step_counter += 1

        # Remove temporal log files
        removeTmpLogs(MD_paths)
    
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

    main_wf(args)

    # NOTE: test with MPI + GPU - improve performance. Then test with longer simulation times
    
    # NOTE: add limit to the number of failed steps from a given structure -> generate new starting conditions

    # NOTE: Add system preparation with AnteChamber / GROMACS and analysis of the trajectory. Recall that MMPBSA needs mol2 file from AnteChamber!!

    # github.com/Valdes-Tresanco-MS/gmx_MMPBSA

