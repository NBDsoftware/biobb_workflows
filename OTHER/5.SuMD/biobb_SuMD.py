##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import sys
import yaml
import time
import shutil
import argparse
import numpy as np
import mdtraj as md
from glob import glob
from scipy import stats
from pathlib import PurePath
from zipfile import ZipFile

# Import biobb
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.grompp_mdrun import grompp_mdrun

def checkFilesExist(*file_path, logger):
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

            logger.warning(" {} file was not found".format(path))
            
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

def analyze_cv(short_MD_paths, short_MD_prop, sumd_prop, sumd_log):
    """
    Analyze user-defined CV. The CV is defined as the distance between two atom selections from
    the YAML configuration. Fits a line to the evolution of the CV and returns the CV average value and the slope of the fitted line

    Inputs
    ------

    short_MD_paths (dict): dictionary with paths to short MD files
    short_MD_prop  (dict): dictionary with properties of the short MD protocol
    sumd_prop      (dict): dictionary with properties of the SuMD protocol
    sumd_log       (obj):  logger for the SuMD protocol
    
    Output
    ------

    cv_mean  (float): average value of CV during short MD
    cv_slope (float): slope of the fitted line to the CV evolution 
    """

    sumd_log.info('  Analyzing CV...')

    # Read trajectory and topology
    step_traj = md.load(short_MD_paths["output_xtc_path"], top=short_MD_paths["output_gro_path"])

    # Find the COM of each group selection for each frame (nm)
    com_group1 = md.compute_center_of_mass(traj = step_traj, select = sumd_prop['colvar']['group1_selection'])
    com_group2 = md.compute_center_of_mass(traj = step_traj, select = sumd_prop['colvar']['group2_selection'])

    # Find the vector joining the COMs for each frame
    distance_vectors = com_group1 - com_group2

    # Find the norm of this vector for each frame (nm)
    CV = np.linalg.norm(distance_vectors, axis=1)

    sumd_log.info('  CV: {}'.format(CV))

    # Time array (index of frames - arbitrary units)
    time_au = np.arange(0, len(CV))

    # Mean value of the CV
    cv_mean = np.mean(CV)

    sumd_log.info('  CV mean: {}'.format(cv_mean))

    # Compute linear regression 
    res = stats.linregress(time_au, CV)

    # Find slope in nm/frame
    cv_slope = res.slope

    # Trajectory saving period in ns: time step (ps) * saving frequency (time steps) / 1000 (ps/ns)
    saving_period = short_MD_prop['mdp']['dt']*short_MD_prop['mdp']['nstxout-compressed']/1000

    # Find slope in nm/ns
    cv_slope = cv_slope/saving_period

    sumd_log.info('  CV slope: {}'.format(cv_slope))

    return cv_mean, cv_slope

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

def remove_tmp_logs(MD_paths):
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


def main_wf(configuration_path, input_structure, topology_path, output_path):
    '''
    Main workflow that takes as input the GROMACs topology and structure files together with a certain configuration as global variables to
    perform Supervised Molecular Dynamics. 

    The main idea is to run short MDs and track the evolution of a certain CV defined by the user. If the CV value approaches the 
    target value during the short MD simulation, then the final state of the system is accepted and used to perform the next short MD.
    However, if the CV gets away from the target value, then the next simulation re-starts from the last accepted state generating 
    new velocities to try to explore a different direction in conformational space.

    Inputs
    ------

    configuration_path (str): path to YAML configuration file
    input_structure    (str): path to GROMACS structure file (.gro)
    topology_path      (str): path to GROMACS topology file (.zip of .top)
    output_path        (str): path to output folder
    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initialize a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True, prefix='global')

    # Get global properties and paths 
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Paths and Properties of short MD (generating new velocities)
    new_vel_MD_paths = global_paths['short_MD'].copy()
    new_vel_MD_prop = global_prop['short_MD'].copy()
    _ = new_vel_MD_paths.pop('input_cpt_path', None)

    # Paths and Properties of short MD (continuation from last accepted state)
    continuation_MD_paths = global_paths['short_MD'].copy()
    continuation_MD_prop = global_prop['short_MD'].copy()
    continuation_MD_prop['mdp'].update({'continuation' : 'yes'})
    continuation_MD_prop['mdp'].update({'gen-vel' : 'no'})
    _ = continuation_MD_prop['mdp'].pop('gen-temp', None)

    # Properties and paths of concatenation step 
    concat_prop = global_prop["trajectory_cat"]
    concat_paths = global_paths["trajectory_cat"] 

    # Initialize a SuMD log file
    sumd_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True, prefix='sumd')

    # Get SuMD properties
    sumd_prop = conf.properties['sumd']

    if output_path is None:
        # Default output path
        output_path = conf.get_working_dir_path()

    # Create short MD folder
    if not os.path.exists(new_vel_MD_prop["path"]):
        os.makedirs(new_vel_MD_prop["path"])

    # Create concatenation folder
    if not os.path.exists(concat_prop["path"]):
        os.makedirs(concat_prop["path"])
    
    # Copy input structure and topology to short MD folder
    shutil.copy(input_structure, continuation_MD_paths["input_gro_path"])
    shutil.copy(topology_path, continuation_MD_paths["input_top_zip_path"])

    # Initialize 'cv within threshold' condition
    within_threshold = False

    # Initialize 'last step accepted' condition
    last_step_accepted = False

    # Initialize total step counter 
    step_counter = 0

    # Find max number of steps (short MDs) to do
    max_num_steps = sumd_prop['num_steps']['max_total']

    # Find target value of CV
    cv_target = sumd_prop['colvar']['target']

    # Find threshold in the CV within which we start normal MD 
    cv_threshold = sumd_prop['colvar']['threshold']

    # Find threshold for the slope fitting CV evolution
    slope_threshold = sumd_prop['slope_threshold']

    # While all conditions True, keep running steps
    while (step_counter < max_num_steps) and not within_threshold:

        sumd_log.info('STEP {}'.format(step_counter))

        if last_step_accepted:

            # Continue MD with same velocities 
            sumd_log.info('  Continuation of previous step...')

            # Execute a short MD simulation, restarting from last checkpoint
            grompp_mdrun(**continuation_MD_paths, properties=continuation_MD_prop)

            last_step_accepted = False

        else:
            
            # Restart from last accepted structure with new velocities 
            sumd_log.info('  Generating new velocities...')

            # Execute a short MD simulation, generating new velocities
            grompp_mdrun(**new_vel_MD_paths, properties=new_vel_MD_prop)

        # Use MDTraj to analyze the CV (distance between 2 user-defined groups) 
        cv_mean, cv_slope = analyze_cv(continuation_MD_paths, continuation_MD_prop, sumd_prop, sumd_log)

        # Check if distance to target is smaller than CV threshold
        if abs(cv_mean - cv_target) < cv_threshold:

            sumd_log.info('  CV is within threshold of target! :)')

            within_threshold = True
        
        # Check if slope is higher than CV slope threshold
        if abs(cv_slope) > slope_threshold:

            cv_increasing_towards_target = (cv_mean < cv_target) and (cv_slope > 0)
            cv_decreasing_towards_target = (cv_mean > cv_target) and (cv_slope < 0)

            # Accept last step
            if cv_increasing_towards_target or cv_decreasing_towards_target:
                
                sumd_log.info('  Accepting step!')
                
                # First accepted step
                if not os.path.exists(concat_paths["output_trj_path"]):
                    
                    # Copy short MD trajectory to concatenation folder
                    shutil.copyfile(continuation_MD_paths["output_xtc_path"], concat_paths["output_trj_path"])
                
                # Subsequent accepted steps
                else:

                    # remove previous zip trajectory bundle if it exists
                    removeFiles(concat_paths["input_trj_zip_path"])

                    # Create a ZipFile object
                    zipObject = ZipFile(concat_paths["input_trj_zip_path"], 'w')

                    # Add previously concatenated trajectory to zip file
                    zipObject.write(concat_paths["output_trj_path"])

                    # Add short MD trajectory to zip file
                    zipObject.write(continuation_MD_paths["output_xtc_path"])

                    # Close zip file
                    zipObject.close()

                    # Concatenate
                    trjcat(**concat_paths, properties=concat_prop)

                # Replace last accepted structure and checkpoint
                shutil.copyfile(continuation_MD_paths["output_gro_path"], continuation_MD_paths["input_gro_path"])
                shutil.copyfile(continuation_MD_paths["output_cpt_path"], continuation_MD_paths["input_cpt_path"])

                last_step_accepted = True

        # Increase total counter
        step_counter += 1

        # Remove temporal log files
        # remove_tmp_logs(MD_paths)

        """
        # log.err and log.out temporal file paths
        MD_paths.update({"std_err_path" : os.path.join(conf.get_working_dir_path(), "log*.err")})
        MD_paths.update({"std_out_path" : os.path.join(conf.get_working_dir_path(), "log*.out")})
        """
    
    # Move final trajectory to output path
    shutil.move(concat_paths["output_trj_path"], output_path)

    # Print timing information to the log file
    elapsed_time = time.time() - start_time

    sumd_log.info('Elapsed time: {} minutes'.format(round(elapsed_time/60, 1)))

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('-input', dest='input_path',
                        help="Path to input structure (.gro)", 
                        required=True)

    parser.add_argument('-topology', dest='topology_path',
                        help="Path to compressed input topology (.zip of .top and .itp files)", 
                        required=True)

    parser.add_argument('-config', dest='config_path',
                        help="Path to configuration file (YAML)", 
                        required=True)

    parser.add_argument('-output', dest='output_path',
                        help="Path where results will be dumped (./output by default)", 
                        required=False)

    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, input_structure = args.input_path, topology_path = args.topology_path, output_path = args.output_path)

    # NOTE: test with MPI + GPU - improve performance. Then test with longer simulation times
    
    # NOTE: add limit to the number of failed steps from a given structure -> generate new starting conditions

    # NOTE: Add system preparation with AnteChamber / GROMACS and analysis of the trajectory. Recall that MMPBSA needs mol2 file from AnteChamber!!

    # github.com/Valdes-Tresanco-MS/gmx_MMPBSA
