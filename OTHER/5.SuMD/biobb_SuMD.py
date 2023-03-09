##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import time
import shutil
import argparse
import numpy as np
import mdtraj as md
from glob import glob
from scipy import stats
from zipfile import ZipFile

# Import biobb
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.grompp_mdrun import grompp_mdrun
from biobb_analysis.gromacs.gmx_image import gmx_image


def remove_files(*file_path):
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

def remove_tmp_files(*paths):
    """
    Remove temporal files created by biobb_gromacs if any

    Inputs
    ------

    paths   (str): variable number of paths to files including filename
    """

    tmp_files = []

    for path in paths:
        matching_files = glob(path)
        tmp_files.extend(matching_files)

    for file in tmp_files:
        os.remove(file)

    return

def get_tmp_folders(a_dir):
    """
    Temporal fix!
    Get all the unique temporal directories created by grompp_mdrun. Exclude all other directories
    """

    all_subdirs = [name for name in os.listdir(a_dir) if os.path.isdir(os.path.join(a_dir, name))]

    paths = []

    for name in all_subdirs:
        
        # Remove all dirs with short names
        if len(name) < 20:
            all_subdirs.remove(name)

        # Convert the rest into paths
        else:
            paths.append(os.path.join(a_dir, name))

    return paths

def get_tpr_file(paths):
    """
    Temporal fix!
    Finds the internal.tpr file in the temporal directories created by grompp_mdrun
    Returns the path to that file or None if it does not exist
    """

    for path in paths:
        matching_files = glob(os.path.join(path, 'internal.tpr'))
        if matching_files:
            return matching_files[0]
    
    return None

def main_wf(configuration_path, input_structure, topology_path, index_path):
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
    index_path         (str): path to GROMACS index file (.ndx)
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
    continuation_MD_prop['mdp']['continuation'] = 'yes'
    continuation_MD_prop['mdp']['gen-vel'] = 'no'
    _ = continuation_MD_prop['mdp'].pop('gen-temp', None)

    # Paths and properties of trajectory imaging step
    traj_imaging_paths = global_paths['trajectory_imaging'].copy()
    traj_imaging_prop = global_prop['trajectory_imaging'].copy()

    # Paths and properties of trajectory fitting step
    traj_fitting_paths = global_paths['trajectory_fitting'].copy()
    traj_fitting_prop = global_prop['trajectory_fitting'].copy()

    # Add index path to paths if provided
    if index_path is not None:
        new_vel_MD_paths['input_ndx_path'] = index_path
        continuation_MD_paths['input_ndx_path'] = index_path
        traj_imaging_paths['input_index_path'] = index_path
        traj_fitting_paths['input_index_path'] = index_path

    # Paths and properties of concatenation step 
    concat_prop = global_prop["trajectory_cat"]
    concat_paths = global_paths["trajectory_cat"] 

    # Initialize a SuMD log file
    sumd_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True, prefix='sumd')

    # Get SuMD properties
    sumd_prop = conf.properties['sumd']

    # Create short MD folder -> to copy input structure and topology inside before step execution
    if not os.path.exists(new_vel_MD_prop["path"]):
        os.makedirs(new_vel_MD_prop["path"])

    # Create concatenation folder -> to put Zip file inside before step execution
    if not os.path.exists(concat_prop["path"]):
        os.makedirs(concat_prop["path"])

    # Create gmx imaging folder -> temporal fix 
    if not os.path.exists(traj_imaging_prop["path"]):
        os.makedirs(traj_imaging_prop["path"])
    
    # Copy input structure and topology to short MD folder
    shutil.copy(input_structure, continuation_MD_paths["input_gro_path"])
    shutil.copy(topology_path, continuation_MD_paths["input_top_zip_path"])

    # Initialize 'last step accepted' condition
    last_step_accepted = False

    # Initialize total step counter 
    step_counter = 0

    # Failed steps counter
    failed_steps = 0

    # Find max number of steps (short MDs) to do
    max_num_steps = sumd_prop['num_steps']['max_total']

    # Find target value of CV
    cv_target = sumd_prop['colvar']['target']

    # Find threshold in the CV within which we start normal MD 
    cv_threshold = sumd_prop['colvar']['threshold']

    # Find threshold for the slope fitting CV evolution
    slope_threshold = sumd_prop['slope_threshold']

    # Initialize list for temporal folders with tpr files
    tmp_folders = []

    while (step_counter < max_num_steps):

        sumd_log.info('STEP {}'.format(step_counter))

        if last_step_accepted or (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            if failed_steps > sumd_prop['num_steps']['max_failed']:
                sumd_log.info('  Maximum number of failed steps reached, running non-supervised MD :(')
                
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

            sumd_log.info('  CV is within threshold of target :)')
            sumd_log.info('  Accepting step')
            last_step_accepted = True
        
        # Check if slope is larger than CV slope's threshold
        if abs(cv_slope) > slope_threshold:
            
            # Check if CV is approaching target
            cv_increasing_towards_target = (cv_mean < cv_target) and (cv_slope > 0)
            cv_decreasing_towards_target = (cv_mean > cv_target) and (cv_slope < 0)

            if cv_increasing_towards_target or cv_decreasing_towards_target:

                sumd_log.info('  CV is approaching target :)')
                sumd_log.info('  Accepting step')
                last_step_accepted = True
        
        if last_step_accepted:

            # First accepted step
            if not os.path.exists(concat_paths["output_trj_path"]):
                
                # Copy short MD trajectory to concatenation folder
                shutil.copyfile(continuation_MD_paths["output_xtc_path"], concat_paths["output_trj_path"])
            
            # Subsequent accepted steps
            else:

                # remove previous zip trajectory bundle if it exists
                remove_files(concat_paths["input_trj_zip_path"])

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

            # Reset failed steps counter
            failed_steps = 0

        else:

            # Increase failed steps counter
            failed_steps += 1

        if last_step_accepted or (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            # Replace last accepted structure and checkpoint for new ones
            shutil.copyfile(continuation_MD_paths["output_gro_path"], continuation_MD_paths["input_gro_path"])
            shutil.copyfile(continuation_MD_paths["output_cpt_path"], continuation_MD_paths["input_cpt_path"])

        # Increase total counter
        step_counter += 1

        # Find generic paths to log.err and log.out files
        short_MD_stdout_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.out")
        short_MD_stderr_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.err")
        concat_stdout_path = os.path.join(concat_prop["path"], concat_prop["step"] + "_log*.out")
        concat_stderr_path = os.path.join(concat_prop["path"], concat_prop["step"] + "_log*.err")

        # Remove temporal log files from grompp_mdrun and trjcat 
        remove_tmp_files(short_MD_stdout_path, short_MD_stderr_path, concat_stdout_path, concat_stderr_path)
    
    # Get all temporal unique folders with internal.tpr from grompp_mdrun 
    tmp_folders=get_tmp_folders(os.getcwd())

    # Find internal.tpr file 
    tpr_path = get_tpr_file(tmp_folders)

    # Copy internal.tpr file to input top path
    shutil.copy(tpr_path, traj_imaging_paths["input_top_path"])

    # Image the trajectory 
    gmx_image(**traj_imaging_paths, properties=traj_imaging_prop)

    # Fit the trajectory
    gmx_image(**traj_fitting_paths, properties=traj_fitting_prop)

    # Move all temporal unique folders with internal.tpr from grompp_mdrun to common folder
    for tmp_folder in tmp_folders:
        shutil.move(tmp_folder, os.path.join(os.getcwd(), "tmp"))

    # Print timing information to the log file
    elapsed_time = time.time() - start_time

    sumd_log.info('Elapsed time: {} minutes'.format(round(elapsed_time/60, 1)))

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('-structure', dest='structure_path',
                        help="Path to input structure (.gro)", 
                        required=True)

    parser.add_argument('-topology', dest='topology_path',
                        help="Path to compressed input topology (.zip of .top and .itp files)", 
                        required=True)

    parser.add_argument('-config', dest='config_path',
                        help="Path to configuration file (YAML)", 
                        required=True)
    
    parser.add_argument('-index', dest='index_path',
                        help="Path to index file (.ndx)", 
                        required=False)

    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, input_structure = args.structure_path, topology_path = args.topology_path, index_path = args.index_path)

    # NOTE: one could estimate the CV fluctuations from the first step, and adjust the CV's slope threshold and the output frequency accordingly

    # NOTE: refine the criteria for accepting a step: the slope with a few points is an approximation, we could check if the last value of the CV is closer to the target or not (taking into account the CV fluctuations) 

    # NOTE: analysis, github.com/Valdes-Tresanco-MS/gmx_MMPBSA or SuMD-analyzer
