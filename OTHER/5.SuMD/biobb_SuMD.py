##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import csv
import time
import shutil
import argparse
import numpy as np
from glob import glob
from scipy import stats
import MDAnalysis as mda
from zipfile import ZipFile
from MDAnalysis.lib.distances import distance_array

# Import biobb
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.grompp_mdrun import grompp_mdrun
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str


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

def analyze_cv(traj_path, topology_path, MD_properties, accepted_steps, sumd_properties, sumd_log, output_path):
    """
    Analyze user-defined CV. The CV is defined as the distance between two atom selections from
    the YAML configuration. Fits a line to the evolution of the CV, computes the slope and the  
    mean value of the CV. If the CV is within a threshold of the target value or if the CV is
    approaching the target value, the step is accepted.

    Inputs
    ------

    traj_path        (str):  path to the trajectory file
    topology_path    (str):  path to the topology file
    MD_properties   (dict):  dictionary with properties of the short MD protocol
    accepted_steps   (int):  number of previously accepted steps
    sumd_properties (dict):  dictionary with properties of the SuMD protocol
    sumd_log         (obj):  logger for the SuMD protocol
    output_path      (str):  path to the output folder to save the CV time series file
    
    Output
    ------

    accept_step (bool): True if the step is accepted, False otherwise
    """

    # Initialize accept_step to False
    accept_step = False

    # Find target value of CV
    cv_target = sumd_properties['colvar']['target']

    # Find threshold in the CV within which we start normal MD 
    cv_threshold = sumd_properties['colvar']['threshold']

    # Find threshold for the slope fitting CV evolution
    slope_threshold = sumd_properties['slope_threshold']

    sumd_log.info('  Analyzing CV...')

    # Read trajectory and topology
    mda_universe = mda.Universe(topology_path, traj_path)

    # Select the atoms
    group_1_selection = mda_universe.select_atoms(sumd_properties['colvar']['group1_selection'])
    group_2_selection = mda_universe.select_atoms(sumd_properties['colvar']['group2_selection'])

    # Calculate the distance between them taking into account periodic boundary conditions (box)
    cv_time_series = []
    for ts in mda_universe.trajectory: 
        distance = distance_array(group_1_selection.center_of_mass(), group_2_selection.center_of_mass(), box=mda_universe.dimensions)
        cv_time_series.append(distance[0][0])
    
    # Convert to nm
    cv_time_series = np.array(cv_time_series)/10

    sumd_log.info('  CV time series: {}'.format(cv_time_series))

    # Time array (index of frames)
    time = np.arange(0, len(cv_time_series))

    # Trajectory saving period in ns: time step (ps) * saving frequency in time steps / 1000 
    saving_period = MD_properties['dt']*MD_properties['nstxout-compressed']/1000

    # Time array in ns
    time = time*saving_period

    # Find length of short MD trajectory
    short_MD_length = (len(cv_time_series)-1)*saving_period

    # Add time of previously accepted steps
    time = time + accepted_steps*short_MD_length

    # Mean value of the CV
    cv_mean = np.mean(cv_time_series)

    sumd_log.info('  CV mean: {}'.format(cv_mean))

    # Compute linear regression 
    res = stats.linregress(time, cv_time_series)

    # Find slope in nm/ns
    cv_slope = res.slope

    sumd_log.info('  CV slope: {}'.format(cv_slope))

    # Find if the last CV value is closer than the first to the target value

    # Check if CV is within threshold of target
    if (abs(cv_mean - cv_target) < cv_threshold) and (abs(cv_time_series[-1] - cv_target) < cv_threshold):

        sumd_log.info('  CV is within threshold of target :)')
        sumd_log.info('  Accepting step')
        accept_step = True
    
    # Check if CV is approaching target
    if abs(cv_slope) > slope_threshold:
        
        # Slope sign 
        positive_slope = cv_slope > 0
        negative_slope = cv_slope < 0

        # Relative position of CV to target value 
        cv_below_target = cv_mean < cv_target
        cv_above_target = cv_mean > cv_target

        # CV evolution
        cv_increasing = cv_time_series[0] < cv_time_series[-1]
        cv_decreasing = cv_time_series[0] > cv_time_series[-1]

        if (cv_below_target and positive_slope and cv_increasing) or (cv_above_target and negative_slope and cv_decreasing):

            sumd_log.info('  CV is approaching target')
            sumd_log.info('  Accepting step')
            accept_step = True
    
    if not accept_step:
        sumd_log.info('  CV is not approaching target')
        sumd_log.info('  Rejecting step')

    if accept_step:
        # Append CV data to CV evolution file
        save_cv_data(cv_time_series, time, os.path.join(output_path, 'cv_evolution.dat'))

    return accept_step

def save_cv_data(cv_time_series, time, path):
    """
    Save CV time series and time array to a file
    """

    # combine time and values into a list of tuples
    data = list(zip(time, cv_time_series))

    # check if file exists
    if os.path.isfile(path):
        # append data to existing file
        with open(path, mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)
    else:
        # create new file and write header and data
        with open(path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['time (ns)', 'CV (nm)'])
            writer.writerows(data)

    return

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
    Temporal fix! - Deprecated
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

    # Paths and properties of short MD imaging step
    imaging_MD_paths = global_paths['short_MD_image'].copy()
    imaging_MD_prop = global_prop['short_MD_image'].copy()

    # Paths and properties of short MD structure extraction step
    str_MD_paths = global_paths['short_MD_str'].copy()
    str_MD_prop = global_prop['short_MD_str'].copy()

    # Add index path to steps if provided
    if index_path is not None:
        new_vel_MD_paths['input_ndx_path'] = index_path
        continuation_MD_paths['input_ndx_path'] = index_path
        imaging_MD_paths['input_index_path'] = index_path
        str_MD_paths['input_index_path'] = index_path

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

    # Copy input structure and topology to short MD folder
    shutil.copy(input_structure, continuation_MD_paths["input_gro_path"])
    shutil.copy(topology_path, continuation_MD_paths["input_top_zip_path"])

    # Create concatenation folder -> to put Zip file inside before step execution
    if not os.path.exists(concat_prop["path"]):
        os.makedirs(concat_prop["path"])

    # Create imaging folder -> temporal fix 
    if not os.path.exists(imaging_MD_prop["path"]):
        os.makedirs(imaging_MD_prop["path"])
    
    # Create structure extraction folder -> temporal fix 
    if not os.path.exists(str_MD_prop["path"]):
        os.makedirs(str_MD_prop["path"])

    # Initialize 'last step accepted' condition
    last_step_accepted = False

    # Initialize counters
    total_steps, failed_steps, accepted_steps = 0,0,0

    # Find max number of steps (short MDs) to do
    max_num_steps = sumd_prop['num_steps']['max_total']

    # Initialize list for temporal folders with tpr files -> will be fixed in next grompp_mdrun release
    tmp_folders = []

    while (total_steps < max_num_steps):

        sumd_log.info('STEP {}'.format(total_steps))

        if last_step_accepted:
                
            # Continue MD with same velocities 
            sumd_log.info('  Continuation of previous step...')
            # Execute a short MD simulation, restarting from last checkpoint
            grompp_mdrun(**continuation_MD_paths, properties=continuation_MD_prop)
        
        elif (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            # Continue MD with same velocities 
            sumd_log.info('  Maximum number of failed steps reached, continuation of previous step...')
            # Execute a short MD simulation, restarting from last checkpoint
            grompp_mdrun(**continuation_MD_paths, properties=continuation_MD_prop)

        else:
            
            # Restart from last accepted structure with new velocities 
            sumd_log.info('  Generating new velocities...')
            # Execute a short MD simulation, generating new velocities
            grompp_mdrun(**new_vel_MD_paths, properties=new_vel_MD_prop)

        # Image short MD trajectory (to avoid PBC issues when analyzing the CV)
        gmx_image(**imaging_MD_paths, properties=imaging_MD_prop)

        # Extract structure from short MD trajectory
        gmx_trjconv_str(**str_MD_paths, properties=str_MD_prop)

        # NOTE: which topology should we use? the step.gro? or a .gro from the imaged trajectory?

        # Use MDAnalysis to analyze the CV (distance between 2 user-defined groups)
        last_step_accepted = analyze_cv(imaging_MD_paths["output_traj_path"], str_MD_paths["output_str_path"], 
                                        continuation_MD_prop['mdp'], accepted_steps, sumd_prop, sumd_log, output_path = conf.get_working_dir_path())
        
        if last_step_accepted:

            # First accepted step
            if not os.path.exists(concat_paths["output_trj_path"]):
                
                # Copy short MD trajectory to concatenation folder
                shutil.copyfile(imaging_MD_paths["output_traj_path"], concat_paths["output_trj_path"])
            
            # Subsequent accepted steps
            else:

                # remove previous zip trajectory bundle if it exists
                remove_files(concat_paths["input_trj_zip_path"])

                # Create a ZipFile object
                zipObject = ZipFile(concat_paths["input_trj_zip_path"], 'w')

                # Add previously concatenated trajectory to zip file
                zipObject.write(concat_paths["output_trj_path"])

                # Add short MD trajectory to zip file
                zipObject.write(imaging_MD_paths["output_traj_path"])

                # Close zip file
                zipObject.close()

                # Concatenate
                trjcat(**concat_paths, properties=concat_prop)

            # Reset failed steps counter
            failed_steps = 0

            # Increase accepted steps counter
            accepted_steps += 1

        else:

            # Increase failed steps counter
            failed_steps += 1

        if last_step_accepted or (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            # Replace last accepted structure and checkpoint for newly accepted ones
            shutil.copyfile(continuation_MD_paths["output_gro_path"], continuation_MD_paths["input_gro_path"])
            shutil.copyfile(continuation_MD_paths["output_cpt_path"], continuation_MD_paths["input_cpt_path"])

        # Increase total counter
        total_steps += 1

        # Find generic paths to temporal log files from steps
        short_MD_stdout_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.out")
        short_MD_stderr_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.err")
        image_MD_stdout_path = os.path.join(imaging_MD_prop["path"], imaging_MD_prop["step"] + "_log*.out")
        image_MD_stderr_path = os.path.join(imaging_MD_prop["path"], imaging_MD_prop["step"] + "_log*.err")
        concat_stdout_path = os.path.join(concat_prop["path"], concat_prop["step"] + "_log*.out")
        concat_stderr_path = os.path.join(concat_prop["path"], concat_prop["step"] + "_log*.err")

        # Remove temporal log files from steps
        remove_tmp_files(short_MD_stdout_path, short_MD_stderr_path,
                         image_MD_stdout_path, image_MD_stderr_path, 
                         concat_stdout_path, concat_stderr_path)
    
    # Get all temporal unique folders with internal.tpr from grompp_mdrun -> will be fixed in next grompp_mdrun release
    tmp_folders=get_tmp_folders(os.getcwd())

    # Move all temporal unique folders with internal.tpr from grompp_mdrun to common folder -> will be fixed in next grompp_mdrun release
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

    # TODO: estimate CV's local fluctuations and relaxation time from the first step, and adjust the CV's slope threshold and the output frequency accordingly
    # TODO: as we are keeping only the dry trajectory - we could increase the output frequency of the MD simulation
    # TODO: add biasing potential to RC - several approaches (start from the most simple)
