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
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.mdrun import mdrun
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj

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

    # Round time and CV values to 3 decimal places
    time = np.round(time, 3)
    cv_time_series = np.round(cv_time_series, 3)

    # combine time and values into a list of tuples
    data = list(zip(time, cv_time_series))

    # check if file exists
    if os.path.isfile(path):
        # append data to existing file
        with open(path, mode='a', newline='') as file:
            # Append data using whitespace as separator 
            writer = csv.writer(file, delimiter=' ')
            writer.writerows(data)
    else:
        # create new file and write header and data
        with open(path, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=' ')
            writer.writerow(['time (ns)', 'CV (nm)'])
            writer.writerows(data)

    return

def remove_logs(continuation_MD_prop, global_prop):
    """
    Remove all log files in previous steps
    """

    # Find generic paths to log files from steps
    short_MD_stdout_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.out")
    short_MD_stderr_path = os.path.join(continuation_MD_prop["path"], continuation_MD_prop["step"] + "_log*.err")
    mdrun_stdout_path = os.path.join(global_prop["short_MD_mdrun"]["path"], global_prop["short_MD_mdrun"]["step"] + "_log*.out")
    mdrun_stderr_path = os.path.join(global_prop["short_MD_mdrun"]["path"], global_prop["short_MD_mdrun"]["step"] + "_log*.err")
    pbc_1_stdout_path = os.path.join(global_prop["pbc_1_whole"]["path"], global_prop["pbc_1_whole"]["step"] + "_log*.out")
    pbc_1_stderr_path = os.path.join(global_prop["pbc_1_whole"]["path"], global_prop["pbc_1_whole"]["step"] + "_log*.err")
    pbc_2_stdout_path = os.path.join(global_prop["pbc_2_cluster"]["path"], global_prop["pbc_2_cluster"]["step"] + "_log*.out")
    pbc_2_stderr_path = os.path.join(global_prop["pbc_2_cluster"]["path"], global_prop["pbc_2_cluster"]["step"] + "_log*.err")
    pbc_3_stdout_path = os.path.join(global_prop["pbc_3_extract_frame"]["path"], global_prop["pbc_3_extract_frame"]["step"] + "_log*.out")
    pbc_3_stderr_path = os.path.join(global_prop["pbc_3_extract_frame"]["path"], global_prop["pbc_3_extract_frame"]["step"] + "_log*.err")
    pbc_4_stdout_path = os.path.join(global_prop["pbc_4_nojump"]["path"], global_prop["pbc_4_nojump"]["step"] + "_log*.out")
    pbc_4_stderr_path = os.path.join(global_prop["pbc_4_nojump"]["path"], global_prop["pbc_4_nojump"]["step"] + "_log*.err")
    pbc_5_stdout_path = os.path.join(global_prop["pbc_5_center"]["path"], global_prop["pbc_5_center"]["step"] + "_log*.out")
    pbc_5_stderr_path = os.path.join(global_prop["pbc_5_center"]["path"], global_prop["pbc_5_center"]["step"] + "_log*.err")
    pbc_6_stdout_path = os.path.join(global_prop["pbc_6_fit"]["path"], global_prop["pbc_6_fit"]["step"] + "_log*.out")
    pbc_6_stderr_path = os.path.join(global_prop["pbc_6_fit"]["path"], global_prop["pbc_6_fit"]["step"] + "_log*.err")
    concat_stdout_path = os.path.join(global_prop["original_trajectory_cat"]["path"], global_prop["original_trajectory_cat"]["step"] + "_log*.out")
    concat_stderr_path = os.path.join(global_prop["original_trajectory_cat"]["path"], global_prop["original_trajectory_cat"]["step"] + "_log*.err")

    # Remove log files
    remove_tmp_files(short_MD_stdout_path, short_MD_stderr_path,
                        mdrun_stdout_path, mdrun_stderr_path, 
                        pbc_1_stdout_path, pbc_1_stderr_path,
                        pbc_2_stdout_path, pbc_2_stderr_path,
                        pbc_3_stdout_path, pbc_3_stderr_path,
                        pbc_4_stdout_path, pbc_4_stderr_path,
                        pbc_5_stdout_path, pbc_5_stderr_path,
                        pbc_6_stdout_path, pbc_6_stderr_path,
                        concat_stdout_path, concat_stderr_path)

    return

def remove_tmp_files(*paths):
    """
    Remove temporal files if they exist

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
    Get all the unique temporal directories created by biobb. Exclude all other directories
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
    new_MD_paths = global_paths['short_MD_grompp'].copy()
    new_MD_prop = global_prop['short_MD_grompp'].copy()
    _ = new_MD_paths.pop('input_cpt_path', None)

    # Paths and Properties of short MD (continuation from last accepted state)
    continuation_MD_paths = global_paths['short_MD_grompp'].copy()
    continuation_MD_prop = global_prop['short_MD_grompp'].copy()
    continuation_MD_prop['mdp']['continuation'] = 'yes'
    continuation_MD_prop['mdp']['gen-vel'] = 'no'
    _ = continuation_MD_prop['mdp'].pop('gen-temp', None)

    # Add index path to steps if provided
    if index_path is not None:
        global_paths['dry_structure']['input_index_path'] = index_path
        new_MD_paths['input_ndx_path'] = index_path
        continuation_MD_paths['input_ndx_path'] = index_path
        global_paths['pbc_1_whole']['input_index_path'] = index_path
        global_paths['pbc_2_cluster']['input_index_path'] = index_path
        global_paths['pbc_3_extract_frame']['input_index_path'] = index_path
        global_paths['pbc_4_nojump']['input_index_path'] = index_path
        global_paths['pbc_5_center']['input_index_path'] = index_path
        global_paths['pbc_6_fit']['input_index_path'] = index_path

    # Initialize a SuMD log file
    sumd_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True, prefix='sumd')

    # Get SuMD properties
    sumd_prop = conf.properties['sumd']

    # Create short MD grompp folder -> to copy input structure and topology inside before first step
    if not os.path.exists(new_MD_prop["path"]):
        os.makedirs(new_MD_prop["path"])

    # Copy input structure and topology to short MD grompp folder
    shutil.copy(input_structure, continuation_MD_paths["input_gro_path"])
    shutil.copy(topology_path, continuation_MD_paths["input_top_zip_path"])

    # Create concatenation folders 
    if not os.path.exists(global_prop["original_trajectory_cat"]["path"]):
        os.makedirs(global_prop["original_trajectory_cat"]["path"])

    # Initialize 'last step accepted' condition
    last_step_accepted = False

    # Initialize step counters
    total_steps, failed_steps, accepted_steps = 0,0,0

    # Find max number of steps (short MDs) to do
    max_num_steps = sumd_prop['num_steps']['max_total']

    # Extract dry structure from input structure
    global_paths['dry_structure']['input_structure_path'] = input_structure
    global_paths['dry_structure']['input_top_path'] = input_structure
    gmx_trjconv_str(**global_paths['dry_structure'], properties=global_prop['dry_structure'])

    while (total_steps < max_num_steps):

        sumd_log.info('STEP {}'.format(total_steps))

        if last_step_accepted:
                
            # Continue MD with same velocities 
            sumd_log.info('  Continuation of previous step...')
            grompp(**continuation_MD_paths, properties=continuation_MD_prop)
            mdrun(**global_paths['short_MD_mdrun'], properties=global_prop['short_MD_mdrun'])
        
        elif (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            # Continue MD with same velocities 
            sumd_log.info('  Maximum number of failed steps reached, continuation of previous step...')
            grompp(**continuation_MD_paths, properties=continuation_MD_prop)
            mdrun(**global_paths['short_MD_mdrun'], properties=global_prop['short_MD_mdrun'])

        else:
            
            # Restart from last accepted structure with new velocities 
            sumd_log.info('  Generating new velocities...')
            grompp(**new_MD_paths, properties=new_MD_prop)
            mdrun(**global_paths['short_MD_mdrun'], properties=global_prop['short_MD_mdrun'])

        # Use MDAnalysis to analyze the CV (distance between 2 user-defined groups)
        last_step_accepted = analyze_cv(global_paths['short_MD_mdrun']["output_xtc_path"], global_paths['short_MD_mdrun']["output_gro_path"], 
                                        continuation_MD_prop['mdp'], accepted_steps, sumd_prop, sumd_log, output_path = conf.get_working_dir_path())
        
        if last_step_accepted:

            # First accepted step
            if not os.path.exists(global_paths["original_trajectory_cat"]["output_trj_path"]):
                
                # Copy short MD trajectory to concatenation folders
                shutil.copyfile(global_paths['short_MD_mdrun']['output_xtc_path'], global_paths["original_trajectory_cat"]["output_trj_path"])
            
            # Subsequent accepted steps
            else:

                # remove previous zip trajectory bundles if they exists
                remove_tmp_files(global_paths["original_trajectory_cat"]["input_trj_zip_path"])

                # Create a ZipFile objects
                zipObject = ZipFile(global_paths["original_trajectory_cat"]["input_trj_zip_path"], 'w')

                # Add previously concatenated trajectories to zip files
                zipObject.write(global_paths["original_trajectory_cat"]["output_trj_path"])

                # Add short MD trajectory to zip files
                zipObject.write(global_paths['short_MD_mdrun']['output_xtc_path'])

                # Close zip files
                zipObject.close()

                # Concatenate
                trjcat(**global_paths["original_trajectory_cat"], properties=global_prop["original_trajectory_cat"])

            # Reset failed steps counter
            failed_steps = 0

            # Increase accepted steps counter
            accepted_steps += 1

        else:

            # Increase failed steps counter
            failed_steps += 1

        if last_step_accepted or (failed_steps > sumd_prop['num_steps']['max_failed']):
            
            # Replace last accepted structure and checkpoint for newly accepted ones
            shutil.copyfile(global_paths["short_MD_mdrun"]["output_gro_path"], continuation_MD_paths["input_gro_path"])
            shutil.copyfile(global_paths["short_MD_mdrun"]["output_cpt_path"], continuation_MD_paths["input_cpt_path"])

        # Increase total counter
        total_steps += 1

        # Remove log files from previous steps
        remove_logs(continuation_MD_prop, global_prop)

    # Get all temporal unique folders with -> will be fixed
    tmp_folders=get_tmp_folders(os.getcwd())

    # Move all temporal unique folders to common folder -> will be fixed
    for tmp_folder in tmp_folders:
        shutil.move(tmp_folder, os.path.join(os.getcwd(), "tmp"))

    # Try to image the concatenated trajectory

    # Make molecule whole
    gmx_image(**global_paths['pbc_1_whole'], properties=global_prop['pbc_1_whole'])

    # Cluster molecules - e.g. RNA strands and ligand
    gmx_image(**global_paths['pbc_2_cluster'], properties=global_prop['pbc_2_cluster'])

    # Extract initial frame from trajectory
    gmx_trjconv_trj(**global_paths['pbc_3_extract_frame'], properties=global_prop['pbc_3_extract_frame'])

    # Avoid jumps of the ligand and use extracted frame as reference
    gmx_image(**global_paths['pbc_4_nojump'], properties=global_prop['pbc_4_nojump'])
    
    # Center both the RNA and the ligand in the box
    gmx_image(**global_paths['pbc_5_center'], properties=global_prop['pbc_5_center'])

    # Minimize RMSD of RNA atoms wrt reference (fitting)
    gmx_image(**global_paths['pbc_6_fit'], properties=global_prop['pbc_6_fit'])

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
