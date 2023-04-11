#!/usr/bin/env python3

# Importing all the needed libraries
import time
import argparse
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_pydock.pydock.setup import setup
from biobb_pydock.pydock.ftdock import ftdock
from biobb_pydock.pydock.dockser import dockser

def main_wf(configuration_path):
    '''
    Main ... workflow 

    Inputs
    ------

        configuration_path (str): path to input.yml

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    start_time = time.time()
    
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Sampling protein-protein docking poses

    # STEP 1
    global_log.info("step1_setup: setup receptor and ligand proteins for pyDock")
    setup(**global_paths["step1_setup"], properties=global_prop["step1_setup"])
    
    # STEP 2
    global_log.info("step2_ftdock: sample docking poses using ftdock (FFT-based algorithm)")
    ftdock(**global_paths["step2_ftdock"], properties=global_prop["step2_ftdock"])

    # Scoring protein-protein docking poses

    # STEP 3
    global_log.info("step3_dockser: score docking poses using pyDock")
    dockser(**global_paths["step3_dockser"], properties=global_prop["step3_dockser"])
    
    # Clustering with RMSD

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

    return global_paths, global_prop


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Protein-protein docking workflow")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    args = parser.parse_args()

    main_wf(configuration_path=args.config_path)
