##!/usr/bin/env python3

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing python modules
import os
import re
import glob
import time
import argparse
import shutil
from pathlib import Path

# Import pre-defined workflows

# Import biobb
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

def main_wf():
    '''
   
    Inputs
    ------


    Outputs
    -------


    '''



    start_time = time.time()

    # Receiving the main input configuration file (YAML)
    conf = settings.ConfReader(...)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Launch several short MD trajectories

    # Use PLUMED to analyze a certain colvar - take plumed file as input (the user is responsible to make it work)

    # Find those trajectories that are giving a negative slope in the colvar

    # Restart those from last checkpoint and start the rest from the beginning - several possibilities arise here, 
    # we can focus on the best trajectory or try to explore as much as possible or avoid getting trapped


    # Print timing information to the log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % ... )
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Simple clustering, cavity analysis and docking pipeline using BioExcel Building Blocks")

    parser.add_argument('--main-config', dest='main_config_path',
                        help="Configuration file for pocket VS (YAML)", 
                        required=True)


    args = parser.parse_args()

    main_wf(  )

