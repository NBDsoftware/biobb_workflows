#!/usr/bin/env python3

import time
import argparse
import subprocess
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_chemistry.ambertools.reduce_remove_hydrogens import reduce_remove_hydrogens
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.cat_pdb import cat_pdb
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
from biobb_amber.leap.leap_gen_top import leap_gen_top
from biobb_amber.sander.sander_mdrun import sander_mdrun
from biobb_amber.process.process_minout import process_minout
from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb
from biobb_amber.leap.leap_solvate import leap_solvate
from biobb_amber.leap.leap_add_ions import leap_add_ions
from biobb_amber.process.process_mdout import process_mdout
from biobb_amber.pmemd.pmemd_mdrun import pmemd_mdrun
from biobb_analysis.ambertools.cpptraj_rms import cpptraj_rms
from biobb_analysis.ambertools.cpptraj_rgyr import cpptraj_rgyr
from biobb_analysis.ambertools.cpptraj_image import cpptraj_image


def main(configuration_path, output_path = None, input_pdb_path = None):
    start_time = time.time()
    conf = settings.ConfReader(configuration_path)
    # Enforce output_path if provided
    if output_path is not None:
        output_path = fu.get_working_dir_path(output_path, restart = conf.properties.get('restart', 'False'))
        conf.properties['working_dir_path'] = output_path
    else:
        output_path = conf.get_working_dir_path()

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    
    global_paths["step00_reduce_remove_hydrogens"]["input_path"]=input_pdb_path
    global_log.info("step00_reduce_remove_hydrogens: Removing Hydrogens")
    reduce_remove_hydrogens(**global_paths["step00_reduce_remove_hydrogens"], properties=global_prop["step00_reduce_remove_hydrogens"])

    global_log.info("step0_extract_molecule: Extracting Protein")
    extract_molecule(**global_paths["step0_extract_molecule"], properties=global_prop["step0_extract_molecule"])

    global_log.info("step000_cat_pdb: Concatenating protein with included ions")
    cat_pdb(**global_paths["step000_cat_pdb"], properties=global_prop["step000_cat_pdb"])

    global_log.info("step1_pdb4amber_run: Preparing PDB file for AMBER")
    pdb4amber_run(**global_paths["step1_pdb4amber_run"], properties=global_prop["step1_pdb4amber_run"])

    global_log.info("step2_leap_gen_top: Create protein system topology")
    leap_gen_top(**global_paths["step2_leap_gen_top"], properties=global_prop["step2_leap_gen_top"])

    global_log.info("step3_sander_mdrun_minH: Minimize Hydrogens")
    sander_mdrun(**global_paths["step3_sander_mdrun_minH"], properties=global_prop["step3_sander_mdrun_minH"])

    global_log.info("step4_process_minout_minH: Checking Energy Minimization results")
    process_minout(**global_paths["step4_process_minout_minH"], properties=global_prop["step4_process_minout_minH"])

    global_log.info("step5_sander_mdrun_min: Minimize the system")
    sander_mdrun(**global_paths["step5_sander_mdrun_min"], properties=global_prop["step5_sander_mdrun_min"])

    global_log.info("step6_process_minout_min: Checking Energy Minimization results")
    process_minout(**global_paths["step6_process_minout_min"], properties=global_prop["step6_process_minout_min"])

    global_log.info("step7_amber_to_pdb: Getting minimized structure")
    amber_to_pdb(**global_paths["step7_amber_to_pdb"], properties=global_prop["step7_amber_to_pdb"])

    global_log.info("step8_leap_solvate: Create water box")
    leap_solvate(**global_paths["step8_leap_solvate"], properties=global_prop["step8_leap_solvate"])

    global_log.info("step9_leap_add_ions: Adding ions")
    leap_add_ions(**global_paths["step9_leap_add_ions"], properties=global_prop["step9_leap_add_ions"])

    global_log.info("step10_sander_mdrun_energy: Running Energy Minimization")
    sander_mdrun(**global_paths["step10_sander_mdrun_energy"], properties=global_prop["step10_sander_mdrun_energy"])

    global_log.info("step11_process_minout_energy: Checking Energy Minimization results")
    process_minout(**global_paths["step11_process_minout_energy"], properties=global_prop["step11_process_minout_energy"])

    global_log.info("step12_pmemd_mdrun_warm: Warming up the system")
    pmemd_mdrun(**global_paths["step12_pmemd_mdrun_warm"], properties=global_prop["step12_pmemd_mdrun_warm"])

    global_log.info("step13_process_mdout_warm: Checking results from the system warming up")
    process_mdout(**global_paths["step13_process_mdout_warm"], properties=global_prop["step13_process_mdout_warm"])

    global_log.info("step14_pmemd_mdrun_nvt: Equilibrating the system (NVT)")
    pmemd_mdrun(**global_paths["step14_pmemd_mdrun_nvt"], properties=global_prop["step14_pmemd_mdrun_nvt"])

    global_log.info("step15_process_mdout_nvt: Checking NVT Equilibration results")
    process_mdout(**global_paths["step15_process_mdout_nvt"], properties=global_prop["step15_process_mdout_nvt"])

    global_log.info("step16_pmemd_mdrun_npt: Equilibrating the system (NPT)")
    pmemd_mdrun(**global_paths["step16_pmemd_mdrun_npt"], properties=global_prop["step16_pmemd_mdrun_npt"])

    global_log.info("step17_process_mdout_npt: Checking NPT Equilibration results")
    process_mdout(**global_paths["step17_process_mdout_npt"], properties=global_prop["step17_process_mdout_npt"])


    global_log.info("step18_pmemd_free_mdrun: Creating portable binary run file to run a free MD simulation")
    pmemd_mdrun(**global_paths["step18_pmemd_free_mdrun"], properties=global_prop["step18_pmemd_free_mdrun"])



    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    if system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="ABC MD Setup pipeline using BioExcel Building Blocks")
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)
    parser.add_argument('--system', required=False)
    parser.add_argument('--input_pdb', dest='input_pdb_path',
                        help="Input PDB file (default: input_structure_path in step 1 of configuration file)",
                        required=False)
    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    args = parser.parse_args()
    main(configuration_path=args.config_path, input_pdb_path=args.input_pdb_path, output_path=args.output_path)
