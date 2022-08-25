#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import sys
import os
import time
import argparse
import shutil
from pathlib import Path, PurePath

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_io.api.pdb import Pdb
from biobb_model.model.fix_side_chain import FixSideChain
from biobb_model.model.mutate import Mutate
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
from biobb_md.gromacs.editconf import Editconf
from biobb_md.gromacs.solvate import Solvate
from biobb_md.gromacs.grompp import Grompp
from biobb_md.gromacs.genion import Genion
from biobb_md.gromacs.mdrun import Mdrun
from biobb_analysis.gromacs.gmx_rms import GMXRms
from biobb_analysis.gromacs.gmx_rgyr import GMXRgyr
from biobb_analysis.gromacs.gmx_energy import GMXEnergy
from biobb_analysis.gromacs.gmx_image import GMXImage
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr

def prep_output_file(input_file_path, output_file_path):
    wdir = PurePath(output_file_path).parent
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    shutil.copy(input_file_path, output_file_path)

def write_pdb_from_gro(output_pdb_path, input_gro_path):
    prop = {'selection': 'Protein'}
    GMXTrjConvStr(
        input_structure_path=input_gro_path,
        input_top_path=input_gro_path,
        output_str_path=output_pdb_path,
        properties=prop
    ).launch()
    
def run_wf(args):
    if 'to_do' not in args:
        args.to_do = 'free'
        
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(args.config_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Declaring the steps of the workflow, one by one 
    # Using as inputs the global paths and global properties
    # identified by the corresponding step name
    # Writing information about each step to the global log 

    if 'pdb:' in args.input_pdb_path:
        # Downloading desired PDB file
        global_log.info("step1_pdb: Downloading from PDB")
        
        pdbCode = args.input_pdb_path.split(':')[1]
        
        global_log.info("step1_pdb: Downloading {} from PDB".format(pdbCode))
        prop = {
            'path': PurePath(global_paths["step1_pdb"]["output_pdb_path"]).parent,
            'pdb_code': pdbCode
        }
        Pdb(**global_paths["step1_pdb"], properties=prop).launch()
    else:
        global_log.info("step1_pdb: Adding input PDB ({}) to working dir".format(args.input_pdb_path))
        prep_output_file(args.input_pdb_path, global_paths["step1_pdb"]["output_pdb_path"])
        pdbCode = os.path.splitext(args.input_pdb_path)[0]

    if args.mut_list:
        global_log.info("step1_mutations: Preparing mutated structure")
        prop = {
            'path': PurePath(global_paths["step1_mutations"]["output_pdb_path"]).parent,
            'mutation_list': args.mut_list
        }
        Mutate(**global_paths["step1_mutations"], properties=prop).launch()
    else:
        prep_output_file(global_paths["step1_pdb"]["output_pdb_path"], global_paths["step1_mutations"]["output_pdb_path"])
    
    global_log.info("step2_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
    FixSideChain(**global_paths["step2_fixsidechain"], properties=global_prop["step2_fixsidechain"]).launch()

    if args.to_do == 'fix':
        shutil.copy(global_paths["step2_fixsidechain"]["output_pdb_path"], args.output_pdb_path)
        global_log.info("Fix completed. Final structure saved on " + args.output_pdb_path)
        return 0
        
    global_log.info("step3_pdb2gmx: Generate the topology")
    Pdb2gmx(**global_paths["step3_pdb2gmx"], properties=global_prop["step3_pdb2gmx"]).launch()

    global_log.info("step4_editconf: Create the solvent box")
    Editconf(**global_paths["step4_editconf"], properties=global_prop["step4_editconf"]).launch()

    global_log.info("step5_solvate: Fill the solvent box with water molecules")
    Solvate(**global_paths["step5_solvate"], properties=global_prop["step5_solvate"]).launch()

    global_log.info("step6_grompp_genion: Preprocess ion generation")
    Grompp(**global_paths["step6_grompp_genion"], properties=global_prop["step6_grompp_genion"]).launch()

    global_log.info("step7_genion: Ion generation")
    Genion(**global_paths["step7_genion"], properties=global_prop["step7_genion"]).launch()

    global_log.info("step8_grompp_min: Preprocess energy minimization")
    Grompp(**global_paths["step8_grompp_min"], properties=global_prop["step8_grompp_min"]).launch()

    global_log.info("step9_mdrun_min: Execute energy minimization")
    Mdrun(**global_paths["step9_mdrun_min"], properties=global_prop["step9_mdrun_min"]).launch()
    
    global_log.info("step10_energy_min: Compute potential energy during minimization")
    GMXEnergy(**global_paths["step10_energy_min"], properties=global_prop["step10_energy_min"]).launch()

    if args.to_do == 'min':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step9_mdrun_min"]["output_gro_path"])
        global_log.info("Minimization completed. Final structure saved on " + args.output_pdb_path)
        return 0
       
    global_log.info("step11_grompp_nvt: Preprocess system temperature equilibration")
    Grompp(**global_paths["step11_grompp_nvt"], properties=global_prop["step11_grompp_nvt"]).launch()

    global_log.info("step12_mdrun_nvt: Execute system temperature equilibration")
    Mdrun(**global_paths["step12_mdrun_nvt"], properties=global_prop["step12_mdrun_nvt"]).launch()

    global_log.info("step13_energy_nvt: Compute temperature during NVT equilibration")
    GMXEnergy(**global_paths["step13_energy_nvt"], properties=global_prop["step13_energy_nvt"]).launch()

    if args.to_do == 'nvt':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step12_mdrun_nvt"]["output_gro_path"])
        global_log.info("NVT Equilibration completed. Final structure saved on " + args.output_pdb_path)
        return 0
    
    global_log.info("step14_grompp_npt: Preprocess system pressure equilibration")
    Grompp(**global_paths["step14_grompp_npt"], properties=global_prop["step14_grompp_npt"]).launch()

    global_log.info("step15_mdrun_npt: Execute system pressure equilibration")
    Mdrun(**global_paths["step15_mdrun_npt"], properties=global_prop["step15_mdrun_npt"]).launch()

    global_log.info("step16_energy_npt: Compute Density & Pressure during NPT equilibration")
    GMXEnergy(**global_paths["step16_energy_npt"], properties=global_prop["step16_energy_npt"]).launch()

    if args.to_do == 'npt':
        write_pdb_from_gro(args.output_pdb_path, global_paths["step15_mdrun_npt"]["output_gro_path"])
        global_log.info("NPT Equilibration completed. Final structure saved on " + args.output_pdb_path)
        return 0
    
    global_log.info("step17_grompp_md: Preprocess free dynamics")
    Grompp(**global_paths["step17_grompp_md"], properties=global_prop["step17_grompp_md"]).launch()

    global_log.info("step18_mdrun_md: Execute free molecular dynamics simulation")
    Mdrun(**global_paths["step18_mdrun_md"], properties=global_prop["step18_mdrun_md"]).launch()

    global_log.info("step19_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")
    GMXRms(**global_paths["step19_rmsfirst"], properties=global_prop["step19_rmsfirst"]).launch()

    global_log.info("step20_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
    GMXRms(**global_paths["step20_rmsexp"], properties=global_prop["step20_rmsexp"]).launch()

    global_log.info("step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
    GMXRgyr(**global_paths["step21_rgyr"], properties=global_prop["step21_rgyr"]).launch()

    global_log.info("step22_image: Imaging the resulting trajectory")
    GMXImage(**global_paths["step22_image"], properties=global_prop["step22_image"]).launch()

    global_log.info("step23_dry: Removing water molecules and ions from the resulting structure")
    GMXTrjConvStr(**global_paths["step23_dry"], properties=global_prop["step23_dry"]).launch()

    write_pdb_from_gro(args.output_pdb_path, global_paths["step18_mdrun_md"]["output_gro_path"])
    global_log.info("Free Equilibration completed. Final structure saved on " + args.output_pdb_path)
    return 0
    
def main():
    parser = argparse.ArgumentParser("Simple MD Protein Setup")
    parser.add_argument('-i', dest='input_pdb_path',
                        help="Input pdb file or id (as pdb:id)", required=True)
    parser.add_argument('-o', dest='output_pdb_path',
                        help="Output pdb file", required=True)
    parser.add_argument(
        '--op', dest='to_do', help="Extent of the pipeline to execute (fix, min, nvt, npt, free")
    parser.add_argument('--mut_list', dest='mut_list',
                        help='Mutations list as *|A:V45W[,...]')
    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)")
    
    args = parser.parse_args()

    run_wf(args)


if __name__ == "__main__":
    main()
