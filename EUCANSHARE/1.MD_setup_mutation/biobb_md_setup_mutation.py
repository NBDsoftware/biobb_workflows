#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import time
import argparse
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
# from biobb_io.api.pdb import pdb
# from biobb_io.api.canonical_fasta import canonical_fasta
# from biobb_model.model.fix_backbone import fix_backbone
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_altlocs import fix_altlocs
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.fix_chirality import fix_chirality
from biobb_model.model.mutate import mutate
from biobb_gromacs.gromacs.trjcat import trjcat
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.genion import genion
from biobb_gromacs.gromacs.mdrun import mdrun
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure


def main_wf(configuration_path, last_step, mutation_list, num_trajs):
    '''
    Main MD Setup, mutation and run workflow. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system to simulate, minimize it, equilibrate it and finally do N production runs.

    Inputs
    ------

        configuration_path (str): path to input.yml
        input_pdb          (str): either the PDB code in pdb:id format or the path the the pdb file
        last_step          (str): last step of the workflow to execute ('pdb', 'fix', 'min', 'nvt', 'npt' or 'all')
        mutation_list      (str): Mutations as comma-separated list with the format: 'chain_id : Old_residue_code Residue_number New_residue_code'.
                                Examples: 'A:G34T' or 'A:F38C,A:N39W,A:T40G'
        num_trajs          (int): number of trajectories to launch one after another

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    start_time = time.time()

    # Set default value for 'last_step' arg
    if last_step is None:
        last_step = 'all'

    # Set default value for 'n_trajs' arg
    if num_trajs is None:
        num_trajs = 1
    else:
        num_trajs = int(num_trajs)

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # STEP 1: extracts molecule of interest: protein
    global_log.info("step1_extractMolecule: extract molecule of interest (protein)")
    extract_molecule(**global_paths["step1_extractMolecule"], properties=global_prop["step1_extractMolecule"])

    # STEP 2 (A): Fix alternative locations
    global_log.info("step2A_fixaltlocs: Fix alternative locations")
    fix_altlocs(**global_paths["step2A_fixaltlocs"], properties=global_prop["step2A_fixaltlocs"])

    # STEP 2 (B): Add mutations if requested
    global_log.info("step2B_mutations: Preparing mutated structure")
    mutate(**global_paths["step2B_mutations"], properties=global_prop["step2B_mutations"])

    # STEP 2 (C): Get canonical FASTA
    # global_log.info("step2C_canonical_fasta: Get canonical FASTA")
    # canonical_fasta(**global_paths["step2C_canonical_fasta"], properties=global_prop["step2C_canonical_fasta"])

    # STEP 2 (D): Model missing heavy atoms of backbone
    # global_log.info("step2D_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
    # fix_backbone(**global_paths["step2D_fixbackbone"], properties=global_prop["step2D_fixbackbone"])

    # STEP 2 (E): model missing heavy atoms of side chains
    global_log.info("step2E_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
    fix_side_chain(**global_paths["step2E_fixsidechain"], properties=global_prop["step2E_fixsidechain"])

    # STEP 2 (F): model SS bonds (CYS -> CYX) if necessary
    global_log.info("step2F_fixssbonds: Fix SS bonds")
    fix_ssbonds(**global_paths["step2F_fixssbonds"], properties=global_prop["step2F_fixssbonds"])

    # STEP 2 (G): Fix amides
    global_log.info("step2G_fixamides: fix clashing amides")
    fix_amides(**global_paths["step2G_fixamides"], properties=global_prop["step2G_fixamides"])

    # STEP 2 (H): Fix chirality
    global_log.info("step2H_fixchirality: fix chirality of residues")
    fix_chirality(**global_paths["step2H_fixchirality"], properties=global_prop["step2H_fixchirality"])

    # STEP 2 (I): renumber structure atoms and residues
    global_log.info("step2I_renumberstructure: renumber structure")
    renumber_structure(**global_paths["step2I_renumberstructure"], properties=global_prop["step2I_renumberstructure"])

    # STEP 3: add H atoms, generate coordinate (.gro) and topology (.top) file
    global_log.info("step3_pdb2gmx: Generate the topology")
    pdb2gmx(**global_paths["step3_pdb2gmx"], properties=global_prop["step3_pdb2gmx"])

    # STEP 4: Create simulation box
    global_log.info("step4_editconf: Create the solvent box")
    editconf(**global_paths["step4_editconf"], properties=global_prop["step4_editconf"])

    # STEP 5: Add solvent molecules
    global_log.info("step5_solvate: Fill the solvent box with water molecules")
    solvate(**global_paths["step5_solvate"], properties=global_prop["step5_solvate"])

    # STEP 6: ion generation pre-processing
    global_log.info("step6_grompp_genion: Preprocess ion generation")
    grompp(**global_paths["step6_grompp_genion"], properties=global_prop["step6_grompp_genion"])

    # STEP 7: ion generation
    global_log.info("step7_genion: Ion generation")
    genion(**global_paths["step7_genion"], properties=global_prop["step7_genion"])

    # STEP 8: minimization pre-processing
    global_log.info("step8_grompp_min: Preprocess energy minimization")
    grompp(**global_paths["step8_grompp_min"], properties=global_prop["step8_grompp_min"])

    # STEP 9: minimization
    global_log.info("step9_mdrun_min: Execute energy minimization")
    mdrun(**global_paths["step9_mdrun_min"], properties=global_prop["step9_mdrun_min"])

    # STEP 10: dump potential energy evolution NOTE: create plots automatically
    global_log.info("step10_energy_min: Compute potential energy during minimization")
    gmx_energy(**global_paths["step10_energy_min"], properties=global_prop["step10_energy_min"])

    # STEP 11: NVT equilibration pre-processing
    global_log.info("step11_grompp_nvt: Preprocess system temperature equilibration")
    grompp(**global_paths["step11_grompp_nvt"], properties=global_prop["step11_grompp_nvt"])

    # STEP 12: NVT equilibration
    global_log.info("step12_mdrun_nvt: Execute system temperature equilibration")
    mdrun(**global_paths["step12_mdrun_nvt"], properties=global_prop["step12_mdrun_nvt"])

    # STEP 13: dump temperature evolution
    global_log.info("step13_energy_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step13_energy_nvt"], properties=global_prop["step13_energy_nvt"])

    # STEP 14: NPT equilibration pre-processing
    global_log.info("step14_grompp_npt: Preprocess system pressure equilibration")
    grompp(**global_paths["step14_grompp_npt"], properties=global_prop["step14_grompp_npt"])

    # STEP 15: NPT equilibration
    global_log.info("step15_mdrun_npt: Execute system pressure equilibration")
    mdrun(**global_paths["step15_mdrun_npt"], properties=global_prop["step15_mdrun_npt"])

    # STEP 16: dump density and pressure evolution
    global_log.info("step16_energy_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step16_energy_npt"], properties=global_prop["step16_energy_npt"])

    # Loop over trajectories
    global_log.info(f"Number of trajectories: {num_trajs}")
    trj_list = []
    for traj in (f"traj_{i}" for i in range(num_trajs)):
        traj_prop = conf.get_prop_dic(prefix=traj)
        traj_paths = conf.get_paths_dic(prefix=traj)

        # STEP 17: free NPT production run pre-processing
        global_log.info(f"{traj} >  step17_grompp_md: Preprocess free dynamics")
        traj_paths['step17_grompp_md']['input_gro_path'] = global_paths["step17_grompp_md"]['input_gro_path']
        traj_paths['step17_grompp_md']['input_top_zip_path'] = global_paths["step17_grompp_md"]['input_top_zip_path']
        traj_paths['step17_grompp_md']['input_cpt_path'] = global_paths["step17_grompp_md"]['input_cpt_path']
        grompp(**traj_paths['step17_grompp_md'], properties=traj_prop["step17_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{traj} >  step18_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step18_mdrun_md'], properties=traj_prop['step18_mdrun_md'])

        # STEP 19: dump RMSD with respect to equilibrated structure (first frame) NOTE: add computation of RMSF, PCA and projection onto PC for visualization
        global_log.info(f"{traj} >  step19_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")
        gmx_rms(**traj_paths['step19_rmsfirst'], properties=traj_prop['step19_rmsfirst'])

        # STEP 20: dump RMSD with respect to minimized structure
        global_log.info(f"{traj} >  step20_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
        traj_paths['step20_rmsexp']['input_structure_path'] = global_paths["step20_rmsexp"]['input_structure_path']
        gmx_rms(**traj_paths['step20_rmsexp'], properties=traj_prop['step20_rmsexp'])

        # STEP 21: dump Radius of gyration NOTE: add computation of RMSF, PCA and projection onto PC for visualization
        global_log.info(f"{traj} >  step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
        gmx_rgyr(**traj_paths['step21_rgyr'], properties=traj_prop['step21_rgyr'])

        # STEP 22: image (correct for PBC) the trajectory centering the protein and dumping only protein atoms
        global_log.info(f"{traj} >  step22_image: Imaging the resulting trajectory")
        gmx_image(**traj_paths['step22_image'], properties=traj_prop['step22_image'])
        trj_list.append(traj_paths['step22_image']['output_traj_path'])

        # STEP 23: remove water and ions from structure obtained after equilibration, before production run
        global_log.info(f"{traj} >  step23_dry: Removing water molecules and ions from the equilibrated structure")
        traj_paths['step23_dry']['input_structure_path'] = global_paths["step23_dry"]['input_structure_path']
        gmx_trjconv_str(**traj_paths['step23_dry'], properties=traj_prop['step23_dry'])

    # STEP 24: concatenate trajectories
    global_log.info("step24_trjcat: Concatenate trajectories")
    fu.zip_list(zip_file=global_paths["step24_trjcat"]['input_trj_zip_path'], file_list=trj_list)
    trjcat(**global_paths["step24_trjcat"], properties=global_prop["step24_trjcat"])

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

    parser = argparse.ArgumentParser("Simple MD Protein Setup")

    parser.add_argument('--config', dest='config_path',
                        help="Configuration file (YAML)",
                        required=True)

    # Execute workflow until 'last_step' -> all executes all steps (all is default)
    parser.add_argument('--until', dest='last_step',
                        help="Extent of the pipeline to execute (pdb, fix, min, nvt, npt or all; default: all)",
                        required=False)

    parser.add_argument('--mut_list', dest='mut_list',
                        help="Mutations as comma-separated list with the format: 'chain_id : Old_residue_code Residue_number New_residue_code'. Examples: 'A:G34T' or 'A:F38C,A:N39W,A:T40G' (default: None)",
                        required=False)

    parser.add_argument('--num_trajs', dest='num_trajs',
                        help="Number of trajectories (default: 1)",
                        required=False)

    args = parser.parse_args()

    _, _ = main_wf(configuration_path=args.config_path,
                   last_step=args.last_step,
                   mutation_list=args.mut_list,
                   num_trajs=args.num_trajs)
