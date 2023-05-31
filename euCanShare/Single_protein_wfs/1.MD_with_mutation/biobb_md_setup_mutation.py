#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import time
import random
import argparse

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_io.api.canonical_fasta import canonical_fasta
from biobb_model.model.fix_backbone import fix_backbone
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
from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_analysis.gromacs.gmx_rms import gmx_rms
from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure


def main_wf(configuration_path, setup_only, num_trajs, output_path = None, input_pdb_path = None, pdb_chains = None,
            mutation_list = None, input_gro_path = None, input_top_path = None, fix_backbn = None, fix_ss = None):
    '''
    Main MD Setup, mutation and run workflow. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system, minimize it, equilibrate it and finally do N production runs.

    Inputs
    ------

        configuration_path (str): path to YAML configuration file
        setup_only        (bool): whether to only setup the system or also run the simulations
        num_trajs          (int): number of trajectories to launch in serial
        output_path        (str): (Optional) path to output folder
        input_pdb_path     (str): (Optional) path to input PDB file
        pdb_chains         (str): (Optional) chains to be extracted from the input PDB file
        mutation_list      (str): (Optional) list of mutations to be introduced in the structure
        input_gro_path     (str): (Optional) path to input structure file (.gro)
        input_top_path     (str): (Optional) path to input topology file (.zip)
        fix_backbn        (bool): (Optional) wether to add missing backbone atoms
        fix_ss            (bool): (Optional) wether to add disulfide bonds

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    start_time = time.time()

    # Receiving the input configuration file (YAML)
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

    # If prepared structure is not provided
    if input_gro_path is None:

        # If input PDB is given as argument
        if input_pdb_path is not None:
            global_paths["step1_extractMolecule"]["input_structure_path"] = input_pdb_path

        # If chains are given as argument
        if pdb_chains is not None:
            global_prop["step1_extractMolecule"]["chains"] = pdb_chains

        # STEP 1: extracts molecule of interest
        global_log.info("step1_extractMolecule: extract molecule of interest (protein)")
        extract_molecule(**global_paths["step1_extractMolecule"], properties=global_prop["step1_extractMolecule"])

        # STEP 2 (A): Fix alternative locations
        global_log.info("step2A_fixaltlocs: Fix alternative locations")
        fix_altlocs(**global_paths["step2A_fixaltlocs"], properties=global_prop["step2A_fixaltlocs"])

        if mutation_list is not None:
            global_prop["step2B_mutations"]["mutation_list"] = ",".join(mutation_list)

        # STEP 2 (B): Add mutations if requested
        global_log.info("step2B_mutations: Preparing mutated structure")
        mutate(**global_paths["step2B_mutations"], properties=global_prop["step2B_mutations"])

        if fix_backbn:
            # STEP 2 (C): Get canonical FASTA
            global_log.info("step2C_canonical_fasta: Get canonical FASTA")
            canonical_fasta(**global_paths["step2C_canonical_fasta"], properties=global_prop["step2C_canonical_fasta"])

            # STEP 2 (D): Model missing heavy atoms of backbone
            global_log.info("step2D_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
            fix_backbone(**global_paths["step2D_fixbackbone"], properties=global_prop["step2D_fixbackbone"])
        else:
            global_paths['step2E_fixsidechain']['input_pdb_path'] = global_paths['step2B_mutations']['output_pdb_path']

        # STEP 2 (E): model missing heavy atoms of side chains
        global_log.info("step2E_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        fix_side_chain(**global_paths["step2E_fixsidechain"], properties=global_prop["step2E_fixsidechain"])

        if fix_ss:
            # STEP 2 (F): model SS bonds (CYS -> CYX)
            global_log.info("step2F_fixssbonds: Fix SS bonds")
            fix_ssbonds(**global_paths["step2F_fixssbonds"], properties=global_prop["step2F_fixssbonds"])
        else:
            global_paths['step2G_fixamides']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

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

        if setup_only:
            global_log.info("setup_only: setup_only flag is set to True, exiting...")
            return

    else:

        global_paths['step8_grompp_min']['input_gro_path'] = input_gro_path
        global_paths['step8_grompp_min']['input_top_zip_path'] = input_top_path

    # STEP 8: minimization pre-processing
    global_log.info("step8_grompp_min: Preprocess energy minimization")
    grompp(**global_paths["step8_grompp_min"], properties=global_prop["step8_grompp_min"])

    # STEP 9: minimization
    global_log.info("step9_mdrun_min: Execute energy minimization")
    mdrun(**global_paths["step9_mdrun_min"], properties=global_prop["step9_mdrun_min"])

    # STEP 9B: create index file
    global_log.info("step9B_make_ndx: Create index file")
    make_ndx(**global_paths["step9B_make_ndx"], properties=global_prop["step9B_make_ndx"])

    # STEP 10: dump potential energy evolution
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

    # Production trajectories
    global_log.info(f"Number of trajectories: {num_trajs}")
    traj_list = []
    for traj in (f"traj_{i}" for i in range(num_trajs)):

        traj_prop = conf.get_prop_dic(prefix=traj)
        traj_paths = conf.get_paths_dic(prefix=traj)

        # Change seed of temperature coupling algorithm (V-rescale)
        new_seed = {'ld-seed': random.randint(1, 1000000)}
        traj_prop['step17_grompp_md']['mdp'].update(new_seed)

        # STEP 17: free NPT production run pre-processing
        global_log.info(f"{traj} >  step17_grompp_md: Preprocess free dynamics")
        # Update common paths to all trajectories
        traj_paths['step17_grompp_md']['input_gro_path'] = global_paths["step17_grompp_md"]['input_gro_path']
        traj_paths['step17_grompp_md']['input_top_zip_path'] = global_paths["step17_grompp_md"]['input_top_zip_path']
        if global_paths["step17_grompp_md"].get('input_ndx_path'):
            traj_paths['step17_grompp_md']['input_ndx_path'] = global_paths["step17_grompp_md"]['input_ndx_path']
        if global_paths["step17_grompp_md"].get('input_cpt_path'):
            traj_paths['step17_grompp_md']['input_cpt_path'] = global_paths["step17_grompp_md"]['input_cpt_path']
        grompp(**traj_paths['step17_grompp_md'], properties=traj_prop["step17_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{traj} >  step18_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step18_mdrun_md'], properties=traj_prop['step18_mdrun_md'])

        # STEP 19: dump RMSD with respect to equilibrated structure (first frame)
        global_log.info(f"{traj} >  step19_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")
        gmx_rms(**traj_paths['step19_rmsfirst'], properties=traj_prop['step19_rmsfirst'])

        # STEP 20: dump RMSD with respect to minimized structure
        global_log.info(f"{traj} >  step20_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
        # Update common paths to all trajectories
        traj_paths['step20_rmsexp']['input_structure_path'] = global_paths["step20_rmsexp"]['input_structure_path']
        gmx_rms(**traj_paths['step20_rmsexp'], properties=traj_prop['step20_rmsexp'])

        # STEP 21: dump Radius of gyration
        global_log.info(f"{traj} >  step21_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
        gmx_rgyr(**traj_paths['step21_rgyr'], properties=traj_prop['step21_rgyr'])

        # STEP 22: dump RMSF
        global_log.info(f"{traj} >  step22_rmsf: Compute Root Mean Square Fluctuation to measure the protein flexibility during the free MD simulation")
        # Update common paths to all trajectories
        traj_paths['step22_rmsf']['input_top_path'] = global_paths["step22_rmsf"]['input_top_path']
        cpptraj_rmsf(**traj_paths['step22_rmsf'], properties=traj_prop['step22_rmsf'])

        # STEP 23: image the trajectory
        traj_paths['step23_image_traj']['input_index_path'] = global_paths["step23_image_traj"]['input_index_path']
        global_log.info(f"{traj} >  step23_image_traj: Imaging the trajectory")
        gmx_image(**traj_paths['step23_image_traj'], properties=traj_prop['step23_image_traj'])

        # STEP 24: fit the trajectory
        traj_paths['step24_fit_traj']['input_index_path'] = global_paths["step24_fit_traj"]['input_index_path']
        global_log.info(f"{traj} >  step24_fit_traj: Fit the trajectory")
        gmx_image(**traj_paths['step24_fit_traj'], properties=traj_prop['step24_fit_traj'])

        traj_list.append(traj_paths['step24_fit_traj']['output_traj_path'])

    # STEP 25: obtain dry structure
    global_log.info("step25_dry_str: Obtain dry structure")
    gmx_trjconv_str(**global_paths["step25_dry_str"], properties=global_prop["step25_dry_str"])

    # STEP 26: concatenate trajectories
    global_log.info("step26_trjcat: Concatenate trajectories")
    fu.zip_list(zip_file=global_paths["step26_trjcat"]['input_trj_zip_path'], file_list=traj_list)
    trjcat(**global_paths["step26_trjcat"], properties=global_prop["step26_trjcat"])

    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % output_path)
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

    parser.add_argument('--setup_only', action='store_true',
                        help="Only setup the system (default: False)",
                        required=False)

    parser.add_argument('--num_trajs', dest='num_trajs',
                        help="Number of trajectories (default: 1)",
                        required=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)

    parser.add_argument('--input_pdb', dest='input_pdb_path',
                        help="Input PDB file (default: input_structure_path in step 1 of configuration file)",
                        required=False)

    parser.add_argument('--pdb_chains', nargs='+', dest='pdb_chains',
                        help="PDB chains to be extracted from PDB file (default: chains in properties of step 1)",
                        required=False)

    parser.add_argument('--mutation_list', nargs='+', dest='mutation_list',
                        help="List of mutations to be introduced in the protein (default: None, ex: A:Arg220Ala)",
                        required=False)

    parser.add_argument('--input_gro', dest='input_gro_path',
                        help="Input structure file ready to minimize (.gro). To provide an externally prepared system, use together with --input_top (default: None)",
                        required=False)

    parser.add_argument('--input_top', dest='input_top_path',
                        help="Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, use together with --input_gro (default: None)",
                        required=False)

    parser.add_argument('--fix_ss', action='store_true',
                        help="Add disulfide bonds to the protein. Use carefully! (default: False)",
                        required=False)

    parser.add_argument('--fix_backbone', action='store_true',
                        help="Add missing backbone atoms. Requires internet connection, PDB code and a MODELLER license key (default: False)",
                        required=False)

    args = parser.parse_args()

    # Default number of trajectories
    if args.num_trajs is None:
        num_trajs = 1
    else:
        num_trajs = int(args.num_trajs)

    # Check .gro structure and .zip topology are given together
    if (args.input_gro_path is None and args.input_top_path is not None) or (args.input_gro_path is not None and args.input_top_path is None):
        raise Exception("Both --input_gro and --input_top must be provided together")

    # Check .pdb structure and .gro/.zip topology are not given together
    if (args.input_pdb_path is not None and args.input_gro_path is not None):
        raise Exception("Both --input_pdb and --input_gro/--input_top are provided. Please provide only one of them")

    main_wf(configuration_path=args.config_path, setup_only=args.setup_only, num_trajs=num_trajs, output_path=args.output_path,
            input_pdb_path=args.input_pdb_path, pdb_chains=args.pdb_chains, mutation_list=args.mutation_list,
            input_gro_path=args.input_gro_path, input_top_path=args.input_top_path,
            fix_backbn=args.fix_backbone, fix_ss=args.fix_ss)