#!/usr/bin/env python3

# Conversion of the BioExcel building blocks Protein MD Setup Jupyter Notebook tutorial
# to a command line workflow with two files: Python Script and YAML input configuration file
# Example of Python Script (should be accompanied by a YAML input configuration file)

# Importing all the needed libraries
import time
import random
import argparse
import os

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
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_analysis.ambertools.cpptraj_rmsf import cpptraj_rmsf
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure

def set_gromacs_path(properties: dict, binary_path: str) -> None:
    """
    Set the path to the GROMACS binary for all steps using GROMACS.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        binary_path (str): Path to the GROMACS binary.
    """

    list_of_steps = ['step3_pdb2gmx', 'step4_editconf', 'step5_solvate', 'step6_grompp_genion', 'step7_genion',
                        'step8_grompp_min', 'step9_mdrun_min', 'step11_grompp_nvt', 'step12_mdrun_nvt',
                        'step14_grompp_npt', 'step15_mdrun_npt', 'step17_grompp_md', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['binary_path'] = binary_path

def set_mpi_path(properties: dict, mpi_bin: str, mpi_np: int) -> None:
    """
    Set the path to the MPI binary for all steps using MPI.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        mpi_bin (str): Path to the MPI binary.
        mpi_np (int): Number of processors to be used.
    """

    list_of_steps = ['step9_mdrun_min', 'step12_mdrun_nvt', 'step15_mdrun_npt', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['mpi_bin'] = mpi_bin
        properties[step]['mpi_np'] = mpi_np

def set_gpu_use(properties: dict, gpu_use: bool) -> None:
    """
    Set the use of GPU for all steps using GROMACS that support it.

    Inputs
    ------

        properties (dict): Dictionary containing the properties.
        gpu_use (bool): Whether to use GPU or not.
    """

    list_of_steps = ['step12_mdrun_nvt', 'step15_mdrun_npt', 'step18_mdrun_md']

    for step in list_of_steps:
        properties[step]['use_gpu'] = gpu_use

def set_general_properties(properties, conf, global_log) -> None:
    """
    Set all the additional global properties of this workflow, i.e. those properties included at the beginning of the YAML configuration file that
    are general to all steps and are not included already when the global properties are parsed.
    
    Inputs
    ------
    
        properties (dict): Dictionary containing the properties.
        conf (class settings.ConfReader): Configuration file reader.
    """
    
    # Enforce gromacs binary path for all steps using gromacs
    if conf.properties.get('binary_path'):
        global_log.info(f"Using GROMACS binary path: {conf.properties['binary_path']}")
        set_gromacs_path(properties, conf.properties['binary_path'])

    # Enforce mpi binary path for all steps using mpi
    if conf.properties.get('mpi_bin'):
        global_log.info(f"Using MPI binary path: {conf.properties['mpi_bin']}")
        set_mpi_path(properties, conf.properties['mpi_bin'], conf.properties.get('mpi_np'))

    # Enforce gpu use for all steps using gromacs that support it
    if conf.properties.get('use_gpu'):
        global_log.info(f"Using GPU for GROMACS steps")
        set_gpu_use(properties, conf.properties['use_gpu'])
    

def main_wf(configuration_path, setup_only, num_parts, num_replicas, output_path = None, input_pdb_path = None, pdb_chains = None,
            mutation_list = None, input_gro_path = None, input_top_path = None, fix_backbn = None, fix_ss = None, fix_amide_clashes = None, 
            his = None, nsteps = None, analysis = None):
    '''
    Main setup, mutation and MD run workflow with GROMACS. Can be used to retrieve a PDB, fix some defects of the structure,
    add specific mutations, prepare the system, minimize it, equilibrate it and finally do N production runs (either replicas or parts).

    Inputs
    ------

        configuration_path (str): path to YAML configuration file
        setup_only        (bool): whether to only setup the system or also run the simulations
        num_parts          (int): number of parts of the trajectory to launch in serial
        num_replicas       (int): number of replicas of the trajectory
        output_path        (str): (Optional) path to output folder
        input_pdb_path     (str): (Optional) path to input PDB file
        pdb_chains         (str): (Optional) chains to be extracted from the input PDB file
        mutation_list      (str): (Optional) list of mutations to be introduced in the structure
        input_gro_path     (str): (Optional) path to input structure file (.gro)
        input_top_path     (str): (Optional) path to input topology file (.zip)
        fix_backbn        (bool): (Optional) wether to add missing backbone atoms
        fix_ss            (bool): (Optional) wether to add disulfide bonds
        fix_amide_clashes (bool): (Optional) wether to flip clashing amides to relieve the clashes

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
    
    # Check num parts and num of replicas
    if (num_replicas is None) and (num_parts is None):
        # Default behavior: 1 replica
        num_replicas=1
    elif (num_replicas is not None) and (num_parts is not None):
        # Cannot set both num parts and num replicas
        global_log.error("Number of trajectories and replicas cannot be set at the same time")

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Set general properties for all steps
    set_general_properties(global_prop, conf, global_log)

    # If prepared structure is not provided
    if input_gro_path is None:

        # If input PDB is given as argument
        if input_pdb_path is not None:
            global_paths["step1_extractMolecule"]["input_structure_path"] = input_pdb_path

        # If chains are given as argument
        if pdb_chains is not None:
            global_prop["step1_extractMolecule"]["chains"] = pdb_chains

        # STEP 1: extract molecules of interest
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

        if fix_amide_clashes:
            # STEP 2 (G): Fix amides
            global_log.info("step2G_fixamides: fix clashing amides")
            fix_amides(**global_paths["step2G_fixamides"], properties=global_prop["step2G_fixamides"])
        else:
            if fix_ss:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2F_fixssbonds']['output_pdb_path']
            else:
                global_paths['step2H_fixchirality']['input_pdb_path'] = global_paths['step2E_fixsidechain']['output_pdb_path']

        # STEP 2 (H): Fix chirality
        global_log.info("step2H_fixchirality: fix chirality of residues")
        fix_chirality(**global_paths["step2H_fixchirality"], properties=global_prop["step2H_fixchirality"])

        # STEP 2 (I): renumber structure atoms and residues
        global_log.info("step2I_renumberstructure: renumber structure")
        renumber_structure(**global_paths["step2I_renumberstructure"], properties=global_prop["step2I_renumberstructure"])

        # NOTE: Check RMSF analysis step, not working - we need amber top
        
        # NOTE: Histidine protonation states come from pdb4amber which should be used within the WF!
        # STEP 3: add H atoms, generate coordinate (.gro) and topology (.top) file
        global_log.info("step3_pdb2gmx: Generate the topology")
        if his:
            global_prop["step3_pdb2gmx"]["his"]=his
 
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
    global_log.info("step13_temp_nvt: Compute temperature during NVT equilibration")
    gmx_energy(**global_paths["step13_temp_nvt"], properties=global_prop["step13_temp_nvt"])

    # STEP 14: NPT equilibration pre-processing
    global_log.info("step14_grompp_npt: Preprocess system pressure equilibration")
    grompp(**global_paths["step14_grompp_npt"], properties=global_prop["step14_grompp_npt"])

    # STEP 15: NPT equilibration
    global_log.info("step15_mdrun_npt: Execute system pressure equilibration")
    mdrun(**global_paths["step15_mdrun_npt"], properties=global_prop["step15_mdrun_npt"])

    # STEP 16: dump density and pressure evolution
    global_log.info("step16_density_npt: Compute Density & Pressure during NPT equilibration")
    gmx_energy(**global_paths["step16_density_npt"], properties=global_prop["step16_density_npt"])

    if num_replicas:
        # Folder names for replicas
        simulation_folders = [f"replica_{i}" for i in range(int(num_replicas))]
    elif num_parts:
        # Folder names for parts
        simulation_folders = [f"parts_{i}" for i in range(int(num_parts))]

    global_log.info(f"Number of parts: {num_parts}")
    global_log.info(f"Number of replicas: {num_replicas}")

    # Run each simulation (replica or part)
    traj_list = []
    for simulation in simulation_folders:
        
        traj_prop = conf.get_prop_dic(prefix=simulation)
        traj_paths = conf.get_paths_dic(prefix=simulation)

        # Set general properties for all steps
        set_general_properties(traj_prop, conf, global_log)

        # Update previous global paths
        traj_paths['step17_grompp_md']['input_gro_path'] = global_paths["step17_grompp_md"]['input_gro_path']
        traj_paths['step17_grompp_md']['input_cpt_path'] = global_paths["step17_grompp_md"]['input_cpt_path']
        traj_paths['step17_grompp_md']['input_top_zip_path'] = global_paths["step17_grompp_md"]['input_top_zip_path']
        traj_paths['step17_grompp_md']['input_ndx_path'] = global_paths["step17_grompp_md"]['input_ndx_path']
        
        # Enforce nsteps if provided
        if nsteps is not None:
            traj_prop['step17_grompp_md']['mdp']['nsteps']=int(nsteps)
            
        # Simulations are replicas
        if num_replicas:
            
            # Change seed and velocities for each replica
            traj_prop['step17_grompp_md']['mdp']['ld-seed'] = random.randint(1, 1000000)
            traj_prop['step17_grompp_md']['mdp']['continuation'] = 'no'
            traj_prop['step17_grompp_md']['mdp']['gen-vel'] = 'yes'

        # Simulations are parts of a single trajectory
        if num_parts:
            
            # Divide the number of steps by the number of parts
            traj_prop['step17_grompp_md']['mdp']['nsteps']=int(traj_prop['step17_grompp_md']['mdp']['nsteps']/int(num_parts))
            
            # For all parts except the first one, use the previous gro and cpt files
            if simulation != simulation_folders[0]:
                traj_paths['step17_grompp_md']['input_gro_path'] = previous_gro_path
                traj_paths['step17_grompp_md']['input_cpt_path'] = previous_cpt_path

        # STEP 17: free NPT production run pre-processing
        global_log.info(f"{simulation} >  step17_grompp_md: Preprocess free dynamics")
        grompp(**traj_paths['step17_grompp_md'], properties=traj_prop["step17_grompp_md"])

        # STEP 18: free NPT production run
        global_log.info(f"{simulation} >  step18_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**traj_paths['step18_mdrun_md'], properties=traj_prop['step18_mdrun_md'])
        
        # Append the trajectory to the list
        traj_list.append(traj_paths['step18_mdrun_md']['output_trr_path'])
        
        # Update the previous gro and cpt files
        previous_gro_path = traj_paths['step18_mdrun_md']['output_gro_path']
        previous_cpt_path = traj_paths['step18_mdrun_md']['output_cpt_path']

    # NOTE: If the workflow has to support both replicas and parts, this analysis should be done for each replica and part! 
    if int(analysis) == 1:
        
        # STEP 19: concatenate trajectories
        global_log.info("step19_trjcat: Concatenate trajectories")
        fu.zip_list(zip_file=global_paths["step19_trjcat"]['input_trj_zip_path'], file_list=traj_list)
        trjcat(global_paths["step19_trjcat"]["input_trj_zip_path"], global_paths["step19_trjcat"]["output_trj_path"], properties=global_prop["step19_trjcat"])

        # STEP 20: obtain dry trajectory
        global_log.info("step20_dry_trj: Obtain dry trajectory")
        gmx_trjconv_trj(**global_paths["step20_dry_trj"], properties=global_prop["step20_dry_trj"])

        # STEP 21: obtain dry structure
        global_log.info("step21_dry_str: Obtain dry structure")
        gmx_trjconv_str(**global_paths["step21_dry_str"], properties=global_prop["step21_dry_str"])

        #Remove the concatenated file that contains waters
        os.remove(global_paths["step19_trjcat"]["output_trj_path"])

        # STEP 22: dump RMSD with respect to equilibrated structure (first frame)
        global_log.info("step22_rmsfirst: Compute Root Mean Square deviation against equilibrated structure (first)")
        gmx_rms(**global_paths['step22_rmsfirst'], properties=global_prop['step22_rmsfirst'])

        # STEP 23: dump RMSD with respect to minimized structure
        global_log.info("step23_rmsexp: Compute Root Mean Square deviation against minimized structure (exp)")
        gmx_rms(**global_paths['step23_rmsexp'], properties=global_prop['step23_rmsexp'])

        # STEP 24: dump Radius of gyration
        global_log.info("step24_rgyr: Compute Radius of Gyration to measure the protein compactness during the free MD simulation")
        gmx_rgyr(**global_paths['step24_rgyr'], properties=global_prop['step24_rgyr'])

        # STEP 25: dump RMSF
        global_log.info("step25_rmsf: Compute Root Mean Square Fluctuation to measure the protein flexibility during the free MD simulation")
        cpptraj_rmsf(**global_paths['step25_rmsf'], properties=global_prop['step25_rmsf'])

        # STEP 26: image the trajectory
        global_log.info("step26_image_traj: Imaging the trajectory")
        gmx_image(**global_paths['step26_image_traj'], properties=global_prop['step26_image_traj'])

        # STEP 27: fit the trajectory
        global_log.info("step27_fit_traj: Fit the trajectory")
        gmx_image(**global_paths['step27_fit_traj'], properties=global_prop['step27_fit_traj'])

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

    parser.add_argument('--num_parts', dest='num_parts',
                        help="Number of parts of the trajectorie (default: 1)",
                        required=False)
    
    parser.add_argument('--num_replicas', dest='num_replicas',
                        help="Number of replicas (default: 1)",
                        required=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)

    parser.add_argument('--input_pdb', dest='input_pdb_path',
                        help="Input PDB file (default: input_structure_path in step 1 of configuration file)",
                        required=False)

    parser.add_argument('--his', dest='his',
                        help="Histidine protonation states list.",
                        required=False, type = str)

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

    parser.add_argument('--fix_amides', action='store_true', dest='fix_amide_clashes',
                        help="Flip clashing amides to relieve the clashes (default: False)",
                        required=False)

    parser.add_argument('--nsteps', dest='nsteps',
                        help="Number of steps of the simulation",
                        required=False)

    parser.add_argument('--analysis', dest='analysis',
                        help="Analysis",
                        required=False)

    args = parser.parse_args()

    # Check .gro structure and .zip topology are given together
    if (args.input_gro_path is None and args.input_top_path is not None) or (args.input_gro_path is not None and args.input_top_path is None):
        raise Exception("Both --input_gro and --input_top must be provided together")

    # Check .pdb structure and .gro/.zip topology are not given together
    if (args.input_pdb_path is not None and args.input_gro_path is not None):
        raise Exception("Both --input_pdb and --input_gro/--input_top are provided. Please provide only one of them")

    main_wf(configuration_path=args.config_path, setup_only=args.setup_only, num_parts=args.num_parts, num_replicas=args.num_replicas, output_path=args.output_path,
            input_pdb_path=args.input_pdb_path, pdb_chains=args.pdb_chains, mutation_list=args.mutation_list,
            input_gro_path=args.input_gro_path, input_top_path=args.input_top_path,
            fix_backbn=args.fix_backbone, fix_ss=args.fix_ss, fix_amide_clashes=args.fix_amide_clashes, his=args.his, nsteps=args.nsteps, analysis=args.analysis)