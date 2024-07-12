#!/shared/work/BiobbWorkflows/envs/biobb_sp_virtual_screening/bin/python

# Importing all the needed libraries
import os
import re
import time
import shutil
import argparse
from openbabel import pybel
from typing import Dict, List, Pattern, Tuple, Union, Optional

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_vs.utils.box import box
from biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_vs.vina.autodock_vina_run import autodock_vina_run
from biobb_chemistry.babelm.babel_convert import babel_convert
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens
from biobb_structure_utils.utils.extract_residues import extract_residues

def find_matching_str(pattern: Union[str, Pattern[str]], filepath: str) -> Optional[str]:
    '''
    Finds the first match of a regular expression pattern in a file.
    Returns the matching string or None if there is no match.

    Parameters
    ----------
    pattern : Union[str, Pattern[str]]
        Regular expression pattern to search in file lines.
    filepath : str
        File path to search in.

    Returns
    -------
    Optional[str]
        String matching the pattern or None if there is no match.
    '''
    try:
        with open(filepath, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    return match.group(1)
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    except IOError as e:
        print(f"An I/O error occurred: {e}")
    
    return None

def get_affinity(pdbqt_path: str) -> float:
    '''
    Find best binding affinity from pdbqt file from AutoDock Vina. 

    The first line with the following structure contains the best affinity:

    REMARK VINA RESULT:    -4.259      0.000      0.000

    Where the first number is the binding affinity.
    
    Inputs
    ------

        pdbqt_path: path to pdbqt file

    Outputs
    -------

        affinity: best binding affinity
    '''

    affinity_pattern=r'REMARK VINA RESULT:\s+(-?\d+\.\d+)'

    affinity = find_matching_str(affinity_pattern, pdbqt_path)

    if affinity:
        return float(affinity)
    else:
        return None

def create_summary(ranking: List[Tuple], properties: Dict, output_path: str) -> List[str]:
    '''
    Create file with ranking of ligands according to affinity

    Inputs
    ------

        ranking          : list with tuples -> (affinity, ligand name, smiles) ordered by affinity
        properties       : dictionary with properties of the step
        output_path      : path to output directory
    
    Output  
    ------

        top_ligand_names : list with names of top ligands
    '''

    # Define summary file name
    summary_filename = "top_ligands.csv"

    # Define summary file path
    summary_path = os.path.join(output_path, summary_filename)

    # Find number of top ligands to save
    ranking_length = min(properties['num_top_ligands'], len(ranking)) 

    # Extract top ligands (ranking is already sorted)
    top_ligands = ranking[:ranking_length]
    top_ligand_names = []

    # Create step folder if it does not exist
    if not os.path.exists(properties['path']):
        os.makedirs(properties['path'])
    
    # Create summary file with top ligands
    with open(summary_path, 'w') as file:

        # Write header
        file.write("Rank,Affinity,Name,Identifier \n")

        # For each ligand
        for rank, affinity_tuple in enumerate(top_ligands):

            affinity, ligand_name, ligand_ID = affinity_tuple

            # Write line
            file.write(f"{rank+1},{affinity},{ligand_name},{ligand_ID}\n")

            # Add ligand name to list
            top_ligand_names.append(ligand_name)

    return top_ligand_names

def validate_step(*output_paths: str) -> bool:
    '''
    Check all output files exist and are not empty
    
    Inputs
    ------

        *output_paths (str): variable number of paths to output file/s

    Output
    ------

        validation_result (bool): result of validation
    '''

    # Initialize value 
    validation_result = True

    # Check existence of files
    for path in output_paths:
        validation_result = validation_result and os.path.exists(path)

    # Check files are not empty if they exist
    if (validation_result):

        for path in output_paths:
            file_not_empty = os.stat(path).st_size > 0
            validation_result = validation_result and file_not_empty

    return validation_result

def read_ligand_lib(ligand_lib_path: str) -> Tuple[List[str], List[str]]:
    '''
    Read all ligand identifiers from ligand library file. The expected format is one of the following:

    - Format 1:

            ligand1_id \n
            ligand2_id

    - Format 2:

            ligand1_id  ligand1_name \n
            ligand2_id  ligand2_name

    Where ligand_id is a unique identifier (e.g. SMILES) and name_ligand is a string with the ligand name
    
    Inputs
    ------

        ligand_lib_path (str): path to ligand library file

    Output
    ------

        ligand_ids     : list of ligand identifiers
        ligand_names   : list with ligand names
    '''

    ligand_ids = []
    ligand_names = []

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Process every line
        for index, line in enumerate(ligand_lib):

            line = line.split()

            # Append ligand ID to list
            ligand_ids.append(line[0])

            # If there is no name, use index as name
            if len(line)>1:
                ligand_names.append(line[1])
            else:
                ligand_names.append(str(index))

        # If there are no ligands, raise an error
        if len(ligand_ids) == 0:
            raise ValueError(f"No ligands found in ligand library file {ligand_lib_path}")
        
    return ligand_ids, ligand_names

def write_smiles(smiles: str, smiles_path: str):
    '''
    Writes a SMILES code into a file. If the file exists, it will be overwritten.

    Inputs
    ------
        smiles         :  SMILES code
        smiles_path    :  smiles file path
    '''

    # Save SMILES in tmp file inside step_path, overwrite if exists
    smiles_tmp_file = open(smiles_path, 'w')
    smiles_tmp_file.write(smiles)
    smiles_tmp_file.close()

def get_ranking(ligand_ids: List, ligand_names: List, autodock_vina_paths: Dict, output_path: Dict) -> List[Tuple]:   
    """
    Takes the name of each ligand and finds the best affinity given by AutoDock Vina from the output pdbqt file.
    Returns a list of tuples ordered by affinity: (affinity, ligand_name, ligand_ids) 

    Inputs
    ------

        ligand_ids          : list of ligand ids
        ligand_names        : list of ligand names
        autodock_vina_paths : paths dictionary for autodock vina step
        output_path         : path to output directory
    
    Output
    ------

        ranking         : list of tuples ordered by affinity: (affinity, ligand_name, ligand_ids) 
    """

    # Generic path to pdbqt file with poses
    output_pdbqt_path = autodock_vina_paths['output_pdbqt_path']

    # Name of the pdbqt file
    pdbqt_filename = os.path.basename(output_pdbqt_path)

    # Name of the autodock vina step folder
    autodock_step_name = os.path.basename(os.path.dirname(output_pdbqt_path))

    # List where best affinity for each ligand will be stored
    ranking = []
    
    for ligand_id, ligand_name in zip(ligand_ids, ligand_names):

        # Find path to ligand subfolder
        ligand_folder = os.path.join(output_path, ligand_name)

        # Find path to autodock vina step
        vina_folder = os.path.join(ligand_folder, autodock_step_name)

        # Find path to output pdbqt file with poses
        pdbqt_path = os.path.join(vina_folder, pdbqt_filename)

        # Find best affinity among different poses
        affinity = get_affinity(pdbqt_path = pdbqt_path)

        if affinity:
            ranking.append((affinity, ligand_name, ligand_id))

    # Sort list according to affinity
    ranking = sorted(ranking)

    return ranking

def clean_output(ligand_names: List, output_path: str):
    """
    Removes all ligand sub folders in the output folder

    Inputs
    ------

        ligand_names    : list of ligand names
        output_path     : path to output directory
    """

    # Remove all ligand subdirectories
    for name in ligand_names:
        ligand_path = os.path.join(output_path, name)

        if os.path.exists(ligand_path):
            shutil.rmtree(ligand_path)
    
def check_arguments(global_log, global_paths, global_prop, ligand_lib_path, structure_path, input_pockets_zip, dock_to_residues):
    """
    Check the arguments provided by the user and values of configuration file
    """

    # Check the ligand library path exists and it's a file
    if not os.path.exists(ligand_lib_path):
        global_log.warning(f"WARNING: Ligand library file {ligand_lib_path} does not exist")
    elif not os.path.isfile(ligand_lib_path):
        global_log.warning(f"WARNING: Ligand library path {ligand_lib_path} is not a file")

    # Check we have a structure file
    config_structure_path = global_paths['step3_str_check_add_hydrogens']['input_structure_path']
    if structure_path is None and not os.path.exists(config_structure_path):
        global_log.warning(f"WARNING: Structure file {config_structure_path} does not exist")
    elif structure_path is not None and not os.path.exists(structure_path):
        global_log.warning(f"WARNING: Structure file {structure_path} does not exist")

    # Check we have a pockets zip file
    config_pockets_zip = global_paths['step1_fpocket_select']['input_pockets_zip']
    if input_pockets_zip is None and not os.path.exists(config_pockets_zip):
        global_log.warning(f"WARNING: Pockets zip file {config_pockets_zip} does not exist")
    elif input_pockets_zip is not None and not os.path.exists(input_pockets_zip):
        global_log.warning(f"WARNING: Pockets zip file {input_pockets_zip} does not exist")

    # Check size of box
    if dock_to_residues:
        if global_prop['step2_box']['offset'] > 5:
            global_log.warning(f"WARNING: box offset is {global_prop['step2_box']['offset']} angstroms. This may be unnecessarily large when docking to residues surrounding the binding site. Consider using a smaller value to improve performance.")

def main_wf(configuration_path: str, ligand_lib_path: str, structure_path: str, input_pockets_zip: str, pocket: str, output_path: str, num_top_ligands: int, keep_poses: bool, dock_to_residues: bool, cpus: int, exhaustiveness: int, debug: bool) -> Tuple[Dict, Dict]:
    '''
    Main VS workflow. This workflow takes a ligand library, a pocket (defined by the output of a cavity analysis or some residues) 
    and a receptor to screen the cavity using the ligand library (with AutoDock).

    Inputs
    ------

        configuration_path   : path to input.yml 
        ligand_lib_path      : path to ligand library. Either a SMILES file or a SDF file
        structure_path       : path to receptor structure
        input_pockets_zip    : path to input pockets zip file
        pocket               : pocket name
        output_path          : path to output directory
        num_top_ligands      : number of top ligands to be saved
        keep_poses           : keep poses of top ligands
        dock_to_residues     : dock to residues instead of cavity
        cpus                 : number of cpus to use for each docking
        exhaustiveness       : exhaustiveness of the docking
        debug                : keep intermediate files for debugging

    Outputs
    -------

        global_paths    : dictionary with all workflow paths
        global_prop     : dictionary with all workflow properties
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
    # Dividing it in global properties and global paths
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Check arguments
    check_arguments(global_log, global_paths, global_prop, ligand_lib_path, structure_path, input_pockets_zip, dock_to_residues)
                    
    # Enforce input_pockets_zip if provided
    if input_pockets_zip is not None:
        global_paths['step1_fpocket_select']['input_pockets_zip'] = input_pockets_zip
    
    # Enforce pocket selection if provided
    if pocket is not None:
        global_prop['step1_fpocket_select']['pocket'] = pocket

    # Enforce structure_path if provided
    if structure_path is not None:
        global_paths['step1b_extract_residues']['input_structure_path'] = structure_path
        global_paths['step3_str_check_add_hydrogens']['input_structure_path'] = structure_path
    else:
        structure_path = global_paths['step3_str_check_add_hydrogens']['input_structure_path']
    
    # Enforce num_top_ligands if specified
    if num_top_ligands is not None:
        global_prop['step6_top_ligands']['num_top_ligands'] = int(num_top_ligands)

    # STEP 1: Select pocket or extract residues
    if dock_to_residues:
        # Extract residues from structure
        global_log.info("step1b_extract_residues: Extracting residues from structure")
        extract_residues(**global_paths["step1b_extract_residues"], properties=global_prop["step1b_extract_residues"])

        # Modify step2_box paths to use residues
        global_paths['step2_box']['input_pdb_path'] = global_paths['step1b_extract_residues']['output_residues_path']

    else:
        # Pocket selection from filtered list 
        global_log.info("step1_fpocket_select: Extract pocket cavity")
        fpocket_select(**global_paths["step1_fpocket_select"], properties=global_prop["step1_fpocket_select"])

    # STEP 2: Generate box around selected cavity or residues
    global_log.info("step2_box: Generating cavity box")
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])

    # STEP 3: Prepare target protein for docking 
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 

    docking_start_time = time.time()

    # STEP 4-5: Prepare ligand pdbqt and dock with AutoDock Vina

    # Option 1: SDF library with protonated ligands
    if ligand_lib_path.endswith('.sdf'):

        ligand_names = []
        ligand_ids = []

        global_log.info(f"Reading ligand library in SDF format")
        ligand_supplier = pybel.readfile('sdf', ligand_lib_path)

        for index, ligand in enumerate(ligand_supplier, start=0):

            # Create unique name from title and index
            ligand_name = f"{ligand.title}_{index}" if ligand.title else f"ligand_{index}"
            ligand_id = ligand.title if ligand.title else ""

            # Add ligand name to properties and paths
            ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            ligand_paths = conf.get_paths_dic(prefix=ligand_name)

            ligand_names.append(ligand_name)
            ligand_ids.append(ligand_id)

            # Create ligand subfolder
            ligand_folder = os.path.join(output_path, ligand_name)
            if not os.path.exists(ligand_folder):
                os.makedirs(ligand_folder)

            # Create step folder
            if not os.path.exists(ligand_prop['step4b_babel_convert']['path']):
                os.makedirs(ligand_prop['step4b_babel_convert']['path'])

            # Write ligand to sdf file # NOTE: writing the pdbqt file directly with pybel discards hydrogens
            ligand.write(format='sdf', filename=ligand_paths['step4b_babel_convert']['input_path'])

            # STEP 4: Convert ligand from sdf to pdbqt without adding hydrogens # NOTE format output option -xh is needed
            global_log.info("step4b_babel_convert: Convert ligand to pdbqt format")
            try:
                babel_convert(**ligand_paths['step4b_babel_convert'], properties = ligand_prop["step4b_babel_convert"])
                lastStep_successful = validate_step(ligand_paths['step4b_babel_convert']['output_path'])
            except:
                global_log.info(f"step4b_babel_convert: Open Babel failed to convert ligand {ligand_name} to pdbqt format")
                lastStep_successful = False

            # STEP 5: AutoDock vina
            if lastStep_successful:
            
                # Enforce cpus if specified
                if cpus is not None:
                    ligand_prop['step5_autodock_vina_run']['cpu'] = int(cpus)
                
                # Enforce exhaustiveness if specified
                if exhaustiveness is not None:
                    ligand_prop['step5_autodock_vina_run']['exhaustiveness'] = int(exhaustiveness)
                
                # Update common paths
                ligand_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path'] = global_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path']
                ligand_paths['step5_autodock_vina_run']['input_box_path'] = global_paths['step5_autodock_vina_run']['input_box_path']
                ligand_paths['step5_autodock_vina_run']['input_ligand_pdbqt_path'] = ligand_paths['step4b_babel_convert']['output_path']

                try:
                    global_log.info("step5_autodock_vina_run: Docking the ligand")
                    autodock_vina_run(**ligand_paths['step5_autodock_vina_run'], properties=ligand_prop["step5_autodock_vina_run"])
                    lastStep_successful = validate_step(ligand_paths['step5_autodock_vina_run']['output_log_path'],
                                                        ligand_paths['step5_autodock_vina_run']['output_pdbqt_path'])
                except:
                    global_log.info(f"step5_autodock_vina_run: Autodock Vina failed to dock ligand {ligand_name}")

    # Option 2: SMILES library with ligands to be prepared
    elif ligand_lib_path.endswith('.smi'):

        global_log.info(f"Reading ligand library in SMILES format")
        ligand_ids, ligand_names = read_ligand_lib(ligand_lib_path)

        for ligand_id, ligand_name in zip(ligand_ids, ligand_names):

            # Add ligand name to properties and paths
            ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            ligand_paths = conf.get_paths_dic(prefix=ligand_name)

            # Create ligand subfolder
            ligand_folder = os.path.join(output_path, ligand_name)
            if not os.path.exists(ligand_folder):
                os.makedirs(ligand_folder)

            # Create step folder
            if not os.path.exists(ligand_prop['step4_babel_protonate']['path']):
                os.makedirs(ligand_prop['step4_babel_protonate']['path'])

            # Write smiles to file
            write_smiles(smiles = ligand_id, smiles_path = ligand_paths['step4_babel_protonate']['input_path'])

            # STEP 4: Convert ligand from smiles to pdbqt adding hydrogens at a certain pH
            global_log.info("step4_babel_protonate: Prepare ligand for docking")
            try:
                babel_convert(**ligand_paths['step4_babel_protonate'], properties = ligand_prop["step4_babel_protonate"])
                lastStep_successful = validate_step(ligand_paths['step4_babel_protonate']['output_path'])
            except:
                global_log.info(f"step4_babel_protonate: Open Babel failed to convert ligand {ligand_name} to pdbqt format")
                lastStep_successful = False
        
            # STEP 5: AutoDock vina
            if lastStep_successful:

                # Enforce cpus if specified
                if cpus is not None:
                    ligand_prop['step5_autodock_vina_run']['cpu'] = int(cpus)

                # Enforce exhaustiveness if specified
                if exhaustiveness is not None:
                    ligand_prop['step5_autodock_vina_run']['exhaustiveness'] = int(exhaustiveness)

                # Update common paths
                ligand_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path'] = global_paths['step5_autodock_vina_run']['input_receptor_pdbqt_path']
                ligand_paths['step5_autodock_vina_run']['input_box_path'] = global_paths['step5_autodock_vina_run']['input_box_path']

                try:
                    global_log.info("step5_autodock_vina_run: Docking the ligand")            
                    autodock_vina_run(**ligand_paths['step5_autodock_vina_run'], properties=ligand_prop["step5_autodock_vina_run"])
                    lastStep_successful = validate_step(ligand_paths['step5_autodock_vina_run']['output_log_path'], 
                                                        ligand_paths['step5_autodock_vina_run']['output_pdbqt_path'])
                except:
                    global_log.info(f"step5_autodock_vina_run: Autodock Vina failed to dock ligand {ligand_name}")

    else:

        global_log.error(f"Ligand library file {ligand_lib_path} should be in SDF or SMILES format")
        return

    # Rank ligands: find the best affinity for each ligand
    ranking = get_ranking(ligand_ids, ligand_names, global_paths['step5_autodock_vina_run'], output_path)

    # STEP 6: Find top ligands and create csv file with ranking
    global_log.info("step6_top_ligands: create ranking and save poses for top ligands")  
    top_ligand_names = create_summary(ranking, global_prop['step6_top_ligands'], output_path)

    # STEP 7: extract poses for top ligands if requested
    if keep_poses:

        # Iterate over top ligands
        for ligand_name in top_ligand_names:

            # Add ligand name to properties and paths
            top_ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            top_ligand_paths = conf.get_paths_dic(prefix=ligand_name)

            try:
                # Convert pose from pdbqt to pdb
                global_log.info("step7_babel_prepare_pose: Converting ligand pose to PDB format")    
                babel_convert(**top_ligand_paths['step7_babel_prepare_pose'], properties=top_ligand_prop["step7_babel_prepare_pose"])

                # Move pose to final location
                # Pose path inside ligand subfolder
                pose_path = top_ligand_paths['step7_babel_prepare_pose']['output_path']
                # New pose path in output folder
                new_pose_path = os.path.join(global_prop['step6_top_ligands']['path'], f"{ligand_name}_poses.pdb")
                # Move pose to new location 
                shutil.move(pose_path, new_pose_path)

            except:
                global_log.info(f"step7_babel_prepare_pose: Open Babel failed to convert pose for ligand {ligand_name} to PDB format")
    
    # Show success rate of screening
    success_rate = round(len(ranking)/len(ligand_names)*100, 2)
    global_log.info(f"Success rate: {success_rate}%")

    if not debug:
        # Clean up the output folder 
        clean_output(ligand_names, output_path)

    # Save structure path in output_path
    shutil.copy(structure_path, os.path.join(output_path, 'receptor.pdb'))

    # Timing information
    elapsed_time = time.time() - start_time
    docking_elapsed_time = time.time() - docking_start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow name: Virtual Screening')
    global_log.info('  Output path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('  Ligand library: %s' % ligand_lib_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('Docking time: %.1f minutes' % (docking_elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")
    
    parser.add_argument('-c', '--config', dest='config_path', type=str,
                        help="Configuration file (YAML)",
                        required=True)

    parser.add_argument('-lib', '--ligand_lib', dest='ligand_lib', type=str,
                        help="Path to file with the ligand library. The format should be SMILES (.smi) or SDF (.sdf). For .smi files, one ligand per line is expected: 'smiles name'. For sdf files, the file may contain one or more ligands.",
                        required=True)
    
    parser.add_argument('-s', '--structure_path', dest='structure_path', type=str,
                        help="Path to file with target structure (PDB format)",
                        required=False)
    
    parser.add_argument('-pz', '--input_pockets_zip', dest='input_pockets_zip', type=str,
                        help="Path to file with pockets in a zip file",
                        required=False)

    parser.add_argument('-p', '--pocket', dest='pocket', type=int,
                        help="Pocket number to be used (default: 1)",
                        required=False)

    parser.add_argument('-o', '--output', dest='output_path', type=str,
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    parser.add_argument('-nl', '--num_top_ligands', dest='num_top_ligands', type=int,
                        help="Number of top ligands to be saved (default: corresponding value in YAML config file)",
                        required=False)

    parser.add_argument('-kp', '--keep_poses', dest='keep_poses', action='store_true',
                        help="Save docking poses for top ligands (default: False)",
                        required=False)

    parser.add_argument('-dr', '--dock_to_residues', dest='dock_to_residues', action='store_true',
                        help="Dock to residues instead of pocket. Define the docking box using a set of residues instead of a pocket. (default: False)",
                        required=False)
    
    parser.add_argument('-cpus', '--cpus', dest='cpus', type=int, 
                        help="Number of CPUs to use for each docking (default: 1)",
                        required=False)
    
    parser.add_argument('-ex', '--exhaustiveness', dest='exhaustiveness', type=int,
                        help="Exhaustiveness of the docking (default: 8)",
                        required=False)
    
    parser.add_argument('-d', '--debug', dest='debug', action='store_true',
                        help="Keep intermediate files for debugging (default: False)",
                        required=False)
    
    args = parser.parse_args()

    main_wf(configuration_path = args.config_path, 
            ligand_lib_path = args.ligand_lib,
            structure_path = args.structure_path,
            input_pockets_zip = args.input_pockets_zip,
            pocket = args.pocket,
            output_path = args.output_path,
            num_top_ligands = args.num_top_ligands,
            keep_poses = args.keep_poses,
            dock_to_residues = args.dock_to_residues,
            cpus = args.cpus,
            exhaustiveness = args.exhaustiveness,
            debug = args.debug)
