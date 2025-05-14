#!/usr/bin/env python3
import os
import time
import shutil
import argparse
from pathlib import Path
from Bio.PDB import PDBParser
from typing import List, Union, Dict, Literal

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_amber.leap.leap_gen_top import leap_gen_top
from biobb_chemistry.babelm.babel_minimize import babel_minimize
from biobb_chemistry.acpype.acpype_params_ac import acpype_params_ac
from biobb_chemistry.acpype.acpype_params_gmx import acpype_params_gmx
from biobb_chemistry.babelm.babel_add_hydrogens import babel_add_hydrogens
from biobb_chemistry.ambertools.reduce_add_hydrogens import reduce_add_hydrogens
from biobb_structure_utils.utils.extract_heteroatoms import extract_heteroatoms
from biobb_chemistry.acpype.acpype_convert_amber_to_gmx import acpype_convert_amber_to_gmx

def get_selected_ligands(pdb_path: str, selected_ligand_names:  Union[List[str], None], selected_chains: List[str], selected_model: int, global_log) -> List[Dict[str, str]]:
    """
    Retrieve the selected ligands from the chains and model of interest of a PDB file. If selected_ligand_names is None, all ligands are extracted.
    
    Inputs
    ------
    
        pdb_path               : Path to the PDB file.
        selected_ligand_names  : List of ligand names to be extracted. If None, all ligands are extracted.
        selected_chains        : List of chain IDs to extract ligands from.
        selected_model         : Model number to extract ligands from.
        global_log             : Logger object for logging messages.
    
    Returns
    -------
    
        list: List with dictionaries defining the ligands present in the chains and model of interest, 
        empty if no parameterized ligand is found. Each dictionary has the following keys:
        
            Example 
            
                ligands = [
                    {
                        'name': 'ZZ7',
                        'res_id': '302',
                        'chain': 'A',
                        'model': 1
                    },
                    {
                        'name': 'FLP',
                        'res_id': '402',
                        'chain': 'B',
                        'model': 1
                    }
                ]
    """
    
    # Initialize the list of ligands
    ligands = []
    
    # Check if the PDB file exists
    if not os.path.exists(pdb_path):
        global_log.error(f"File {pdb_path} not found")
        return ligands
    
    # Search for ligands in the chains and model of interest
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)
    for model in structure:
        for chain in model:
            if (chain.get_id() in selected_chains) and (model.id == selected_model):
                for residue in chain:
                    # If residue is made of heteroatoms (atoms of ligands), e.g. of residue.id: ('H_HEM', 500, ' ')
                    # Heteroatom residues have a hetero flag in their residue.id[0] (e.g. H_XXX for a ligand XXX)
                    heteroatom_id = residue.id[0]
                    if heteroatom_id.startswith('H_'):
                        
                        ligand_name = residue.get_resname()
                        
                        if (selected_ligand_names is None) or (ligand_name in selected_ligand_names):
                            
                            res_id = str(residue.id[1])
                            
                            ligand_info = {
                                'name': ligand_name,
                                'res_id': str(res_id),
                                'chain': chain.get_id(),
                                'model': str(model.id+1)
                            }
                            
                            ligands.append(ligand_info)
    return ligands

def get_parameter_sets(custom_parameters_path: Union[str, None], global_log) -> Dict[str, Dict[str, str]]:
    """
    Retrieve the custom parameter sets for one or more ligands.
    
    Inputs
    ------
    
        custom_parameters_path  : Path to the folder with custom parameter sets for one or more ligands.
        global_log              : Logger object for logging messages.
    
    Returns
    -------
    
        dict: Dictionary with the custom parameter sets for one or more ligands. Each key is the ligand name and each value is a dictionary with the following keys:
        
            Example 
            
                parameter_sets = {
                    'ZZ7': {
                        'frcmod': 'path/to/ZZ7.frcmod',
                        'prep': 'path/to/ZZ7.prep'
                    },
                    'FLP': {
                        'frcmod': 'path/to/FLP.frcmod',
                        'prep': 'path/to/FLP.prep'
                    }
                }
    """
    
    # Initialize the dictionary of parameter sets
    parameter_sets = {}
    
    # Check if custom parameters are provided
    if custom_parameters_path is None:
        return parameter_sets
    
    # Check if the custom parameters folder exists
    if not os.path.exists(custom_parameters_path):
        global_log.error(f"Folder {custom_parameters_path} not found")
        return parameter_sets
    
    # Search for custom parameter sets in the folder
    for file_name in os.listdir(custom_parameters_path):  

        # Check if the file is a .frcmod or .prep file
        if file_name.endswith('.frcmod'):
            ligand_name = file_name.replace('.frcmod', '')
            if ligand_name not in parameter_sets:
                parameter_sets[ligand_name] = {}
            parameter_sets[ligand_name]['frcmod'] = os.path.join(custom_parameters_path, file_name)
        elif file_name.endswith('.prep'):
            ligand_name = file_name.replace('.prep', '')
            if ligand_name not in parameter_sets:
                parameter_sets[ligand_name] = {}
            parameter_sets[ligand_name]['prep'] = os.path.join(custom_parameters_path, file_name)

    # Check every ligand has both .frcmod and .prep files
    for ligand_name, parameter_set in parameter_sets.items():
        if ('frcmod' not in parameter_set) or ('prep' not in parameter_set):
            global_log.error(f"Missing .frcmod or .prep file for ligand {ligand_name}")
            return {}

    # Check if there are no custom parameter sets
    if len(parameter_sets) == 0:
        global_log.error("No custom parameter sets found in the folder")
        return {}

    # Log the custom parameter sets found
    global_log.info(f"Found {len(parameter_sets)} custom parameter sets for ligands: {', '.join(parameter_sets.keys())}")

    return parameter_sets

def gmx_top2itp(top_path: str, itp_path: str):
    """
    Read the input GROMACS topology from top_path. 
    Removes any [ defaults ] directive present, as they can only be present in the main topology file.
    Removes any [ molecules ] directive present, as they can only be present in the main topology file.
    Copies the rest of the topology file to itp_path.
    
    Inputs
    ------
    
        top_path : Path to the input GROMACS topology file.
        itp_path : Path to the output GROMACS itp file.
    """
    
    # Read the input topology file
    with open(top_path, 'r') as f:
        lines = f.readlines()
    
    
    # Remove any defaults and molecules directives
    new_lines = []
    reading_defaults = False
    reading_molecules = False
    for line in lines:
        
        # Mark the end of the defaults section
        if reading_defaults and line.startswith("["):
            reading_defaults = False
            
        # Mark the end of the molecules section
        if reading_molecules and line.startswith("["):
            reading_molecules = False
            
        # Mark the beginning of the defaults section
        if line.startswith("[ defaults ]"):
            reading_defaults = True 
        
        # Mark the beginning of the molecules section
        if line.startswith("[ molecules ]"):
            reading_molecules = True

        # Add the line if not reading the defaults section or the molecules section
        if not reading_defaults and not reading_molecules:
            new_lines.append(line)
    
    # Write the new topology file
    with open(itp_path, 'w') as f:
        f.writelines(new_lines)
        
def copy_out_files(file_paths: List[str], ligand_name: str, output_folder: str):
    '''
    Copy coordinate and topology files to the final output folder changing the file names to the ligand name.
    If the file is a GROMACS topology (.top) file, it is converted to an .itp file.
    
    Inputs
    ------
    
        file_paths   : List of file paths to copy.
        ligand_name  : Name of the ligand.
        output_folder: Path to the output folder.
    '''
    
    for path in file_paths:
    
        file_extension = Path(path).suffix
        
        if file_extension == '.top':
            new_file_path = os.path.join(output_folder, f"{ligand_name}.itp")
            gmx_top2itp(path, new_file_path)
        else:
            new_file_path = os.path.join(output_folder, f"{ligand_name}{file_extension}")
            shutil.copyfile(path, new_file_path)

# YML construction
def config_contents() -> str:
    """
    Returns the contents of the YAML configuration file as a string.
    
    The YAML file contains the configuration for the protein preparation workflow.
    
    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    return f""" 
# Global properties (common for all steps)
global_properties:                                                                               # Wether to use GPU support or not
  working_dir_path: output                                                                          # Workflow default output directory
  can_write_console_log: False                                                                      # Verbose writing of log information
  restart: True                                                                                     # Skip steps already performed
  remove_tmp: True  

# Step 1: extract heteroatoms from input PDB file
step1_ligand_extraction:
  tool: extract_heteroatoms
  paths:
    input_structure_path: /path/to/input              # Will be set by the workflow
    output_heteroatom_path: ligand.pdb
  properties:
    heteroatoms: null                                  # Will be set by the workflow

step2A_leap_gen_top:
  tool: leap_gen_top
  paths:
    input_pdb_path: dependency/step1_ligand_extraction/output_heteroatom_path
    input_frcmod_path: path/to/frcmod.zip             # Will be set by the workflow
    input_prep_path: path/to/prep.zip                 # Will be set by the workflow
    output_pdb_path: cofactors.pdb
    output_top_path: cofactors.prmtop
    output_crd_path: cofactors.inpcrd
  properties:
    forcefield: ["protein.ff14SB"]                    # Will be set by the workflow

step3A_amber_to_gmx:
  tool: acpype_convert_amber_to_gmx
  paths:
    input_crd_path: dependency/step2A_leap_gen_top/output_crd_path
    input_top_path: dependency/step2A_leap_gen_top/output_top_path
    output_path_gro: cofactors.gro
    output_path_top: cofactors.top
  properties:
    basename: cofactors

step2B_ambertools_reduce:
  tool: reduce_add_hydrogens
  paths:
    input_path: dependency/step1_ligand_extraction/output_heteroatom_path
    output_path: ligand_reduced.pdb
  properties:
    charges: True

step2B_obabel_reduce:
  tool: babel_add_hydrogens
  paths:
    input_path: dependency/step1_ligand_extraction/output_heteroatom_path
    output_path: ligand_reduced.pdb
  properties:
    coordinates: 3,
    ph: 7.4

step3B_babel_minimize:
  tool: babel_minimize
  paths:
    input_path:  dependency/step2B_ambertools_reduce/output_path
    output_path: ligand_minimized.pdb
  properties:
    criteria: 1e-6
    method: sd
    force_field: GAFF
    hydrogens: False

step4B_acpype_params_gmx:
  tool: acpype_params_gmx
  paths:
    input_path:  dependency/step3B_babel_minimize/output_path
    output_path_gro: output.params.gro
    output_path_itp: output.params.itp
    output_path_top: output.params.top
  properties:
    basename: biobb_GMX_LP

step4B_acpype_params_ac:
  tool: acpype_params_ac
  paths:
    input_path:  dependency/step3B_babel_minimize/output_path
    output_path_frcmod: output.frcmod
    output_path_inpcrd: output.inpcrd
    output_path_lib: output.lib
    output_path_prmtop: output.prmtop
  properties:
    basename: biobb_AC_LP
"""

def create_config_file(config_path: str) -> None:
    """
    Create a YAML configuration file for the workflow if needed.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file to be created.
    
    Returns
    -------
    None
    """
    
    # Check if the file already exists
    if os.path.exists(config_path):
        print(f"Configuration file already exists at {config_path}.")
        return
    
    # Write the contents to the file
    with open(config_path, 'w') as f:
        f.write(config_contents())

# Main workflow
def main(configuration_path: str, input_pdb: str, forcefields: List[str] = ['protein.ff14SB', 'DNA.bsc1', 'gaff'], ligand_names: Union[List[str], None] = None, 
         charges: Union[List[str], None] = None, chains: List[str] = ['A'], model: int = 0, format: Literal['gromacs', 'amber'] = 'gromacs', custom_parameters: Union[str, None] = None,
         protonation_tool: Literal['ambertools', 'obabel', 'none'] = 'ambertools', skip_min: bool = False, 
         output_top_path: Union[str, None] = None, output_path: Union[str, None] = None):
    '''
    Ligand parameterization workflow using BioExcel Building Blocks.

    Inputs
    ------
    
        configuration_path : Path to the configuration file (YAML).
        input_pdb          : Path to the input PDB file with the ligands to parameterize.
        forcefields        : (Optional) List of paths or file names of force fields to use in the parameterization. Only used within LEaP if custom parameters are also provided.
        ligand_names       : (Optional) List of ligand names in the PDB file to parameterize. By default, all ligands are parameterized.
        charges            : (Optional) List of charges for the ligands to parameterize. Only used when no custom parameters are given and the ligand is parameterized using GAFF and antechamber through acpype. By default acpype will guess the charge based on the protonation state.
        chains             : (Optional) Chain ID of the ligands to parameterize. Default: A.
        model              : (Optional) Model number of the ligands to parameterize. Default: 1.
        format             : (Optional) Format of the output topology files. Options: gromacs, amber. Default: gromacs.
        custom_parameters  : (Optional) Path to folder with custom parameter sets for one or more ligands (.prep and .frcmod files with the ligand name).
        protonation_tool   : (Optional) Protonation tool to use. Options: ambertools, obabel. Default: ambertools.
        skip_min           : (Optional) Skip the minimization step. Default: False.
        output_top_path    : (Optional) Output path for the folder with topologies and coordinate files.
        output_path        : (Optional) Output path. Default: working_dir_path in YAML config file.
        
    Outputs
    -------

        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties
    '''
    
    ###########################
    # Workflow initialization #
    ###########################
    
    start_time = time.time()
    
    if configuration_path is None:
        # Create a default configuration file
        configuration_path = "config.yml"
        create_config_file(configuration_path)
    
    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(configuration_path)
    
    # Enforce output_path if provided
    if output_path is not None:
        output_path = fu.get_working_dir_path(output_path, restart = conf.properties.get('restart', 'False'))
        conf.working_dir_path = output_path
    else:
        output_path = conf.get_working_dir_path()
        
    if output_top_path is None:
        output_top_path = os.path.join(output_path, 'topologies')
    os.makedirs(output_top_path, exist_ok=True)
        
    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()
    
    # Extract the ligands from the input PDB file
    global_log.info(f"Searching selected ligands in the input PDB file: {' '.join(ligand_names)}")
    selected_ligands = get_selected_ligands(input_pdb, ligand_names, chains, model, global_log)
    global_log.info(f"Found {len(selected_ligands)} ligands in the input PDB file: {', '.join([ligand['name'] for ligand in selected_ligands])}")
    
    if len(selected_ligands) == 0:
        global_log.error("No ligands found in the input PDB file, review your PDB and selection criteria")
        return
    
    # Check the format of the output topology files
    if format:
        if format not in ['gromacs', 'amber']:
            global_log.error(f"Invalid format {format}. Options: gromacs, amber")
            return
    
    # Check the protonation tool to use
    if protonation_tool:
        if protonation_tool not in ['ambertools', 'obabel', 'none']:
            global_log.error(f"Invalid protonation tool {protonation_tool}. Options: ambertools, obabel, none")
            return
    
    # Find the custom parameters if any
    parameter_sets = get_parameter_sets(custom_parameters, global_log)
    
    # Find the charges of each ligand if any
    ligand_charges = {}
    if charges is not None:
        for q in charges:
            ligand_name, ligand_charge = q.split(':')
            ligand_charges[ligand_name] = int(ligand_charge)
        
    # Process each ligand
    for ligand_info in selected_ligands:
        
        # Ligand specific properties and paths
        ligand_name = ligand_info['name']
        ligand_prop = conf.get_prop_dic(prefix=ligand_name)
        ligand_paths = conf.get_paths_dic(prefix=ligand_name)
        
        # STEP 1: Extract ligand from the PDB file
        ligand_paths["step1_ligand_extraction"]["input_structure_path"] = input_pdb
        ligand_prop["step1_ligand_extraction"]["heteroatoms"] = [ligand_info]
        
        print(ligand_prop["step1_ligand_extraction"]["heteroatoms"])
        extract_heteroatoms(**ligand_paths["step1_ligand_extraction"], properties=ligand_prop["step1_ligand_extraction"])
        
        # Distinguish between ligands with available parameter set (path A) or new ligand (path B)
        
        # Parameter set available 
        if ligand_name in parameter_sets:
            
            if format == 'gromacs':
                
                # If requested format is gromacs, produce the topology in amber format and convert it to gromacs
                global_log.info(f"Parameterizing {ligand_name} with custom parameter set")
                
                # STEP 2A: Use leap to generate the topology from the custom parameter set
                ligand_paths["step2A_leap_gen_top"]["input_frcmod_path"] = parameter_sets[ligand_name]['frcmod']
                ligand_paths["step2A_leap_gen_top"]["input_prep_path"] = parameter_sets[ligand_name]['prep']
                ligand_paths["step2A_leap_gen_top"]["forcefield"] = forcefields
                global_log.info("step2A_leap_gen_top: Generate topology with custom parameter set")
                leap_gen_top(**ligand_paths["step2A_leap_gen_top"], properties=ligand_prop["step2A_leap_gen_top"])
                            
                global_log.info("step3A_amber_to_gmx: Convert topology from AMBER to GROMACS")
                ligand_prop["step3A_amber_to_gmx"]["basename"] = ligand_name
                global_log.info("step3A_amber_to_gmx: Convert topology from AMBER to GROMACS")
                acpype_convert_amber_to_gmx(**ligand_paths["step3A_amber_to_gmx"], properties=ligand_prop["step3A_amber_to_gmx"])
                
                # Output files
                out_files = [ligand_paths["step3A_amber_to_gmx"]["output_path_top"], ligand_paths["step3A_amber_to_gmx"]["output_path_gro"]]
                
            elif format == 'amber':
                
                # If requested format is amber, use the custom parameter set directly
                global_log.info(f"Copying custom parameter set for {ligand_name}")
                out_files = [parameter_sets[ligand_name]['frcmod'], parameter_sets[ligand_name]['prep']]
                
        # Parameter set not available - generate topology from scratch
        else: 
            
            global_log.info(f"Parameterizing {ligand_name} using GAFF through acpype")
    
            # NOTE: the charge should be taken into account by protonation tool
            if protonation_tool == 'ambertools':
                
                # STEP 2B: Add hydrogens to the ligand using ambertools
                global_log.info("step2B_ambertools_reduce: Add hydrogens to the ligand")
                reduce_add_hydrogens(**ligand_paths["step2B_ambertools_reduce"], properties=ligand_prop["step2B_ambertools_reduce"])
            
            elif protonation_tool == 'obabel':
                
                # STEP 2B: Add hydrogens to the ligand using obabel
                global_log.info("step2B_obabel_reduce: Add hydrogens to the ligand")
                babel_add_hydrogens(**ligand_paths["step2B_obabel_reduce"], properties=ligand_prop["step2B_obabel_reduce"])
                
                # Modify minimize path
                ligand_paths["step3B_babel_minimize"]["input_path"] = ligand_paths["step2B_obabel_reduce"]["output_path"]
            
            elif protonation_tool == 'none':
                
                # STEP 2B: Skip ligand protonation
                global_log.info("Skipping step2B: No hydrogens added to the ligand. Make sure the ligand is correctly protonated.")
                
                # Modify minimize path
                ligand_paths["step3B_babel_minimize"]["input_path"] = ligand_paths["step1_ligand_extraction"]["output_heteroatom_path"]
        
        
            # STEP 3B: Minimize the ligand
            if not skip_min:
                global_log.info("step3B_babel_minimize: Energetically minimize the ligand")
                babel_minimize(**ligand_paths["step3B_babel_minimize"], properties=ligand_prop["step3B_babel_minimize"])
            else:
                global_log.info("Skipping step3B: No minimization performed. Make sure the ligand is correctly minimized.")
                ligand_paths["step4B_acpype_params_gmx"]["input_path"] = ligand_paths["step3B_babel_minimize"]["input_path"]
                ligand_paths["step4B_acpype_params_ac"]["input_path"] = ligand_paths["step3B_babel_minimize"]["input_path"]
            
            if ligand_name in ligand_charges:
                ligand_prop["step4B_acpype_params_gmx"]["charge"] = ligand_charges[ligand_name]
                ligand_prop["step4B_acpype_params_ac"]["charge"] = ligand_charges[ligand_name]
            else:
                ligand_prop["step4B_acpype_params_gmx"]["charge"] = None
                ligand_prop["step4B_acpype_params_ac"]["charge"] = None
                
            # Create gromacs topology
            if format == 'gromacs':
                ligand_prop["step4B_acpype_params_gmx"]["basename"] = ligand_name
                global_log.info("step4B_acpype_params_gmx: Generating GROMACS ligand parameters")
                acpype_params_gmx(**ligand_paths["step4B_acpype_params_gmx"], properties=ligand_prop["step4B_acpype_params_gmx"])
                out_files = [ligand_paths["step4B_acpype_params_gmx"]["output_path_itp"], ligand_paths["step4B_acpype_params_gmx"]["output_path_gro"]]
            
            # Create amber topology 
            if format == 'amber':
                ligand_prop["step4B_acpype_params_ac"]["basename"] = ligand_name
                global_log.info("step4B_acpype_params_ac: Generating AMBER ligand parameters")
                acpype_params_ac(**ligand_paths["step4B_acpype_params_ac"], properties=ligand_prop["step4B_acpype_params_ac"])
                out_files = [ligand_paths["step4B_acpype_params_ac"]["output_path_frcmod"], ligand_paths["step4B_acpype_params_ac"]["output_path_lib"]]
                
        # Copy the topology of this ligand to the final output folder
        copy_out_files(out_files, ligand_name, output_top_path)
        
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

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Ligand parameterization for AMBER or GROMACS")
    
    parser.add_argument('--config', dest='config_path', 
                        help='Configuration file (YAML)', 
                        required=True)

    parser.add_argument('--input_pdb', dest='input_pdb',
                        help='Path to the input PDB file with the ligands to parameterize.', 
                        required=True)
    
    parser.add_argument('--forcefields', dest='forcefields', nargs='+',
                        help='List of paths or file names of force fields to use in the parameterization. Only used within LEaP if custom parameters are also provided. Look in the $AMBERHOME/dat/leap/cmd/ directory for available force fields.',
                        required=False, default=['protein.ff14SB', 'DNA.bsc1', 'gaff'])
    
    parser.add_argument('--ligands', dest='ligand_names', nargs='+',
                        help='List of ligand names in the PDB file to parameterize. By default, all ligands are parameterized.', 
                        required=False, default=None)
    
    parser.add_argument('--charges', dest='charges', nargs='+',
                        help='List of charges for the ligands. Format: "ligand_name:charge". Ex: "JZ4:-2 FLP:1". Only used when no custom parameters are given and the ligand is parameterized using GAFF and antechamber through acpype. By default acpype will guess the charge based on the protonation state.',)
    
    parser.add_argument('--chains', dest='chains', nargs='+',
                        help='Chain IDs of the PDB to extract ligands from. Default: ["A"].',
                        required=False, default=['A'])

    parser.add_argument('--model', dest='model',
                        help='Model number of the PDB to extract ligands from. Default: 0.',
                        required=False, default=0)
    
    parser.add_argument('--format', dest='format', 
                        help='Format of the output topology files. Options: gromacs, amber. Default: gromacs.',
                        required=False, default='gromacs')
    
    parser.add_argument('--custom_parameters', dest='custom_parameters',
                        help='Path to folder with custom parameter sets for one or more ligands (.frcmod and .prep files with ligand name). LEaP will be used to generate the topology for these ligands. These ligands will rely on the charge given by the custom parameter set and the protonation state of its template.',
                        required=False, default=None)
    
    parser.add_argument('--protonation_tool', dest='protonation_tool', 
                        help='Protonation tool, use none to avoid ligand protonation. Options: ambertools, obabel, none. Default: ambertools.',
                        required=False, default='ambertools')
    
    parser.add_argument('--skip_min', action='store_true',
                        help='Skip the minimization step.',
                        required=False, default=False)
    
    parser.add_argument('--output_top_path', dest='output_top_path',
                        help='Output path for the folder with ligand topologies (Amber: .frcmod and .prep/.lib files, Gromacs: .gro and .itp files).',
                        required=False, default=None)

    parser.add_argument('--output', dest='output_path',
                        help="Output path. Default: working_dir_path in YAML config file.",
                        required=False, default=None)   

    args = parser.parse_args()
    
    main(configuration_path = args.config_path, input_pdb = args.input_pdb, forcefields=args.forcefields, ligand_names = args.ligand_names, charges = args.charges,  
         chains = args.chains, model = args.model, format = args.format, custom_parameters = args.custom_parameters, protonation_tool = args.protonation_tool, 
         skip_min = args.skip_min, output_top_path = args.output_top_path, output_path=args.output_path)
