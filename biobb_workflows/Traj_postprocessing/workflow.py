#!/usr/bin/env python3

from typing import List, Optional
from pathlib import Path
import numpy as np

import argparse
import time
import os

from Bio.PDB import PDBParser

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str

# Constants
# Possible ion names not recognized by GROMACS default "Ion" group
ions_library : List = ["K+", "CL-", "MG"]

# All other solvent names not recognized by GROMACS default "SOL" group
solvent_library : List = []

# Group names used for T-coupling and trajectory post-processing
solvent_group : str = "Solvent_group"
solute_group : str = "Solute_group"
output_group : str = "Output_group"


def get_residue_types(pdb_path: str, target_resnames: List[str]) -> List[str]:
    """Find residue names from target_resnames present in the PDB structure."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    target_set = set(target_resnames)
    found_residues = set()
    for residue in structure.get_residues():
        res_name = residue.get_resname().strip()
        if res_name in target_set:
            found_residues.add(res_name)
    return list(found_residues)


def get_atom_types(pdb_path: str, target_atom_names: List[str]) -> List[str]:
    """Find atom names from target_atom_names present in the PDB structure."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    target_set = set(target_atom_names)
    found_atoms = set()
    for atom in structure.get_atoms():
        atom_name = atom.get_name().strip()
        if atom_name in target_set:
            found_atoms.add(atom_name)
    return list(found_atoms)


def build_solvent_selection(solvent_names: List[str], ion_names: List[str]) -> str:
    """Build a GROMACS selection string for solvent and ions groups."""
    if solvent_names:
        solvent_selection = f'"SOL" | {" | ".join(f"r {s}" for s in solvent_names)}'
    else:
        solvent_selection = '"SOL"'
    if ion_names:
        ions_selection = f'"Ion" | {" | ".join(f"a {i}" for i in ion_names)}'
    else:
        ions_selection = '"Ion"'
    return f'{solvent_selection} | {ions_selection}'


def rename_last_ndx_group(ndx_path: str, new_name: str) -> None:
    """Rename the last group in a GROMACS index file."""
    with open(ndx_path, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines) - 1, -1, -1):
        line = lines[i].strip()
        if line.startswith("[") and line.endswith("]"):
            lines[i] = f"[ {new_name} ]\n"
            break
    with open(ndx_path, 'w') as f:
        f.writelines(lines)


def get_central_atom_index(pdb_path: str) -> int:
    """Return the 1-based index of the atom closest to the geometric center."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    atoms = list(structure.get_atoms())
    coords = np.array([atom.get_vector().get_array() for atom in atoms])
    center = coords.mean(axis=0)
    distances = np.linalg.norm(coords - center, axis=1)
    return int(np.argmin(distances)) + 1


def add_group(atom_indices: List, group_name: str, old_ndx_path: str, new_ndx_path: str) -> None:
    """Append a new group with given atom indices to an ndx file, saving to a new path."""
    COLUMNS = 15
    group_block = f"\n[ {group_name} ]\n"
    for i, idx in enumerate(atom_indices):
        group_block += f"{idx:>4}"
        if (i + 1) % COLUMNS == 0:
            group_block += "\n"
        else:
            group_block += " " if (i + 1) < len(atom_indices) else ""
    if len(atom_indices) % COLUMNS != 0:
        group_block += "\n"
    with open(old_ndx_path, 'r') as f:
        old_content = f.read()
    with open(new_ndx_path, 'w') as f:
        f.write(old_content)
        f.write(group_block)

def get_input_pdb(input_structure: str, gmx_bin: str, output_path: str, log) -> str:
    """     
    Extract a PDB structure from the input structure file (GRO or TPR) using editconf, and return its path.
    """
    
    # Create a directory for the pre-processing step
    step_dir = os.path.join(output_path, '0_input_pdb')
    os.makedirs(step_dir, exist_ok=True)
    
    # Extract PDB from the input structure (GRO or TPR)
    output_pdb = os.path.join(step_dir, 'input_structure.pdb')
    editconf(
        input_gro_path=input_structure,
        output_gro_path=output_pdb,
        properties={'binary_path': gmx_bin}
    )
    
    return output_pdb

def common_config_contents(
    gmx_bin: str = 'gmx',
    debug: bool = False,
    solvent_selection: str = '"SOL" | "Ion"',
    output_selection: str = f'"{solute_group}"',
    structure_path: str = '',
    input_topology_path: str = '',
    input_traj_path: str = ''
    
) -> str:
    return f"""
# Global properties (common for all steps)
global_properties:
  can_write_console_log: False
  restart: True
  working_dir_path: output
  remove_tmp: {not debug}

###################################################
# Section 1: Index group creation from structure  #
###################################################

# Add solvent and ions to index file
step1_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '{solvent_selection}'

# Add solute group (complement of solvent) to index file
step2_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    input_ndx_path: dependency/step1_make_ndx/output_ndx_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '! "{solvent_group}"'

# Add output group to index file
step3_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    input_ndx_path: dependency/step2_make_ndx/output_ndx_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '{output_selection}'

###################################################
# Section 2: Structure and trajectory processing  #
###################################################

# Center on solute  
step4_dry_str:
  tool: gmx_trjconv_str
  paths: 
    input_structure_path: {structure_path}
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_str_path: dry_structure.pdb
  properties:
    binary_path: {gmx_bin}     
    selection: "{solute_group}"
    center: False
    pbc: none
    
# Step 5 (Center group creation) is done in Python after extracting the 
# dry structure in step 4. The Center group is needed for centering 
# the trajectory in step 7.

# Extract dry (or full) trajectory
step6_dry_traj:
  tool: gmx_trjconv_trj
  paths:
    input_traj_path: {input_traj_path}
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: dry_traj.xtc
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
""" 

def complete_postprocessing(
    gmx_bin: str = 'gmx',
    input_topology_path: str = ''
) -> str:
    return f"""
# Make the molecules whole 
step7_whole:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step6_dry_traj/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: whole_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: False
    pbc: whole

# Cluster the molecules
step8_cluster:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7_whole/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: cluster_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    cluster_selection: "{solute_group}"
    center: False
    pbc: cluster

# Extract the first frame to use as ref
step9_extract_ref:
  tool: gmx_trjconv_trj
  paths:
    input_traj_path: dependency/step8_cluster/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: first_frame.gro
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
    dump: 0

# Use as ref with pbc nojump 
step10_nojump:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step8_cluster/output_traj_path
    input_top_path: dependency/step9_extract_ref/output_traj_path
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: nojump_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: False
    pbc: nojump

# Center the system - whole solute group / central atom
step11_center:
  tool: gmx_image
  paths: 
    input_traj_path: dependency/step10_nojump/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: centered_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center_selection: "Center"
    center: True
    ur: compact
    pbc: none

# Image the trajectory to put all molecules back in the box
step12_image:
  tool: gmx_image
  paths: 
    input_traj_path: dependency/step11_center/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: centered_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: False
    ur: compact
    pbc: mol

# Fit the trajectory by rotation and translation
step13_fit:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step12_image/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_bin}
    fit_selection: "{solute_group}"
    output_selection: "{output_group}"
    center: False
    fit: rot+trans
"""

def fast_postprocessing(
    gmx_bin: str = 'gmx',
    input_topology_path: str = ''
) -> str:
    return f"""
# Center the trajectory on central atom
step7_center:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step6_dry_traj/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: center_traj.xtc
  properties:
    binary_path: {gmx_bin}
    center_selection: "Center"
    output_selection: "{output_group}"
    center: True
    ur: compact
    pbc: none

# Apply periodic boundary conditions imaging
step8_image:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7_center/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: imaged_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: False
    ur: compact
    pbc: mol

# Fit the trajectory by rotation and translation
step9_fit:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step8_image/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_bin}
    fit_selection: "{solute_group}"
    output_selection: "{output_group}"
    center: False
    fit: rot+trans
"""

def config_contents(
    gmx_bin: str = 'gmx',
    debug: bool = False,
    solvent_selection: str = '"SOL" | "Ion"',
    output_selection: str = f'"{solute_group}"',
    structure_path: str = '',
    input_topology_path: str = '',
    input_traj_path: str = '',
    fast: bool = False
) -> str:
    common_contents = common_config_contents(gmx_bin, 
                                            debug, 
                                            solvent_selection, 
                                            output_selection, 
                                            structure_path, 
                                            input_topology_path,
                                            input_traj_path)
    
    if fast:
        return common_contents + fast_postprocessing(gmx_bin, input_topology_path)
    else:
        return common_contents + complete_postprocessing(gmx_bin, input_topology_path)


def create_config_file(output_path: str, **config_args) -> str:
    """Write YAML config to output_path/config.yml and return its path."""
    config_path = os.path.join(output_path, 'config.yml')
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))
    return config_path


def main_wf(
    input_traj_path: str,
    input_topology_path: str,
    input_structure_path: str,
    gmx_bin: str = 'gmx',
    keep_solvent: bool = False,
    residues_to_keep: Optional[List[int]] = None,
    extra_ions: List[str] = [],
    extra_solvents: List[str] = [],
    fast: bool = False,
    debug: bool = False,
    output_path: str = 'output',
    output_traj_path='trajectory.xtc',
    output_str_path='structure.pdb'
):
    '''
    Post-process a GROMACS MD trajectory: strip solvent, center, image, and fit.

    Inputs
    ------
        input_traj_path:
            path to the input trajectory file (.xtc). Required.
        input_topology_path:
            path to the binary run input file (.tpr). Required.
        input_structure_path:
            path to an input structure file (.gro or .pdb).
            Used to define solvent/output index groups and to find the center group for centering.
            Make sure the structure is not broken due to PBC.
        gmx_bin:
            GROMACS binary path. Default: gmx
        keep_solvent:
            include solvent and ions in the output structure and trajectory.
            Default: False (dry output)
        residues_to_keep:
            residue indices to retain in the output besides the solute.
            Default: None (only solute)
        extra_ions:
            additional ion atom names to include in the solvent group (e.g. --ions NA+ CA2+). Default: []
        extra_solvents:
            additional solvent residue names to include in the solvent group (e.g. --solvents TIP3 TIP4). Default: []
        debug:
            keep intermediate files. Default: False
        output_path:
            directory for workflow output. Default: output

    Outputs
    -------
        output/
            dry_structure.pdb   — processed structure
            fitted_traj.xtc     — processed trajectory (dry, centered, imaged, fitted)
        global_paths (dict), global_prop (dict)
    '''

    start_time = time.time()

    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=True)

    # Initialize a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)

    ###########################################
    # Build GROMACS selection strings for ndx #
    ###########################################

    # Convert provided GRO to PDB for solvent/ion detection if needed
    if Path(input_structure_path).suffix.lower() != '.pdb':
        pdb_structure_path = get_input_pdb(input_structure_path, gmx_bin, output_path, global_log)
    else:
        pdb_structure_path = input_structure_path

    # Construct solvent selection, solute selection will be the the rest of the system
    solvent_library.extend(extra_solvents)
    ions_library.extend(extra_ions)
    solvent_names = get_residue_types(pdb_structure_path, solvent_library)
    ion_names = get_atom_types(pdb_structure_path, ions_library)
    solvent_selection = build_solvent_selection(solvent_names, ion_names)

    # Determine output selection based on user options
    if keep_solvent:
        output_selection = '"System"'
    elif residues_to_keep:
        residues_sel = f"ri {' '.join(str(r) for r in residues_to_keep)}"
        output_selection = f'"{solute_group}" | {residues_sel}'
    else:
        output_selection = f'"{solute_group}"'

    #################################
    # Create workflow configuration #
    #################################

    config_path = create_config_file(
        output_path,
        gmx_bin=gmx_bin,
        debug=debug,
        solvent_selection=solvent_selection,
        output_selection=output_selection,
        structure_path=input_structure_path,
        input_topology_path=input_topology_path,
        input_traj_path=input_traj_path,
        fast=fast
    )
    global_log.info(f"Configuration file: {config_path}")

    conf = settings.ConfReader(config=config_path, system=None)
    prop = conf.get_prop_dic()
    paths = conf.get_paths_dic()

    #######################
    # Create index files  #
    #######################

    global_log.info("step1_make_ndx: Create solvent group in index file")
    make_ndx(**paths['step1_make_ndx'], properties=prop['step1_make_ndx'])
    rename_last_ndx_group(paths['step1_make_ndx']['output_ndx_path'], solvent_group)

    global_log.info("step2_make_ndx: Create solute group in index file")
    make_ndx(**paths['step2_make_ndx'], properties=prop['step2_make_ndx'])
    rename_last_ndx_group(paths['step2_make_ndx']['output_ndx_path'], solute_group)

    global_log.info("step3_make_ndx: Create output group in index file")
    make_ndx(**paths['step3_make_ndx'], properties=prop['step3_make_ndx'])
    rename_last_ndx_group(paths['step3_make_ndx']['output_ndx_path'], output_group)

    input_ndx_path = paths['step3_make_ndx']['output_ndx_path']

    ########################################
    # Extract solute and find central atom #
    ########################################
    
    centering_step = 'step7_center' if fast else 'step11_center'

    # Extract the solute from the input structure
    try:
        global_log.info("step4_dry_str: Center dry structure")
        gmx_trjconv_str(**paths['step4_dry_str'], properties=prop['step4_dry_str'])
        dry_str_ok = True
    except SystemExit as e:
        global_log.error(f"Structure post-processing failed (SystemExit, code={e.code})")
        dry_str_ok = False
    except Exception:
        global_log.exception("Structure post-processing failed with unexpected exception")
        dry_str_ok = False

    # Extract central atom index from solute to center the trajectory
    if dry_str_ok:
        os.mkdir(os.path.join(output_path, 'step5_center_group'))
        center_ndx_path = os.path.join(output_path, 'step5_center_group', 'center.ndx')
        central_index = get_central_atom_index(paths['step4_dry_str']['output_str_path'])
        global_log.info(f"Central atom index for centering: {central_index}")
        add_group([central_index], 'Center', input_ndx_path, center_ndx_path)
        paths[centering_step]['input_index_path'] = center_ndx_path
    else:
        # Fall back to centering on the whole solute group
        global_log.warning("Structure post-processing failed falling back to Solute_group for centering")
        prop[centering_step]['center_selection'] = solute_group

    #########################################################################
    # Process trajectory: strip, whole, nojump, cluster, center, image, fit #
    #########################################################################
        
    final_step = 'step9_fit' if fast else 'step13_fit'

    try:
        global_log.info("step6_dry_traj: Extract dry trajectory")
        gmx_trjconv_trj(**paths['step6_dry_traj'], properties=prop['step6_dry_traj'])
        
        if fast:
            global_log.info("Running in fast mode: skipping whole, nojump, and cluster steps")
            
            global_log.info("step7_center: Center the trajectory using the Center group")
            gmx_image(**paths['step7_center'], properties=prop['step7_center'])
            if not debug:
                os.remove(paths['step6_dry_traj']['output_traj_path'])
            
            global_log.info("step8_image: Image the trajectory to put all molecules back in the box")
            gmx_image(**paths['step8_image'], properties=prop['step8_image'])
            if not debug:
                os.remove(paths['step7_center']['output_traj_path'])
                
            global_log.info("step9_fit: Fit the trajectory by rotation and translation")
            gmx_image(**paths['step9_fit'], properties=prop['step9_fit'])
        else:
            global_log.info("Running in complete mode: performing whole, nojump, and cluster steps for better results")
            
            global_log.info("step7_whole: Make the molecules whole in the trajectory")
            gmx_image(**paths['step7_whole'], properties=prop['step7_whole'])
            if not debug:
                os.remove(paths['step6_dry_traj']['output_traj_path'])

            global_log.info("step8_cluster: Cluster the molecules in the trajectory")
            gmx_image(**paths['step8_cluster'], properties=prop['step8_cluster'])
            if not debug:
                os.remove(paths['step7_whole']['output_traj_path'])

            global_log.info("step9_extract_ref: Extract the first frame to use as reference")
            gmx_trjconv_trj(**paths['step9_extract_ref'], properties=prop['step9_extract_ref'])

            global_log.info("step10_nojump: Make the trajectory whole with nojump PBC")
            gmx_image(**paths['step10_nojump'], properties=prop['step10_nojump'])
            if not debug:
                os.remove(paths['step8_cluster']['output_traj_path'])

            global_log.info("step11_center: Center the trajectory using the Center group")
            gmx_image(**paths['step11_center'], properties=prop['step11_center'])
            if not debug: 
                os.remove(paths['step10_nojump']['output_traj_path'])
            
            global_log.info("step12_image: Image the trajectory to put all molecules back in the box")
            gmx_image(**paths['step12_image'], properties=prop['step12_image'])
            if not debug:
                os.remove(paths['step11_center']['output_traj_path'])
            
            global_log.info("step13_fit: Fit the trajectory by rotation and translation")
            gmx_image(**paths['step13_fit'], properties=prop['step13_fit'])

    except SystemExit as e:
        global_log.error(f"Trajectory post-processing failed (SystemExit, code={e.code})")
    except Exception:
        global_log.exception("Trajectory post-processing failed with unexpected exception")
        
    # Move final outputs to user-specified paths
    final_traj_path = paths[final_step]['output_traj_path']
    final_str_path = paths['step4_dry_str']['output_str_path']
    os.rename(final_traj_path, output_traj_path)
    os.rename(final_str_path, output_str_path)

    elapsed = time.time() - start_time
    global_log.info('')
    global_log.info('Execution successful:')
    global_log.info(f'  Workflow path: {output_path}')
    global_log.info(f'  Config file:   {config_path}')
    global_log.info(f'  Elapsed time:  {elapsed/60:.1f} minutes')

    return paths, prop


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Post-process a GROMACS MD trajectory: strip solvent, center, image, fit."
    )

    ###############
    # Input files #
    ###############

    parser.add_argument(
        '--input_traj', dest='input_traj_path', type=str, required=True,
        help="Input trajectory file (.xtc). Required."
    )
    parser.add_argument(
        '--input_top', dest='input_topology_path', type=str, required=True,
        help="Input binary run input file (.tpr). Required."
    )
    parser.add_argument(
        '--input_structure', dest='input_structure_path', type=str, required=True,
        help=("Input structure file (.gro or .pdb). Used to define solvent/output "
              "index groups and to find the center group for centering. "
              "Make sure the structure is not broken due to PBC.")
    )

    #########################
    # Configuration options #
    #########################

    parser.add_argument(
        '--gmx_bin', type=str, required=False, default='gmx',
        help="GROMACS binary path. Default: gmx"
    )
    parser.add_argument(
        '--keep_solvent', action='store_true', required=False, default=False,
        help="Keep solvent and ions in the output structure and trajectory. Default: False"
    )
    parser.add_argument(
        '--keep_residues', dest='residues_to_keep', type=int, nargs='+',
        required=False,
        help=("Residue indices to keep in the output besides the solute "
              "(e.g. --keep_residues 15 23 105). Default: None")
    )
    parser.add_argument(
        '--ions', dest='extra_ions', type=str, nargs='+',
        required=False, default=[],
        help=("Additional ion atom names to include in the solvent group (e.g. --ions NA+ CA2+). Default: []" )
    )
    parser.add_argument(
        '--solvents', dest='extra_solvents', type=str, nargs='+',
        required=False, default=[],
        help=("Additional solvent residue names to include in the solvent group (e.g. --solvents TIP3 TIP4). Default: []" )
    )
    parser.add_argument(
        '--fast', action='store_true', required=False, default=False,
        help="Skip making solute whole, removing jumps and clustering. Default: False"
    )
    parser.add_argument(
        '--debug', action='store_true', required=False, default=False,
        help="Keep intermediate files. Default: False"
    )
    parser.add_argument(
        '--output', dest='output_path', type=str, required=False, default='output',
        help="Output directory path for the workflow where the steps will be written. Default: output"
    )
    parser.add_argument(
        '--output_traj', dest='output_traj_path', type=str, required=False, default='trajectory.xtc',
        help="Output trajectory file name (e.g. processed_traj.xtc). Default: trajectory.xtc"
    )
    parser.add_argument(
        '--output_str', dest='output_str_path', type=str, required=False, default='structure.pdb',
        help="Output structure file name (e.g. processed_structure.pdb). Default: structure.pdb"
    )

    args = parser.parse_args()

    main_wf(
        input_traj_path=args.input_traj_path,
        input_topology_path=args.input_topology_path,
        input_structure_path=args.input_structure_path,
        gmx_bin=args.gmx_bin,
        keep_solvent=args.keep_solvent,
        residues_to_keep=args.residues_to_keep,
        extra_ions=args.extra_ions,
        extra_solvents=args.extra_solvents,
        fast = args.fast,
        debug=args.debug,
        output_path=args.output_path,
        output_traj_path=args.output_traj_path,
        output_str_path=args.output_str_path
    )
