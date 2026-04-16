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
ions_library = ["K+", "CL-", "MG"]

# All other solvent names not recognized by GROMACS default "SOL" group
solvent_library = []

# Group names used for T-coupling and trajectory post-processing
solvent_group = "Solvent_group"
solute_group = "Solute_group"
output_group = "Output_group"


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
    
    
def config_contents(
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

# Extract dry (or full) structure as PDB
step4_dry_str:
  tool: gmx_trjconv_str
  paths:
    input_structure_path: {structure_path}
    input_top_path: {input_topology_path}
    input_index_path: dependency/step3_make_ndx/output_ndx_path
    output_str_path: dry_structure.pdb
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
    center: True
    pbc: mol
    ur: compact

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


def create_config_file(output_path: str, **config_args) -> str:
    """Write YAML config to output_path/config.yml and return its path."""
    config_path = os.path.join(output_path, 'config.yml')
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))
    return config_path


def main_wf(
    input_traj_path: str,
    input_topology_path: str,
    input_structure_path: Optional[str] = None,
    gmx_bin: str = 'gmx',
    keep_solvent: bool = False,
    residues_to_keep: Optional[List[int]] = None,
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
            path to an input structure file (.gro or .pdb). Optional. If not
            provided, coordinates are extracted from the TPR.
        gmx_bin:
            GROMACS binary path. Default: gmx
        keep_solvent:
            include solvent and ions in the output structure and trajectory.
            Default: False (dry output)
        residues_to_keep:
            residue indices to retain in the output besides the solute.
            Default: None (only solute)
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

    ####################################
    # Resolve / extract structure file #
    ####################################

    if input_structure_path is None:
        # NOTE: then extract PDB from TPR for solvent/ion detection and dry structure generation
        pdb_structure_path = get_input_pdb(input_topology_path, gmx_bin, output_path, global_log)
    elif Path(input_structure_path).suffix.lower() != '.pdb':
        # NOTE: convert provided GRO to PDB for solvent/ion detection and dry structure generation
        pdb_structure_path = get_input_pdb(input_structure_path, gmx_bin, output_path, global_log)
    else:
        pdb_structure_path = input_structure_path

    ##########################################
    # Build GROMACS selection strings for ndx #
    ##########################################

    solvent_names = get_residue_types(pdb_structure_path, solvent_library)
    ion_names = get_atom_types(pdb_structure_path, ions_library)
    solvent_selection = build_solvent_selection(solvent_names, ion_names)

    if keep_solvent:
        output_selection = '"System"'
    elif residues_to_keep:
        residues_sel = f"ri {' '.join(str(r) for r in residues_to_keep)}"
        output_selection = f'"{solute_group}" | {residues_sel}'
    else:
        output_selection = f'"{solute_group}"'

    ##########################
    # Create and load config #
    ##########################

    config_path = create_config_file(
        output_path,
        gmx_bin=gmx_bin,
        debug=debug,
        solvent_selection=solvent_selection,
        output_selection=output_selection,
        structure_path=pdb_structure_path,
        input_topology_path=input_topology_path,
        input_traj_path=input_traj_path
    )
    global_log.info(f"Configuration file: {config_path}")

    conf = settings.ConfReader(config=config_path, system=None)
    prop = conf.get_prop_dic()
    paths = conf.get_paths_dic()

    ###########################
    # Step 1-3: NDX creation  #
    ###########################

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

    ###################################
    # Step 4: Dry structure (PDB out) #
    ###################################

    try:
        global_log.info("step4_dry_str: Extract dry structure")
        gmx_trjconv_str(**paths['step4_dry_str'], properties=prop['step4_dry_str'])
        global_log.info("step4_dry_str: Completed successfully")
        dry_str_ok = True
    except SystemExit as e:
        global_log.error(f"step4_dry_str failed (SystemExit, code={e.code})")
        dry_str_ok = False
    except Exception:
        global_log.exception("step4_dry_str failed with unexpected exception")
        dry_str_ok = False

    # Build Center group index for centering step
    if dry_str_ok:
        os.mkdir(os.path.join(output_path, 'step5_center_group'))
        center_ndx_path = os.path.join(output_path, 'step5_center_group', 'center.ndx')
        central_index = get_central_atom_index(paths['step4_dry_str']['output_str_path'])
        add_group([central_index], 'Center', input_ndx_path, center_ndx_path)
        paths['step7_center']['input_index_path'] = center_ndx_path
    else:
        # Fall back to centering on the whole solute group
        global_log.warning("step4_dry_str failed: falling back to Solute_group for centering")
        prop['step7_center']['center_selection'] = solute_group

    ##############################
    # Steps 5-8: Trajectory work #
    ##############################

    try:
        global_log.info("step6_dry_traj: Extract dry trajectory")
        gmx_trjconv_trj(**paths['step6_dry_traj'], properties=prop['step6_dry_traj'])

        global_log.info("step7_center: Center the trajectory")
        gmx_image(**paths['step7_center'], properties=prop['step7_center'])
        if not debug:
            os.remove(paths['step6_dry_traj']['output_traj_path'])

        global_log.info("step8_image: Image the trajectory")
        gmx_image(**paths['step8_image'], properties=prop['step8_image'])
        if not debug:
            os.remove(paths['step7_center']['output_traj_path'])

        global_log.info("step9_fit: Fit the trajectory")
        gmx_image(**paths['step9_fit'], properties=prop['step9_fit'])
        if not debug:
            os.remove(paths['step8_image']['output_traj_path'])

    except SystemExit as e:
        global_log.error(f"Trajectory post-processing failed (SystemExit, code={e.code})")
    except Exception:
        global_log.exception("Trajectory post-processing failed with unexpected exception")
        
    # Move final outputs to user-specified paths
    final_traj_path = paths['step9_fit']['output_traj_path']
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
        '--input_structure', dest='input_structure_path', type=str, required=False,
        help=("Input structure file (.gro or .pdb). Used to define solvent/output "
              "index groups and to generate the dry structure PDB. "
              "If not provided, coordinates are extracted from the TPR. Default: None")
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
        debug=args.debug,
        output_path=args.output_path,
        output_traj_path=args.output_traj_path,
        output_str_path=args.output_str_path
    )
