# Virtual Screening

![alt text](../../img/virtual_screening.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
conda env create -f environment.yml
conda activate single_protein_wf3
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

To run in an HPC environment adapt the run_HPC.sl script and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

By default, the output will be generated in the "working_dir_path" folder selected in the YAML configuration file. However, the "--output" command line option will overwrite "working_dir_path". The global log files will be in "output/log.out" and "output/log.err". Each step will have its own log files and output in a separate folder inside the output folder.

## Inputs

### Configuration file

Take a look at the YAML configuration file to see the different properties that can be set.

```bash
vi input.yml
```

Specially important are: the selection of the cavity or residue defining the pocket, the box size around the residue selection or cavity, the number of top ligands to show in the final summary and the settings for Autodock Vina (cpus and exhaustiveness). Make sure the binary path specified and the module loaded in the run file agree between them.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python biobb_docking_htvs.py --help
```

Specially important are: the configuration file path, the path to the ligand library (SMILES format), the path to the target structure (pdb format), the output path, the number of top ligands to keep in the final summary, the flag to keep the poses of the best ligands and the flag to dock to a box created around a residue selection instead of a pocket found with Fpocket (a zip file containing pockets from an Fpocket analysis). 

## Description

This workflow has several steps. The input for the workflow is a ligand library in SMILES format, a target structure in pdb format and either a pocket from an Fpocket analysis or a selection of residues. The workflow will dock the ligands to the target structure and rank them according to their affinity.

- **Step 1**: selection of cavity that will be used to dock the ligands. Autodock will use a box created surrounding either: a pocket from an input zip file (see cavity analysis workflow) or a selection of residues. Make sure the selected pocket or residues exist in the input files. By default the pocket is used to create the box. To use a residue selection instead add the --dock_to_residues flag in the command line call to the workflow.

- **Step 2**: creation of box surrounding the selected cavity or residues.

- **Step 3**: addition of H atoms (generation of .pdbqt file from .pdb). The charges in the generated pdbqt file are ignored by Autodock Vina. The H atoms though are relevant. To avoid adding H select mode = null in the control file.

- **Step 4**: convert SMILES to .pdbqt format. Here the 3D structure is generated from the SMILES using OpenBabel. The current approach makes use of --gen3D to generate the lowest energy conformer ([see manual](https://open-babel.readthedocs.io/en/latest/3DStructureGen/SingleConformer.html#gen3d)). In the future faster options will be included (still some optimization might be needed as AutoDock needs a reasonably good initial conformation for the ligand).

- **Step 5**: docking of the ligand onto the target structure using [AutoDock Vina](https://vina.scripps.edu/manual/#summary). The target is kept rigid but different ligand conformers are explored.

- **Step 6**: find top scoring ligands, save their rank, names, smiles and affinity in a csv file.

- **Step 7**: save poses of top scoring ligands (only if --keep_poses is used)


