# Virtual Screening

![alt text](../../img/virtual_screening.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
conda env create -f environment.yml
conda activate biobb_sp_virtual_screening
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the input.yml.

To run in an HPC environment adapt the run_HPC.sl script and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

By default, the output will be generated in the `working_dir_path` folder selected in the YAML configuration file. However, the `--output` command line option will overwrite `working_dir_path`. Global log files will be generated in the output folder. Each step will have its own log files and output in a separate subfolder.

## Inputs

### Configuration file

Take a look at the YAML configuration file to see the different properties that can be set.

```bash
vi input.yml
```

Specially important are: the selection of the cavity or residue defining the pocket, the box size around the residue selection or cavity and the settings for Autodock Vina (cpus and exhaustiveness). Make sure the binary path specified and the module loaded in the run file agree between them. Note that, selecting too many residues to define the cavity or choosing a large box size can increase the computational cost of the docking.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python workflow.py --help
```

Specially important are: the configuration file path, the path to the ligand library, the path to the target structure, the output path, the number of top ligands to keep in the final summary, the flag to keep the poses of the best ligands and the flag to dock to a box created around a residue selection instead of a pocket found with Fpocket (a zip file containing pockets from an Fpocket analysis). See al options here:

```
usage: workflow.py [-h] -c CONFIG_PATH -lib LIGAND_LIB [-s STRUCTURE_PATH] [-pz INPUT_POCKETS_ZIP] [-p POCKET] [-o OUTPUT_PATH] [-nl NUM_TOP_LIGANDS] [-kp] [-dr] [-cpus CPUS] [-ex EXHAUSTIVENESS] [-d]

Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks

options:
  -h, --help            show this help message and exit
  -c CONFIG_PATH, --config CONFIG_PATH
                        Configuration file (YAML)
  -lib LIGAND_LIB, --ligand_lib LIGAND_LIB
                        Path to file with the ligand library. The format should be SMILES (.smi) or SDF (.sdf). For .smi files, one ligand per line is expected: 'smiles name'. For sdf files, the file may contain one
                        or more ligands.
  -s STRUCTURE_PATH, --structure_path STRUCTURE_PATH
                        Path to file with target structure (PDB format)
  -pz INPUT_POCKETS_ZIP, --input_pockets_zip INPUT_POCKETS_ZIP
                        Path to file with pockets in a zip file
  -p POCKET, --pocket POCKET
                        Pocket number to be used (default: 1)
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path (default: working_dir_path in YAML config file)
  -nl NUM_TOP_LIGANDS, --num_top_ligands NUM_TOP_LIGANDS
                        Number of top ligands to be saved (default: all successfully docked ligands)
  -kp, --keep_poses     Save docking poses for top ligands (default: False)
  -dr, --dock_to_residues
                        Dock to residues instead of pocket. Define the docking box using a set of residues instead of a pocket. (default: False)
  -cpus CPUS, --cpus CPUS
                        Number of CPUs to use for each docking (default: 1)
  -ex EXHAUSTIVENESS, --exhaustiveness EXHAUSTIVENESS
                        Exhaustiveness of the docking (default: 8)
  -d, --debug           Keep intermediate files for debugging (default: False)
```

## Description

This workflow has several steps. The input for the workflow is a ligand library in SMILES or SDF format, a target structure in pdb format and either a pocket from an Fpocket analysis or a selection of residues. The workflow will dock the ligands to the target structure and rank them according to their affinity. The top scoring ligands will be saved in a csv file with their rank, name, identifier and affinity. The poses of the top scoring ligands can also be saved.

- **Step 1**: selection of cavity that will be used to dock the ligands. Autodock will use a box created surrounding either: a pocket from an input zip file (see cavity analysis workflow) or a selection of residues. Make sure the selected pocket or residues exist in the input files. By default the pocket is used to create the box. To use a residue selection instead add the `--dock_to_residues` flag in the command line call to the workflow.

- **Step 2**: creation of box surrounding the selected cavity or residues.

- **Step 3**: addition of H atoms to the receptor (generation of .pdbqt file from .pdb). The charges in the generated pdbqt file are ignored by Autodock Vina. However the correct protonation of the receptor is important (see [Vina FAQs](https://autodock-vina.readthedocs.io/en/latest/faq.html)). To avoid adding H to the receptor select `mode = null` in the control file.

- **Step 4**: prepare ligands for docking. If the ligand library is in **SMILES format** (.smi) the workflow will use OpenBabel to convert 2D SMILES to protonated 3D ligand conformers in .pdbqt format. If the ligand library is in **SDF format** (.sdf) the workflow will use OpenBabel to convert the ligands to .pdbqt format without changing the protonation state or the conformer. Note that the protonation state of the ligand is important to obtain the correct affinity with Autodock Vina (see [Vina FAQs](https://autodock-vina.readthedocs.io/en/latest/faq.html).

    - To convert 2D SMILES to 3D conformers, the workflow uses the following obabel command: `obabel <input> -O <output> --gen3d -p <ph_value> -xh`. It generates the [lowest energy conformer](https://openbabel.org/docs/3DStructureGen/SingleConformer.html#gen3d) with hydrogen atoms for a given pH ([see manual](https://openbabel.org/docs/Command-line_tools/babel.html)). The protonation algorithm uses an empirical model to determine the pKa of the different atoms ([see transforms in phmodel](https://github.com/openbabel/openbabel/blob/master/data/phmodel.txt)). Thus the results will probably not be accurate if the microenvironment of the ligand is important for the protonation state. The pH value can be set in the control file. The default value is 7.4. For production screening it is recommended to use LigPrep or [gypsum-dl](https://github.com/durrantlab/gypsum_dl) instead.

    - To convert from SDF format to PDBQT format, the workflow uses the following obabel command: `obabel <input> -O <output> -xh`. 

- **Step 5**: docking of the ligand onto the target structure using [AutoDock Vina](https://vina.scripps.edu/manual/#summary). The target is kept rigid but some degrees of freedom of the ligand conformation are explored (rings are rigid for example). 

- **Step 6**: save poses of top scoring ligands (only if `--keep_poses` is used)


