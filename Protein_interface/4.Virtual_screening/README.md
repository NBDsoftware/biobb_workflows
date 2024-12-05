# Virtual Screening

![alt text](../../img/virtual_screening2.png?raw=true)

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate protein_protein_wf4
```

See options for worklow:

```bash
vi input.yml
python workflow.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

To run in an HPC environment adapt the run_HPC.sl script and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

The output will be generated in the "working_dir_path" folder selected in the corresponding YAML input. The global log files will be in "/working_dir_path/log.out" and "/working_dir_path/log.err". Each step will have its own log files and output in a separate folder inside "/working_dir_path".

## Description

This workflow has several steps. The input for the workflow is a zip file containing pockets (from an Fpocket analysis), a pdb structure of the target and a library of ligands in SMILES format. Paths can be provided through the YAML configuration file or through the command line arguments. The YAML configuration file will contain the default settings and paths of the workflow. The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

- **Step 1**: pocket selection from input zip file with all pockets. Make sure the selected pocket exists in the input file.

- **Step 2**: creation of box surrounding the selected cavity.

- **Step 3**: addition of H atoms and charges to the target structure (generation of .pdbqt file from .pdb).

- **Step 4**: convert SMILES to .pdbqt format. Here the 3D structure is generated from the SMILES using OpenBabel. The current approach makes use of --gen3D to generate the lowest energy conformer ([see manual](https://open-babel.readthedocs.io/en/latest/3DStructureGen/SingleConformer.html#gen3d)). In the future faster options will be included (still some optimization might be needed as AutoDock needs a reasonably good initial conformation for the ligand).

- **Step 5**: docking of the ligand onto the target structure using [AutoDock Vina](https://vina.scripps.edu/manual/#summary). The target is kept rigid but different ligand conformers are explored.

- **Step 6**: find top scoring ligands, save their rank, names, smiles and affinity in a csv file.

- **Step 7**: save poses of top scoring ligands (only if --keep_poses is used)


