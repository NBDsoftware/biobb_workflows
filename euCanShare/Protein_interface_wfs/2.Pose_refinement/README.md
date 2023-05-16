# Pose refinement

![alt text](../../img/pose_refinement.png?raw=true)

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate protein_protein_wf2
```

To install it in an HPC environment, do the same after loading the corresponding Conda or Miniconda module.

See options for worklow:

```bash
vi input_HPC.yml
python biobb_pose_refinement.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run in an HPC environment adapt the run_HPC.sl and input_HPC.yml scripts and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

The output will be generated in the "working_dir_path" folder selected in the corresponding YAML input. The global log files will be in "/working_dir_path/log.out" and "/working_dir_path/log.err". Each step will have its own log files and output in a separate folder inside "/working_dir_path".

## Description

This workflow has several steps. The input for the workflow is a zip file containing pdb structures. Its path can be provided through the YAML configuration file or through the command line arguments. The YAML configuration file will contain the default settings and paths of the workflow. The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

- **Step 0**: extraction of all the poses from the input zip file. The next steps are done for each of the poses.

**Steps 1-5**: prepare the pdb file.

- **Step 1**: add missing side chain atoms using biobb_structure_checking (and Modeller suite if a license key is provided)

- **Step 2**: Add missing disulfide bonds. Use carefully. This step is executed only if "--fix_ss" is used in the command line arguments.

- **Step 3**: flip the clashing amide groups to avoid clashes.

- **Step 4**: fix stereochemical errors in residue side-chains changing their chirality.

- **Step 5**: renumber atomic indices.

- **Step 6**: uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the prepared structure file in step 2. Hydrogen atoms will be added in this step. A force field and water model are chosen with the force_field and water_type properties.

- **Step 7**: generate a simulation box around the structure.

- **Step 8**: generate a box of solvent around the structure.

- **Step 9-10**: randomly replace solvent molecules with monoatomic ions

- **Steps 11 - 13**: energy minimization with steepest descend and constraints on the backbone

- **Steps 14 - 16**: energy minimization with steepest descend and without constraints on the backbone

- **Step 17**: convert from .gro format to .pdb

- **Step 18**: copy final pdb to merged results folder and zip all the energy-minimized structures.


