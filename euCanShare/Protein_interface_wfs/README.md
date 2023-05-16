# Protein interface workflows

Protein interface biobb workflows created for euCanShare. They can be run together or as independent workflows. 

![alt text](../img/protein_protein_scheme.png?raw=true)

## Description

**1.PP_docking**: First uses pyDock to sample and score protein-protein poses starting from 2 pdb structures. Then performs a clustering based on the RMSD between poses and gives back the best ranking pose from each cluster (to obtain a set of good distinct poses). See corresponding folder for more details.

**2.Pose_refinement**: Takes in a zip file with docking poses and minimizes the energy of each of them using GROMACS. See corresponding folder for more details.

**3.Cavity_analysis**: Takes in a zip file with pdb structures and analyzes their cavities using Fpocket. It then filters those cavities according to a user-defined criteria. See corresponding folder for more details.

**4.Virtual_screening**: high-throughput virtual screening of selected pocket using library of ligands (SMILES) and Autodock. See corresponding folder for more details.

## General installation steps

Requirements: git, conda

1. Clone or download this repository:

```bash
git clone https://github.com/PabloNA97/biobb_workflows.git
```

2. Go to one of the workflow folders 

3. Use conda to create the corresponding environment. Install Anaconda first. If you are using an HPC environment load the miniconda module if available or install Miniconda.

```bash
conda env create -f environment.yml
```

4. Check the environment is listed and activate it (the environment name is inside environment.yml)

```bash
conda env list
conda activate <environment_name>
```

5. Provide necessary inputs and configuration through the YAML configuration file (see input.yml or input_HPC.yml) or through the command line arguments (see run_local.sh or run_HPC.sl). Then execute either in locally:

```bash
./run.sh
```

Or using an HPC environment:

```bash
sbatch run_HPC.sl
```

