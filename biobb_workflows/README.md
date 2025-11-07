# BioBB Workflows

BioBB workflows are ready-to-use pipelines built using BioExcel Building Blocks (BioBB) to perform common tasks in biomolecular simulations and modeling. Each workflow is contained in its own folder with all the necessary files to run it (input files, configuration files, execution scripts, environment definition, etc).

Different workflows can be combined to create more complex pipelines.

## Single protein pipeline

![alt text](../img/Single_protein_scheme.png?raw=true)

## Protein-protein docking pipeline

![alt text](../img/protein_protein_scheme.png?raw=true)

## Description

Each workflow is defined on a separate folder:

**Protein_preparation**: prepares a protein structure from a PDB file for further simulations or modeling. It includes tasks such as adding missing residues/atoms and selecting the protonation state of residues at a given pH.

**Ligand_parameterization**: parameterizes a small molecule starting from its SMILES representation using Antechamber and GAFF force field. See corresponding folder for more details.

**MD_gromacs**: prepares and launches a GROMACS MD simulation starting from a prepared PDB structure file. See corresponding folder for more details.

**Cavity_analysis**: clustering of MD trajectory and extraction of most representative structures. Followed by a cavity analysis using Fpocket on those structures. See corresponding folder for more details.

**VS_autodock**: high-throughput virtual screening of selected pocket using library of ligands (SMILES) and AutoDock Vina. See corresponding folder for more details.

**PP_docking**: First uses pyDock to sample and score protein-protein poses starting from 2 pdb structures. Then performs a clustering based on the RMSD between poses and gives back the best ranking pose from each cluster (to obtain a set of good distinct poses). See corresponding folder for more details.

**Pose_refinement**: Takes in a zip file with docking poses and minimizes the energy of each of them using GROMACS. See corresponding folder for more details.

## How to install

Requirements: git, conda

Clone or download this repository:

```bash
git clone https://github.com/NBDsoftware/biobb_workflows.git
```

Go to the workflow folder you want to install and use conda to create the corresponding environment:

```bash
conda env create -f environment.yml
```

If you want to use MODELLER, for the Protein preparation workflow export the KEY_MODELLER variable before creating the environment (no commercial license, avoid for production!):

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY" 
```

To run the workflows from the **BioBB Workflows Plugin for Horus**, copy the 'workflow.py' script to the /bin/ folder of your conda environment with the 'biobb_workflow' name:

```bash
cp workflow.py $CONDA_PREFIX/bin/biobb_workflow
chmod +755 $CONDA_PREFIX/bin/biobb_workflow
```

