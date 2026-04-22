# BioBB Workflows

BioBB workflows are ready-to-use pipelines built using BioExcel Building Blocks (BioBB) to perform common tasks in biomolecular simulations and modeling. Each workflow is a CLI entrypoint of the `biobb_workflows` package.

Different workflows can be combined to create more complex pipelines.

## Single protein pipeline

![alt text](../img/Single_protein_scheme.png?raw=true)

## Protein-protein docking pipeline

![alt text](../img/protein_protein_scheme.png?raw=true)

## Workflows

**protein_preparation**: prepares a protein structure from a PDB file for further simulations or modeling. Includes adding missing residues/atoms and selecting the protonation state of residues at a given pH.

**ligand_parameterization**: parameterizes a small molecule using Antechamber and the GAFF force field. Outputs GROMACS or AMBER topology files.

**md_gromacs**: prepares and launches a GROMACS MD simulation starting from a prepared PDB structure file.

**traj_postprocessing**: post-processes a GROMACS MD trajectory: strips solvent, centers, images and fits.

**cavity_analysis**: clusters an MD trajectory and runs a cavity analysis using Fpocket on the representative structures.

**vs_autodock**: high-throughput virtual screening of a selected pocket using a ligand library (SMILES/SDF) and AutoDock Vina.

**pp_docking**: samples and scores protein-protein poses with pyDock, then clusters by RMSD and returns the best pose per cluster.

**pose_refinement**: takes a zip of protein-protein docking poses and minimizes the energy of each with GROMACS.

## Installation

Requirements: `git`, `conda`

```bash
git clone https://github.com/NBDsoftware/biobb_workflows.git
cd biobb_workflows
```

> To use MODELLER in `protein_preparation`, export your key before creating the environment:
> ```bash
> export KEY_MODELLER="YOUR_MODELLER_KEY"
> ```

### MD environment
Covers: `protein_preparation`, `ligand_parameterization`, `md_gromacs`, `traj_postprocessing`

```bash
conda env create -f environment_md.yml
conda activate biobb_md
pip install -e .[gromacs,ambertools]
```

### VS environment
Covers: `vs_autodock`, `cavity_analysis`

```bash
conda env create -f environment_vs.yml
conda activate biobb_vs
pip install -e .[gromacs,vina]
```

## Usage

Once installed, each workflow is available as a CLI command:

| Command | Workflow |
|---------|----------|
| `biobb_md` | md_gromacs |
| `biobb_lp` | ligand_parameterization |
| `biobb_pp` | protein_preparation |
| `biobb_traj` | traj_postprocessing |
| `biobb_ca` | cavity_analysis |
| `biobb_vs` | vs_autodock |

```bash
biobb_md --help
biobb_pp --input_pdb protein.pdb --output output/
```
