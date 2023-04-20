# Single protein workflows

Single protein biobb workflows created for EUCANSHARE. They can be run together or as independent workflows. They can be run entirely or just partially to analyze intermediate results before proceeding.

## Description

**1.MD_with_mutation**: prepares and launches a GROMACS MD simulation starting from a PDB structure file or PDB code. Optionally, several mutations can be added to the original PDB during the preparation. See corresponding folder for more details.

**2.Clustering_and_cavity_analysis**: clustering of MD trajectory and extraction of most representative structures. Followed by a cavity analysis using Fpocket on those structures. See corresponding folder for more details.

**3.Docking_HTVS**: high-throughput virtual screening of selected pocket using library of ligands and Autodock. See corresponding folder for more details.

## General installation steps

Requirements: git, conda

1. Clone or download this repository:

```bash
git clone https://github.com/PabloNA97/biobb_workflows.git
```

2. Go to one of the workflow folders 

3. Use conda to create the corresponding environment (for workflow 1 export KEY_MODELLER variable as well):

```bash
export KEY_MODELLER=MODELIRANJE
```

```bash
conda env create -f environment.yml
```

4. Check the environment is listed and activate it

```bash
conda env list
conda activate <name_of_environment>
```

5. Provide necessary inputs and launch (look for run_local.sh or run_HPC.sl):

```bash
./run.sh
```

```bash
sbatch run_HPC.sl
```

