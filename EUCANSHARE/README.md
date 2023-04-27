# EUCANSHARE workflows

Biobb workflows created for EUCANSHARE. They can be run together or as independent workflows. They can be run entirely or just partially to analyze intermediate results before proceeding.

## Description

**1.MD_setup_mutation**: prepares and launches a GROMACS MD simulation starting from a PDB structure file or PDB code. Optionally, mutations can be added to the original PDB during the preparation. See corresponding folder for more details.

**2.Clustering_CavityAnalysis**: clustering of MD trajectory and cavity analysis using Fpocket for the most representative conformations. See corresponding folder for more details.

**3.Docking_HTVS**: high-throughput virtual screening of selected pocket using library of ligands. See corresponding folder for more details.

## General installation steps

Requirements: git, conda

1. Clone or download this repository:

```bash
git clone https://github.com/PabloNA97/eucanshare_wfs.git
```

2. Go to one of the 3 workflow folders

3. Use conda to create the corresponding environment (for workflow 1 export KEY_MODELLER variable as well):

```bash
export KEY_MODELLER=HERE YOUR MODELLER KEY
```

```bash
conda env create -f environment.yml
```

4. Check environment is listed and activate environment

```bash
conda env list
conda activate eucanshare_wf<number of workflow>
```

5. Provide necessary inputs modifying the input.yml or the command line call in run.sh and launch:

```bash
./run.sh
```
