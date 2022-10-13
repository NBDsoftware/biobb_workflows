# EUCANSHARE workflows

Biobb workflows created for EUCANSHARE. They can be run together or as independent workflows.

Description:

**1.MD_setup_mutation**: prepares and launches a GROMACS MD simulation starting from a PDB structure file or PDB code. Optionally, mutations can be added to the original PDB during the preparation.

**2.Clustering_CavityAnalysis**: clustering of MD trajectory and cavity analysis using Fpocket for the most representative conformations.

**3.Docking_HTVS**: high-throughput virtual screening of selected pocket using library of ligands.

## Installation 

Requirements: git, conda

1. Clone or download this repository

```bash
git clone 
```

2. Go to workflow folder

3. Use conda to create the environment 

```bash
conda env create -f environment.yml
```
4. Provide necessary input and launch, see run.sh 
