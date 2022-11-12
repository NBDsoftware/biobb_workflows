# EUCANSHARE workflows

Biobb workflows created for EUCANSHARE. They can be run together or as independent workflows.

Description:

**1.MD_setup_mutation**: prepares and launches a GROMACS MD simulation starting from a PDB structure file or PDB code. Optionally, mutations can be added to the original PDB during the preparation.

**2.Clustering_CavityAnalysis**: clustering of MD trajectory and cavity analysis using Fpocket for the most representative conformations.

**3.Docking_HTVS**: high-throughput virtual screening of selected pocket using library of ligands.

## Installation 

Requirements: git, conda

1. Clone or download this repository:

```bash
git clone https://github.com/PabloNA97/eucanshare_wfs.git
```

2. Go to one of the 3 workflow folders 

3. Use conda to create the corresponding environment 

```bash
conda env create -f environment.yml
```
4. Provide necessary inputs modifying the input.yml or the command line call in run.sh and launch:

```bash
./run.sh
```
