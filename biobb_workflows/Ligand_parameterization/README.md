# Ligand Parameterization Workflow

This workflow is designed to parameterize ligands for molecular simulations. It automates the process of generating force field parameters for small molecules, which can then be used in various molecular dynamics (MD) simulations.

## Introduction

Ligand parameterization is a crucial step in preparing molecular simulations. This workflow generates parameters for Amber or Gromacs simulations. The protonation can be done with obabel, ambertools or skipped (if the ligand is already protonated). The workflow accepts .frcmod and .prep files from the [Amber Parameter Database from Manchester University](http://amber.manchester.ac.uk/). The output will be a folder containing .gro and .itp files for Gromacs or .frcmod and .prep/.lib files for Amber.

Note that the parameterization for Gromacs will involve producing the final .gro file as well as the .itp file. This is because the .gro and .itp files have to agree between each other. They are merged with the main protein coordinates and topology file after using pdb2gmx. Thus the parameterization workflow will have to be run once for every new PDB, as the coordinates in the .gro file will change for each system. (NOTE: one possible workaround would be to define each ligand as a residue template .rtp, so that pdb2gmx could re-construct any missing atoms when reading a new PDB.). In the case of Amber however, tleap is able to reconstruct the missing atoms in the PDB by using the .prep/.lib files. So there is no need to run the parameterization workflow for each new PDB if one already has the corresponding .frcmod and .prep/.lib files. 

