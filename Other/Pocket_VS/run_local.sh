#!/bin/bash

# MODELS_PATH=/home/nbd-pablo/Documents/2022_ENSEM/3.P53/MD_Trajectories/3ts8/CLUSTER/selection2/low_num_clusters/MDTRAJ_HIERA/clustering  
MODELS_PATH=/home/nbd-pablo/Documents/2023_EUCANSHARE/Javi_pdbs

LIGAND_LIB_PATH=/home/nbd-pablo/repos/biobb_workflows/OTHER/4.Pocket_VS/input/UCI_subset_reduced.smi

POCKET_RESIDUES='19 20 21 22'

python biobb_pocket_vs.py --main-config input.yml --cavity-config cavity_input.yml --htvs-config htvs_input.yml --lig-lib $LIGAND_LIB_PATH --input $MODELS_PATH --pocket-res "$POCKET_RESIDUES" --until filtering 

# python biobb_pocket_vs.py --main-config input.yml --htvs-config htvs_input.yml --lig-lib $LIGAND_LIB_PATH --input $MODELS_PATH --until all

