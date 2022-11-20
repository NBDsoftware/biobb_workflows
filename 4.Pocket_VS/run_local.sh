#!/bin/bash

MODELS_PATH=/home/nbd-pablo/Documents/2022_ENSEM/3.P53/MD_Trajectories/3ts8/CLUSTER/selection2/low_num_clusters/MDTRAJ_HIERA/clustering  
LIGAND_LIB_PATH=/home/nbd-pablo/repos/eucanshare_wfs/4.Pocket_VS/ligand_lib.txt

python biobb_pocket_vs.py --main-config input.yml --htvs-config htvs_input.yml --lig-lib $LIGAND_LIB_PATH --input $MODELS_PATH --until all
