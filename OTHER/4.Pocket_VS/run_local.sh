#!/bin/bash

# MODELS_PATH=/home/nbd-pablo/Documents/2022_ENSEM/3.P53/MD_Trajectories/3ts8/CLUSTER/selection2/low_num_clusters/MDTRAJ_HIERA/clustering  
MODELS_PATH=/home/nbd-pablo/Documents/2022_ENSEM/3.P53/MD_Trajectories/pocket_vs_results/best_affinities_from_trajs/clustering/cpptraj_clust/10_clusters/representative

LIGAND_LIB_PATH=/home/nbd-pablo/repos/eucanshare_wfs/4.Pocket_VS/input/UCI_subset_active.smi
POCKET_RESIDUES='19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 47 48 49 50 51'

python biobb_pocket_vs.py --main-config input.yml --htvs-config htvs_input.yml --lig-lib $LIGAND_LIB_PATH --input $MODELS_PATH --pocket-res "$POCKET_RESIDUES" --until all

# python biobb_pocket_vs.py --main-config input.yml --htvs-config htvs_input.yml --lig-lib $LIGAND_LIB_PATH --input $MODELS_PATH --until all

