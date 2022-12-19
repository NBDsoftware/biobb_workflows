#!/bin/bash

TRAJ_PATH=/home/nbd-pablo/Documents/2022_EUCANSHARE/NBD_cluster_results/1.MD_setup_mutation_mpirun/output/step24_trjcat/all_trajectories.trr
TOP_PATH=/home/nbd-pablo/Documents/2022_EUCANSHARE/NBD_cluster_results/1.MD_setup_mutation_mpirun/output/step23_dry/imaged_structure.gro 
CLUSTERING_PATH=/home/nbd-pablo/Documents/2022_ENSEM/3.P53/MD_Trajectories/3ts8/CLUSTER/selection2/low_num_clusters/HIERAGGLO/reps

python biobb_clustering_cavity_analysis.py --input-traj $TRAJ_PATH --input-top $TOP_PATH --config input.yml

# python biobb_clustering_cavity_analysis.py --input-clust $CLUSTERING_PATH --config input.yml
