#!/bin/bash
 
CLUSTERING_PATH=/home/nbd-pablo/repos/biobb_workflows/euCanShare/Single_protein_wfs/2.Clustering_and_cavity_analysis/input

python biobb_clustering_cavity_analysis.py --clustering_path $CLUSTERING_PATH --config input.yml
