#!/bin/bash
#SBATCH --job-name=single_protein_wf2
#SBATCH --ntasks=1                                      # total number of tasks across all nodes
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3/4.9.2

# Load GROMACS module
module load GROMACS/2022.3-intel-2021b-CUDA.11.6.0

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/single_protein_wf2

# Test 1
# INPUT=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/input
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/output1
# TRAJ_PATH=$INPUT/all_trajectories.xtc
# TOP_PATH=$INPUT/dry_structure.gro
# python biobb_clustering_cavity_analysis.py --config input_HPC.yml --traj_path $TRAJ_PATH --top_path $TOP_PATH --output $OUTPUT_PATH

# Test 2
# INPUT=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/input2
# CLUSTERS=$INPUT/clusters
# OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/output2
# python biobb_clustering_cavity_analysis.py --config input_HPC.yml --clustering_path $CLUSTERS --output $OUTPUT_PATH

# Test 3
INPUT=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/input2
CLUSTERS=$INPUT/clusters
OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/output3
OUTPUT_SUMMARY_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/2.Clustering_and_cavity_analysis/summary3.yml
python biobb_clustering_cavity_analysis.py --config input_HPC.yml --clustering_path $CLUSTERS --output $OUTPUT_PATH --output_summary $OUTPUT_SUMMARY_PATH