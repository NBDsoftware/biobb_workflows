#!/bin/bash
#SBATCH --job-name=single_protein_wf2
#SBATCH --ntasks=1                                      # total number of tasks across all nodes
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3

# Load GROMACS module
module load GROMACS

# Activate previously created conda environment from environment.yml
source activate /shared/work/BiobbWorkflows/envs/biobb_sp_cavity_analysis

# Path to the workflow
REPO_PATH=/path/to/repo/biobb_workflows
WF_PATH=$REPO_PATH/euCanShare/Single_protein_wfs/2.Clustering_and_cavity_analysis

# Input files
INPUT=/path/to/input/folder
OUTPUT_PATH=/path/to/output/folder

# Launch workflow
python $WF_PATH/biobb_clustering_cavity_analysis.py --config input_HPC.yml --clustering_path $INPUT --output $OUTPUT_PATH
