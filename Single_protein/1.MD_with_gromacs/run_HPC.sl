#!/bin/bash
#SBATCH --job-name=single_protein_wf1
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=8                                  # total number of tasks across all nodes
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
# #SBATCH --constraint=gpu
# #SBATCH --gres=gpu:1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3

# Load slurm to call srun in the computation node
module load slurm

# Load GROMACS module
module load GROMACS

# Activate previously created conda environment from environment.yml
source activate /shared/work/BiobbWorkflows/envs/biobb_md

# Path to the workflow
REPO_PATH=/path/to/repo/biobb_workflows
WF_PATH=$REPO_PATH/euCanShare/Single_protein_wfs/1.MD_with_mutation

# Input files
INPUT_PDB=/path/to/input/folder/structure.pdb
OUTPUT_PATH=/path/to/output/folder

# Launch workflow
python $WF_PATH/biobb_md_setup_mutation.py --config input_HPC.yml --num_parts 1 --input_pdb $INPUT_PDB --output $OUTPUT_PATH --pdb_chains A 
