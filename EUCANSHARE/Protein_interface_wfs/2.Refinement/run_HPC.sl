#!/bin/bash
#SBATCH --job-name=pose_refinement_wf
#SBATCH --nodes=1
#SBATCH --ntasks=4                         
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Load conda
module load Miniconda3/4.9.2

# Load SLURM
module load slurm/slurm/21.08.6 

# Load GROMACS
module load GROMACS/2022.3-intel-2021b

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/single_protein_wf2

PYTHON_INTERPRETER_PATH='/home/pnavarro/.conda/envs/single_protein_wf2/bin'

# Launch workflow
python biobb_pose_refinement.py --config input_HPC.yml 