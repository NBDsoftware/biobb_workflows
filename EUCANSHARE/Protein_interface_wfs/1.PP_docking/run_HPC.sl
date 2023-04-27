#!/bin/bash
#SBATCH --job-name=PP_docking_wf
#SBATCH --nodes=1
#SBATCH --ntasks=8                         
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Load conda
module load Miniconda3/4.9.2

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/biobb_pydock

PYTHON_INTERPRETER_PATH='/home/pnavarro/.conda/envs/biobb_pydock/bin'

# Launch workflow
python biobb_pp_docking.py --config input_HPC.yml 
