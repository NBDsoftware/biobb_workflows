#!/bin/bash
#SBATCH --job-name=SuMD_test
#SBATCH --ntasks=4                                      # total number of tasks across all nodes
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --constraint=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Load conda
module load Miniconda3/4.9.2

# Load GROMACS module
module load GROMACS/2022.3-intel-2021b-CUDA.11.6.0 

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
# source activate /home/pnavarro/.conda/envs/biobb_sumd

# Launch workflow

INPUT_FOLDER=./input
STRUCTURE=$INPUT_FOLDER/npt.gro
TOPOLOGY=$INPUT_FOLDER/topology.zip
INDEX=$INPUT_FOLDER/index.ndx

/home/pnavarro/.conda/envs/biobb_sumd/bin/python biobb_SuMD.py -structure $STRUCTURE -topology $TOPOLOGY -index $INDEX -config input_HPC_GPU.yml 


