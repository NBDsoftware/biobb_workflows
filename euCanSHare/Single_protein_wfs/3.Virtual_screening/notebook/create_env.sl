#!/bin/bash
#SBATCH --job-name=env_creation   
#SBATCH --ntasks=4                        # Number of tasks 
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --time=04:00:00             
#SBATCH --mem-per-cpu=3000       # 1 Gb of RAM per cpu
#SBATCH --output=report_%j.out   # Name of file with standard output
#SBATCH --error=report_%j.err    # Name of file with standard input

# Purge loaded modules - erase all previously loaded modules
module purge

# Load conda (Miniconda is the lightweight version)
ml Miniconda3

conda env create -f environment.yml