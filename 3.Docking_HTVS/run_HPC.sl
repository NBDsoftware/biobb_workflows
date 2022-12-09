#!/bin/bash
#SBATCH --job-name=wf3_mp
#SBATCH --nodes=1
#SBATCH --ntasks=1                                     # total number of tasks across all nodes
#SBATCH --cpus-per-task=8                              # number of processes
#SBATCH --time=00:30:00
#SBATCH --constraint=x86_64,intel,skl
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
#SBATCH --mail-type=begin                              # send email when job begins
#SBATCH --mail-type=end                                # send email when job ends
#SBATCH --mail-user=pablo.navarro@nostrumbiodiscovery.com

# Purge loaded modules
# module purge

# Load conda
module load Miniconda3/4.9.2

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/eucanshare_wf3

srun python biobb_docking_htvs.py --config input.yml --lig-lib input/ligand_lib_large.txt 
