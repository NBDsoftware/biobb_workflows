#!/bin/bash
#SBATCH --job-name=single_protein_wf3
#SBATCH --nodes=1
#SBATCH --ntasks=1                                     # total number of tasks across all nodes
#SBATCH --cpus-per-task=8                              # number of processes
#SBATCH --time=01:00:00
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=pablo.navarro@nostrumbiodiscovery.com

# Purge loaded modules
# module purge

# Load conda
module load Miniconda3/4.9.2

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/single_protein_wf3

srun python biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib.smi 
