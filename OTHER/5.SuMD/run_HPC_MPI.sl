#!/bin/bash
#SBATCH --job-name=sumd
#SBATCH --nodes=1
#SBATCH --ntasks=8                                      # total number of tasks across all nodes
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
##SBATCH --mail-type=all          # send email when job begins, ends and fails
##SBATCH --mail-user=user@mail.com

# Purge loaded modules
module purge 

# Load conda
module load Miniconda3/4.9.2

# Load GROMACS module
module load GROMACS/2022.3-intel-2021b

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/biobb_sumd

#LD_LIBRARY_PATH=/home/pnavarro/.conda/envs/eucanshare_wf1/lib

#export LD_LIBRARY_PATH
#echo $LD_LIBRARY_PATH

# Launch workflow

INPUT_FOLDER=./input
STRUCTURE=$INPUT_FOLDER/complex.gro
TOPOLOGY=$INPUT_FOLDER/topology.zip

/home/pnavarro/.conda/envs/biobb_sumd/bin/python biobb_SuMD.py -input $STRUCTURE -topology $TOPOLOGY -config input_HPC.yml 

