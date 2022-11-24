#!/bin/bash
#SBATCH --job-name=wf1_P38
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
#SBATCH --ntasks=4
#SBATCH --time=00:15:00
#SBATCH --constraint=intel
#SBATCH --nodes=1

# Purge loaded modules
module purge

# Load conda
module load Miniconda3/4.9.2

# Load GROMACS module
module load GROMACS/2022.3-intel-2021b

# Load environments
source activate /shared/work/conda_envs/tmap

# Activate previously created conda environment from environment.yml
conda activate eucanshare_wf1

# Launch workflow
python biobb_md_setup_mutation.py --input input/P38alpha4LOO.pdb --config input_HPC.yml --n_trajs 1
