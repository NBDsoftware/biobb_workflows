#!/bin/bash
#SBATCH --job-name=wf1_P38
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
#SBATCH --ntasks=4
#SBATCH --time=00:15:00

# Purge loaded modules
module purge 

# Run miniconda
. /storage/apps/MINICONDA/3/etc/profile.d/conda.sh

# Activate previously created conda environment from environment.yml
conda activate eucanshare_wf1

# Load corresponding GROMACS module
module load gromacs/2019.6

# Launch workflow
python biobb_md_setup_mutation.py -i input/P38alpha3HEC.pdb -o output --op free --config input.yml

