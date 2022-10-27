#!/bin/bash
#SBATCH --job-name=wf1_P38
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
#SBATCH --ntasks=4
#SBATCH --time=00:15:00

module purge 

# Run miniconda
. /storage/apps/MINICONDA/3/etc/profile.d/conda.sh

# Activate previously created conda environment from environment.yml
conda activate eucanshare_wf1

#load mpi libriries in MN
#ml mkl/2019.4.243   impi/2019.4.243   intel/2019.4.243

#load mpi libraries in Tirant
# module load intel mkl impi

#load tirant cluster grommacs enviroment
# /storage/apps/GROMACS/2019.6/INTEL/bin/GMXRC.bash

module load gromacs/2019.6

python biobb_md_setup_mutation.py -i input/P38alpha3HEC.pdb -o output --op free --config input.yml

