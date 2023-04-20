#!/bin/bash
#SBATCH --job-name=single_protein_wf2
#SBATCH --ntasks=1                                      # total number of tasks across all nodes
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3/4.9.2

# Load GROMACS module
module load GROMACS/2022.3-intel-2021b

# Load SLURM
module load slurm/slurm/21.08.6 

# Activate previously created conda environment from environment.yml
source activate /home/pnavarro/.conda/envs/single_protein_wf2

# Launch workflow
python biobb_clustering_cavity_analysis.py --config input_HPC.yml 