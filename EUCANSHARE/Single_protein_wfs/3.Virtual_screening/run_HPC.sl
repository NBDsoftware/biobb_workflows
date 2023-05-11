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

INPUT=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Single_protein_wfs/3.Virtual_screening/input
STRUCTURE_PATH=$INPUT/9/step2_extract_models/cluster.pdb
INPUT_POCKETS_ZIP=$INPUT/9/step4_filter_cavities/filtered_pockets.zip

# python biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib.smi --structure_path $STRUCTURE_PATH --input_pockets_zip $INPUT_POCKETS_ZIP

python biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib.smi --structure_path $STRUCTURE_PATH --input_pockets_zip $INPUT_POCKETS_ZIP --num_top_ligands 5