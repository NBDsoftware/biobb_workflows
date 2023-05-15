#!/bin/bash
#SBATCH --job-name=protein_protein_wf4
#SBATCH --nodes=1
#SBATCH --ntasks=1                                      # total number of tasks across all nodes
#SBATCH --cpus-per-task=12                              # number of processes
#SBATCH --time=24:00:00
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

# Launch workflow
INPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Protein_interface_wfs/4.Virtual_screening/input
STRUCTURE_PATH=$INPUT_PATH/13_refined_pose.pdb
INPUT_POCKETS_ZIP=$INPUT_PATH/filtered_pockets.zip
OUTPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/EUCANSHARE/Protein_interface_wfs/4.Virtual_screening/output2
python biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib_large.smi --structure_path $STRUCTURE_PATH --input_pockets_zip $INPUT_POCKETS_ZIP --output $OUTPUT_PATH --num_top_ligands 10