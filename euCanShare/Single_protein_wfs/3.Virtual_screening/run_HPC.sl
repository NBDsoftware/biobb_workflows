#!/bin/bash
#SBATCH --job-name=vs_autodock
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
module purge

# Load conda
module load Miniconda3/4.9.2

# Activate conda environment
source activate /shared/work/BiobbWorkflows/envs/biobb_sp_virtual_screening

# Path to the workflow python script
WF_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows/euCanShare/Single_protein_wfs/3.Virtual_screening

# Input files
INPUT_PATH=/shared/scratch/jobs/pnavarro/2023_EUCANSHARE/biobb_workflows/euCanShare/Single_protein_wfs/3.Virtual_screening/input
STRUCTURE_PATH=$INPUT_PATH/9/step2_extract_models/cluster.pdb
INPUT_POCKETS_ZIP=$INPUT_PATH/9/step4_filter_cavities/filtered_pockets.zip

# Launch workflow 
python  $WF_PATH/biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib.smi --structure_path $STRUCTURE_PATH --input_pockets_zip $INPUT_POCKETS_ZIP --num_top_ligands 5