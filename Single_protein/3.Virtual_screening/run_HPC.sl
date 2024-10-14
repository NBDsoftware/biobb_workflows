#!/bin/bash
#SBATCH --job-name=vs_autodock
#SBATCH --nodes=1
#SBATCH --ntasks=1                                     # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                              # number of processes
#SBATCH --time=01:00:00
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin                              # send email when job begins
# #SBATCH --mail-type=end                                # send email when job ends
# #SBATCH --mail-user=pablo.navarro@nostrumbiodiscovery.com

# Purge loaded modules
module purge

# Load conda
module load Miniconda3

# Activate conda environment
source activate /shared/work/BiobbWorkflows/envs/biobb_sp_virtual_screening

# Path to the workflow
REPO_PATH=/path/to/repo/biobb_workflows
WF_PATH=$REPO_PATH/euCanShare/Single_protein_wfs/3.Virtual_screening

# Input files
INPUT_PATH=/path/to/input/folder
STRUCTURE_PATH=$INPUT_PATH/path/to/structure.pdb
INPUT_POCKETS_ZIP=$INPUT_PATH/path/to/pockets.zip  # Alternatively, provide the numbers of the residues forming the pocket

# Launch workflow 
python  $WF_PATH/biobb_docking_htvs.py --config input_HPC.yml --ligand_lib input/ligand_lib.smi --structure_path $STRUCTURE_PATH --input_pockets_zip $INPUT_POCKETS_ZIP --num_top_ligands 5