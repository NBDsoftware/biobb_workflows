#!/bin/bash
#SBATCH --job-name=ligand_param   # Job name
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=8                # total number of tasks across all nodes
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
# #SBATCH --constraint=gpu
# #SBATCH --gres=gpu:1
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# #SBATCH --mail-type=begin        # send email when job begins
# #SBATCH --mail-type=end          # send email when job ends
# #SBATCH --mail-user=your@email.com

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3

# Activate previously created conda environment from environment.yml
source activate /shared/work/BiobbWorkflows/envs/biobb_md

# Path to the workflow
REPO_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows
WF_PATH=$REPO_PATH/Single_protein/Ligand_parameterization

# Input files
INPUT_PDB=path/to/input/pdb/file.pdb
OUTPUT_PATH=output
CUSTOM_PARAMETERS=path/to/custom_parameters/folder

# Launch workflow
python $WF_PATH/workflow.py --config input.yml --input_pdb $INPUT_PDB --ligands LIG HEM --format gromacs --custom_parameters $CUSTOM_PARAMETERS --output $OUTPUT_PATH
