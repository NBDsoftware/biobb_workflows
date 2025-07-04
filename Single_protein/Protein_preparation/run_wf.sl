#!/bin/bash
#SBATCH --job-name=md
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16                
#SBATCH --mem=2000
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition=standard-gpu
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Load GROMACS module
ml GROMACS

# Activate conda environment
module load Miniconda3
source activate /shared/work/BiobbWorkflows/envs/biobb_pp

# Unset env variables that might be set in the loaded modules
unset PYTHONPATH
export PATH=/shared/work/BiobbWorkflows/envs/biobb_md/bin:$PATH

# Path to the workflow
REPO_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows
WF_PATH=$REPO_PATH/Single_protein/MD_gromacs/

# Input files
INPUT_PDB=/path/to/pdb/file.pdb
OUTPUT_PATH=output

# Launch workflow
python $WF_PATH/workflow.py --config input.yml --input_pdb $INPUT_PDB --output $OUTPUT_PATH
