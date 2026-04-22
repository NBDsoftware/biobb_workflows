#!/bin/bash
##### Number of tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1-00:00:00
#####Output and error log files
#SBATCH -e param_%j.e
#SBATCH -o param_%j.o

# Purge loaded modules
module purge 

# Load conda / miniconda module
module load Miniconda3

# Activate conda environment
source activate /shared/work/BiobbWorkflows/envs/biobb_md

# Path to the workflow
REPO_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows
WF_PATH=$REPO_PATH/Single_protein/Ligand_parameterization

# Input arguments
INPUT_PDB=$1
PDB_CHAINS=$2
LIG_BASE_FF=$3
CUSTOM_PARAMETERS=$4
PROTONATION_TOOL=$5
LIGANDS=$6
LIGAND_TOP_PATH=$7
OUTPUT_PATH=$8

echo "Input PDB: $INPUT_PDB"
echo "PDB Chains: $PDB_CHAINS"
echo "Ligand Base Forcefield: $LIG_BASE_FF"
echo "Custom Parameters path: $CUSTOM_PARAMETERS"
echo "Protonation Tool: $PROTONATION_TOOL"
echo "Ligands: $LIGANDS"
echo "Ligands Topologies Path: $LIGAND_TOP_PATH"
echo "Output Path: $OUTPUT_PATH"

# Launch workflow
python $WF_PATH/workflow.py --config input_param.yml --input_pdb $INPUT_PDB --forcefields $LIG_BASE_FF --ligands $LIGANDS --chains $PDB_CHAINS --format gromacs --custom_parameters $CUSTOM_PARAMETERS --protonation_tool $PROTONATION_TOOL --output_top_path $LIGAND_TOP_PATH --output $OUTPUT_PATH
