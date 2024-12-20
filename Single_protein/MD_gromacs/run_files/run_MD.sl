#!/bin/bash
##### Number of tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00:00
#####Output and error log files
#SBATCH -e gmx_%j.e
#SBATCH -o gmx_%j.o

# Purge loaded modules
module purge 

# Load GROMACS module
ml GROMACS

# Activate conda environment
module load Miniconda3
source activate /shared/work/BiobbWorkflows/envs/biobb_md

# Path to the workflow
REPO_PATH=/shared/work/BiobbWorkflows/src/biobb_workflows
WF_PATH=$REPO_PATH/Single_protein/MD_gromacs/

# Input arguments
INPUT_PDB=$1
PDB_CHAINS=$2
MD_FF=$3
LIGAND_TOP_PATH=$4
NUM_PARTS=$5
NSTEPS=$6
DO_FINAL_ANALYSIS=$7
OUTPUT_PATH=$8

echo "Input PDB: $INPUT_PDB"
echo "PDB Chains: $PDB_CHAINS"
echo "MD Forcefield: $MD_FF"
echo "Ligand Top Path: $LIGAND_TOP_PATH"
echo "Number of parts: $NUM_PARTS"
echo "Number of steps: $NSTEPS"
echo "Do final analysis: $DO_FINAL_ANALYSIS"

# Read the condition to do the final analysis
if [ "$DO_FINAL_ANALYSIS" == "true" ]; then
    FINAL_ANALYSIS="--final_analysis"
else
    FINAL_ANALYSIS=""
fi

#################
# MD SIMULATION #
#################

python $WF_PATH/workflow.py --config input_md.yml --forcefield $MD_FF --num_parts $NUM_PARTS --input_pdb $INPUT_PDB --output $OUTPUT_PATH --pdb_chains $PDB_CHAINS --nsteps $NSTEPS --ligands_folder $LIGAND_TOP_PATH --his_protonation_tool "pdb4amber" --fix_ss --fix_amides $FINAL_ANALYSIS 
