#!/bin/bash

####################
# Input parameters #
####################

INPUT_PDB="/path/to/input.pdb"
OUTPUT_DIR="/path/to/output"
LIGAND_PARAMETERS="/path/to/ligand_parameters/folder"
SYSTEM_NAME="name"              # System name 
TOTAL_NUM_STEPS=20000           # Number of steps for the MD simulation - unit is defined by the timestep
STEPS_PER_PART=10000            # Number of steps per part - to respect the HPC queue time limit
SLURM_SCRIPT="run_MD_bck.sl"   # SLURM script to be used for the HPC job
PDB_CHAIN="A"                   # Chain to be used from the input PDB file, e.g. "A" or "A, B"

# Review as well the input.yml to check equilibration time for example! - this will be improved in the future

###############
# MAIN SCRIPT #
###############

# Check TOTAL_NUM_STEPS is greater than STEPS_PER_PART
if [ $TOTAL_NUM_STEPS -lt $STEPS_PER_PART ]; then
    echo "ERROR: TOTAL_NUM_STEPS ($TOTAL_NUM_STEPS) must be greater than STEPS_PER_PART ($STEPS_PER_PART)"
    exit 1
fi

# Check if the input PDB file exists
if [ ! -f $INPUT_PDB ]; then
    echo "ERROR: Input PDB file not found ($INPUT_PDB)"
    exit 1
fi

# Find number of parts
NUM_PARTS=$((TOTAL_NUM_STEPS / STEPS_PER_PART))
echo "Number of parts: $NUM_PARTS"

# Control final analysis execution
DO_FINAL_ANALYSIS=false

# Launch part 1 of the MD simulation (no dependency)
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-1 $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR $LIGAND_PARAMETERS 1 $STEPS_PER_PART $DO_FINAL_ANALYSIS | awk '{print $4}')
echo "Job 1 launched with ID: $JOB_ID"
PREVIOUS_JOB_ID=$JOB_ID

# Launch all other parts of the MD simulation
for PART in $(seq 2 $NUM_PARTS)
do
    ACCUMULATED_STEPS=$((STEPS_PER_PART * PART))
    JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$PART --dependency=afterany:$PREVIOUS_JOB_ID $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR $LIGAND_PARAMETERS $PART $ACCUMULATED_STEPS $DO_FINAL_ANALYSIS | awk '{print $4}')
    echo "Job $PART launched with ID: $JOB_ID"
    PREVIOUS_JOB_ID=$JOB_ID
done

# Launch the final analysis
DO_FINAL_ANALYSIS=true
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$(($NUM_PARTS+1)) --dependency=afterany:$PREVIOUS_JOB_ID $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR $LIGAND_PARAMETERS $NUM_PARTS $TOTAL_NUM_STEPS $DO_FINAL_ANALYSIS | awk '{print $4}')
echo "Job $((NUM_PARTS+1)) (analysis) launched with ID: $JOB_ID"