#!/bin/bash

####################
# Input parameters #
####################

INPUT_PDB="/path/to/input.pdb"
OUTPUT_DIR="/path/to/output"
SYSTEM_NAME="name"                # System name 
TOTAL_NUM_STEPS=20000             # Number of steps for the MD simulation - unit is defined by the timestep
STEPS_PER_PART=10000              # Number of steps per part - to respect the HPC queue time limit
SLURM_SCRIPT="run_HPC_HIS.sl"   # SLURM script to be used for the HPC job
PDB_CHAIN="A, B"                  # Chain to be used from the input PDB file, e.g. "A" or "A, B"

# Review as well the input_HPC.yml to check equilibration time for example! - this will be improved in the future

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

# Launch first part of the MD simulation (no dependency)
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-1 $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR 1 $STEPS_PER_PART 0 | awk '{print $4}')
echo "Job 1 launched with ID: $JOB_ID"
PREVIOUS_JOB_ID=$JOB_ID

# Launch the rest of the parts of the MD simulation
for i in $(seq 2 $NUM_PARTS)
do
    ACCUMULATED_STEPS=$((STEPS_PER_PART * i))
    JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$i --dependency=afterany:$PREVIOUS_JOB_ID $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR $i $ACCUMULATED_STEPS 0 | awk '{print $4}')
    echo "Job $i launched with ID: $JOB_ID"
    PREVIOUS_JOB_ID=$JOB_ID
done

# Launch the final analysis
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$(($NUM_PARTS+1)) --dependency=afterany:$PREVIOUS_JOB_ID $SLURM_SCRIPT $INPUT_PDB "$PDB_CHAIN" $OUTPUT_DIR $NUM_PARTS $TOTAL_NUM_STEPS 1 | awk '{print $4}')
echo "Job $((NUM_PARTS+1)) (analysis) launched with ID: $JOB_ID"