#!/bin/bash

####################
# Input parameters #
####################

SYSTEM_NAME=1r9o                            # System name
INPUT_PDB=1r9o.pdb                          # Path to the input PDB file
PDB_CHAINS="A"                              # Chains to be simulated from the input PDB file. Ex: "A B"

LIG_SLURM_SCRIPT=run_LP.sl  
MD_SLURM_SCRIPT=run_MD.sl   

# Ligand parameterization
LIGANDS="HEM FLP"                           # Ligands to be parameterized - the rest are discarded. Leave empty "" if no ligands are present
LIG_BASE_FF=protein.ff14SB                  # Base Force field to be used for the ligand parameterization of those ligands with parameter modifications in custom_parameters
CUSTOM_PARAMETERS=custom_parameters         # Folder with the custom parameters for one or more ligands 
PARAM_OUTPUT_PATH=parameterization_output   # Output path for the ligand parameterization workflow
PROTONATION_TOOL=ambertools                 # Tool to protonate the ligands. Options: ambertools, obabel, none

# MD simulation
MD_FF=amber99sb-ildn                        # Force field to be used for the MD simulation, see the list of available force fields in pdb2gmx                                
MD_OUTPUT_PATH=md_output                    # Output path for the MD simulation workflow
TOTAL_NUM_STEPS=20000                       # Number of steps for the MD simulation - unit is defined by the timestep
STEPS_PER_PART=10000                        # Number of steps per part - to respect the HPC queue time limit

###########################
# LIGAND PARAMETERIZATION #
###########################

if [ -n "$LIGANDS" ]; then
    # If there are any ligands to be parameterized, launch the ligand parameterization
    LIGAND_TOP_PATH=custom_topologies       # Output path for the resulting ligand topologies
    PARAM_JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-parameterization $LIG_SLURM_SCRIPT $INPUT_PDB "$PDB_CHAINS" $LIG_BASE_FF $CUSTOM_PARAMETERS $PROTONATION_TOOL "$LIGANDS" $LIGAND_TOP_PATH $PARAM_OUTPUT_PATH | awk '{print $4}')
    echo "Parameterization job launched with ID: $PARAM_JOB_ID"
    FIRST_DEPENDENCY="--dependency=afterany:$PARAM_JOB_ID"
else
    # If there are no ligands to be parameterized, set the PARAM_JOB_ID to 0
    LIGAND_TOP_PATH="None"
    FIRST_DEPENDENCY=""
fi

#################
# MD SIMULATION #
#################

# Check TOTAL_NUM_STEPS is greater than STEPS_PER_PART
if [ $TOTAL_NUM_STEPS -lt $STEPS_PER_PART ]; then
    echo "ERROR: TOTAL_NUM_STEPS ($TOTAL_NUM_STEPS) must be greater than STEPS_PER_PART ($STEPS_PER_PART)"
    exit 1
fi

# Find number of parts
NUM_PARTS=$((TOTAL_NUM_STEPS / STEPS_PER_PART))
echo "Number of parts: $NUM_PARTS"

# Control final analysis execution
DO_FINAL_ANALYSIS=false

# Launch part 1 of the MD simulation
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-1 $FIRST_DEPENDENCY $MD_SLURM_SCRIPT $INPUT_PDB "$PDB_CHAINS" $MD_FF $LIGAND_TOP_PATH 1 $STEPS_PER_PART $DO_FINAL_ANALYSIS $MD_OUTPUT_PATH | awk '{print $4}')
echo "Job 1 launched with ID: $JOB_ID"
PREVIOUS_JOB_ID=$JOB_ID

# Launch all other parts of the MD simulation
for PART in $(seq 2 $NUM_PARTS)
do
    ACCUMULATED_STEPS=$((STEPS_PER_PART * PART))
    JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$PART --dependency=afterany:$PREVIOUS_JOB_ID $MD_SLURM_SCRIPT $INPUT_PDB "$PDB_CHAINS" $MD_FF $LIGAND_TOP_PATH $PART $ACCUMULATED_STEPS $DO_FINAL_ANALYSIS $MD_OUTPUT_PATH | awk '{print $4}')
    echo "Job $PART launched with ID: $JOB_ID"
    PREVIOUS_JOB_ID=$JOB_ID
done

# Launch the final analysis
DO_FINAL_ANALYSIS=true
JOB_ID=$(sbatch -q normal --job-name=$SYSTEM_NAME-$(($NUM_PARTS+1)) --dependency=afterany:$PREVIOUS_JOB_ID $MD_SLURM_SCRIPT $INPUT_PDB "$PDB_CHAINS" $MD_FF $LIGAND_TOP_PATH $NUM_PARTS $TOTAL_NUM_STEPS $DO_FINAL_ANALYSIS $MD_OUTPUT_PATH | awk '{print $4}')
echo "Job $((NUM_PARTS+1)) (analysis) launched with ID: $JOB_ID"